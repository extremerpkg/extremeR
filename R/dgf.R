#' Simulate data generating process
#'
#' Function to simulate data from a relevant data generating process (DGP) ;
#' currently this function supports the creation of DGPs with 1 layer of area effects (i.e. small-area effects)
#'
#' @param n.sims how many samples of simulated data would you like?;
#' @param n how large (sample size) should each sample be?;
#' @param pi.hat.naive what should be the fraction of cases in the sample?;
#' @param p how many normally-distributed covariates should the DGP have?;
#' @param X_corr what should be the average correlation among these covariates?;
#' @param pi what is the population-level probability of sampling a case?
#' @param Moran.I.corr what degree of global spatial autocorrelation (Moran I) should the underlying DGP have?;
#' @param spatial_structure on which map should the data be simulated ? (scotland_lipcancer, pennsylvania_lungcancer, and newyork_leukemia)
#'
#' @importFrom R2jags jags
#' @importFrom MASS mvrnorm
#' @importFrom stats runif rnorm optimize
#' @importFrom ape Moran.I
#' @return a list object
#'
#' @export
#'
#' @examples
dgf = function(
    n.sims = 1,
    n = 100,
    pi.hat.naive = 0.5,
    p = 1,
    X_corr = 0,
    pi = 0.05,
    Moran.I.corr = 0.8,
    spatial_structure="scotland_lipcancer"
){

  # some errors
  if(p<1){stop('number of regression coefficients excluding the intercept, `p`, must be at least 1')}
  if(any(X_corr < 0 | X_corr > 1)) stop('correlation between covariates, `X_corr`, must be between 0 and 1')
  if(!spatial_structure %in% c("scotland_lipcancer","pennsylvania_lungcancer","newyork_leukemia")){stop('we do not have the map you have specified\npelase choose between scotland_lipcancer, pennsylvania_lungcancer and newyork_leukemia')}


  # define some observable quantities to calculate for the process to be generated
  n1 = n*pi.hat.naive
  n0 = n-n1
  P1 = (n1+pi*n0)/pi
  P0 = n0
  theta_0 = 0
  theta_1 = n1/(n1+pi*n0)

  # simulate from a CAR pocess on an existing map

  sf::sf_use_s2(FALSE)
  if(spatial_structure=="scotland_lipcancer"){
    temp.sp= scotland$spatial.polygon
    temp.sf= sf::st_as_sf(temp.sp)
    temp.sf$small_area = names(scotland$spatial.polygon)
  }
  if(spatial_structure=="pennsylvania_lungcancer"){
    temp.sp= pennLC$spatial.polygon
    temp.sf= sf::st_as_sf(temp.sp)
    temp.sf$small_area = names(pennLC$spatial.polygon)
    temp.sf$y = as.data.table(pennLC$data)[,lapply(.SD,sum),by = "county",.SDcols = c("cases")]$cases
  }
  if(spatial_structure=="newyork_leukemia"){
    temp.sp = NYleukemia$spatial.polygon
    temp.sf= sf::st_as_sf(temp.sp)
    temp.sf$small_area = names(NYleukemia$spatial.polygon)
    temp.sf$y = NYleukemia$data$cases
  }
  # get distance matrix
  lat_lon = st_coordinates(st_centroid(temp.sf))
  colnames(lat_lon) = c("lat","lon")
  D = geodist(lat_lon,measure = "geodesic")
  # ensure graph is fully connexted
  nb = addnbs(sp.sample = temp.sf,ID = temp.sf$small_area,D=D);
  # Have we achieved a fully connected graph ?
  if(isDisconnected(nb)==FALSE){print("Success! The Graph is Fully Connected")}else{print("Failure... Some parts of the graph are still disconnected...")}
  # get neighborhood objects - these will be useful later for the fitting of the model
  nbs = nb2graph(nb);

  nb_objects = list(
    N_small_area = nbs$N,
    node1_small_area = nbs$node1,
    node2_small_area = nbs$node2,
    N_small_area_edges = nbs$N_edges,
    scaling_factor = scale_nb_components(poly2nb(temp.sf,queen = TRUE))[1],
    nb  = nb2mat(neighbours = nb ,style = 'B',zero.policy = TRUE),
    sp.object = temp.sp
  )
  # # # simulate conditionally autoregressive effects
  jags.model = "model{
  phi ~ dmnorm(zero, tau * (D - alpha*W));
  tau ~ dgamma(1,1);
  }"
  data = list(W = nb_objects$nb,
              D = diag(rowSums(nb_objects$nb)),
              zero = rep(0,nb_objects$N_small_area),
              alpha = 0.99999)
  CAR.sample =
    jags(data=data,
         parameters.to.save=c("phi",'tau'),
         model.file=textConnection(jags.model),
         n.iter=100,
         DIC = F)
  # # # pick simulated data closest to the wanted Moran I
  sample.Moran.I = apply(CAR.sample$BUGSoutput$sims.list$phi,1,function(x){ape::Moran.I(x,nb_objects$nb)}$observed)
  temp.sf$y = as.numeric(scale(CAR.sample$BUGSoutput$sims.list$phi[which(abs(sample.Moran.I-Moran.I.corr)==min(abs(sample.Moran.I-Moran.I.corr))),]))

  # scale the effects for cross-simulation consistency
  gamma = as.numeric(scale(temp.sf$y))
  small_area_id = as.integer(as.factor(round(runif(n = n,min = 1,max = nb_objects$N_small_area))))

  # draw covariates and covariate coefficients
  Sigma = array(1,c(p,p))*X_corr
  diag(Sigma) = 1
  X = cbind(1,mvrnorm(n = n,mu = rep(0,p) ,Sigma = Sigma))
  beta = sapply(X = 1:p,function(x){rnorm(n = 1,mean = 0,sd = 1)})
  # calculate the total number of regression coefficients, including the intercept
  p = dim(X)[2]

  # calibrate intercept of linear DGP to match pi.hat.naive
  intercept.optimizer =
    optimize(f = function(beta0_true){

      beta_temp = c(beta0_true,beta)
      mu = log(P1/P0) + X %*% beta_temp +  gamma[small_area_id]
      rho = inv_logit(mu)

      ybar_gen = (1-rho)*theta_0 + rho*theta_1
      temp = pi.hat.naive - mean(ybar_gen)

      return(abs(temp))
    },
    interval = c(-30,30),
    tol =  1/10000000000)

  # after calibration, generate sample
  beta_temp = c(intercept.optimizer $minimum,beta)
  mu = log(P1/P0) + X %*% beta_temp +  gamma[small_area_id]
  rho = inv_logit(mu)

  y_gen = sapply(1:n.sims,function(x){
    r_gen = rbinom(n = n, size = 1, prob = rho);
    y_gen = ifelse(r_gen==0,
                   rbinom(n = n,size = 1,prob = theta_0),
                   rbinom(n = n,size = 1,prob = theta_1));
    return(y_gen)
  } )

  # generate output object with samples and ground-truth for parameters
  output = list(y = y_gen,
                pi = pi,
                n = n,
                p = p, # remember that this now includes intercept
                X = X,
                X_corr = X_corr,
                small_area_id = small_area_id,
                beta = beta_temp,
                gamma = gamma,
                mu = mu,
                moranI.mu = ape::Moran.I(as.numeric(mu),nb_objects$nb[small_area_id,small_area_id]),
                moranI.gamma = ape::Moran.I(as.numeric(gamma),nb_objects$nb)
  )

  output = append(output,nb_objects)

  return(output)
}
