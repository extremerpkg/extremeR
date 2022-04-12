# Analyse prepared data
#'
#' Function to estimate a Bayesian case-control model.
#'
#' @param data # list, survey and shapefile combined list object created with the `data.prep()` function.
#' @param show_code logical, should the stan code print on the screen at the end of the run? Defaults to `TRUE`.
#' @param contamination logical, should the model include the Rota et al. (2013) style contamination layer? Defaults to `TRUE`.
#' @param offset logical, should the model include a King and Zeng (2001) style offset ? if contamination is also specified, this will be a contaminated-offset. Defaults to `TRUE`.
#' @param beta_prior string, what prior should the regression coefficients have? Choice between: "normal" and "cauchy."
#' @param small_area_prior string, what should be the small-area effects type? Choice between: "fixed," "random," "ICAR," and "BYM2." Specify `NA` for no area effects. Defaults to "fixed.
#' @param large_area_prior string, what should be the large-area effects type? Choice between: "fixed" and "random." Specify `NA` for no area effects. Defaults to `NA`.
#' @param intercept_scale integer, scale of intercept; see `?stan` for further information.
#' @param iter integer, number of draws; see `?stan` for further information.
#' @param warmup number of warmup draws; see `?stan` for further information.
#' @param thin integer, number to thin draws; see `?stan` for further information.
#' @param cores integer, number of cores to use; see `?stan` for further information.
#' @param chains integer, number of chains; see `?stan` for further information.
#' @param control list; tree depth and adapt delta; see `?stan` for further information.
#' @param verbose logical, print full stan output? Defaults to `TRUE`.
#'
#' @importFrom stats model.matrix
#' @importFrom rstan stan
#'
#' @return
#' @export
#'
#' @examples
fit = function(data = data_sim,
               contamination = F,
               show_code = T,
               offset = F,
               beta_prior = "normal",
               small_area_prior = "fixed",
               large_area_prior = NA,
               intercept_scale = NA,
               iter = 100,
               warmup = 50,
               thin = 4,
               cores = 4,
               chains = 4,
               control = list(max_treedepth = 7, adapt_delta = 0.8),
               verbose = TRUE ){

  model_skeleton = "data{
int<lower = 1> n; // total number of observations
int<lower = 1> p; // number of covariates in design matrix
int<lower = 0> y[n]; // vector of labels
matrix[n, p] X; // design matrix
}//end_data

parameters{//begin_param

}//end_param

transformed parameters{//begin_transform
vector[n] mu;
mu = X * beta;
}//end_transform

model{//begin_model

// likelihood
y  ~  bernoulli_logit(mu);
}//end_model

generated quantities{//begin_gen
vector[n] log_lik;
vector[n] y_gen;

for (i in 1:n) {
log_lik[i] = bernoulli_logit_lpmf(y[i] | mu[i]);
y_gen[i] = bernoulli_logit_rng(mu[i]);
} }//end_gen
";

  # # # parameters to monitor
  pars = c("log_lik","y_gen","mu")

  # # # define beta prior

  if(beta_prior=="normal"){

    pars = c(pars,"beta")

    model_skeleton =
      gsub("\\}//end_param",
           "vector[p] beta;\n}//end_param",
           model_skeleton)

    if(!is.na(intercept_scale)){
      model_skeleton =
        gsub("// likelihood",
             paste("beta[1] ~ normal(0,",intercept_scale,");\n//beta ~ normal(0,1);\n// likelihood",sep=""),
             model_skeleton)
    }else{
      model_skeleton =
        gsub("// likelihood",
             "beta ~ normal(0,1);\n// likelihood",
             model_skeleton)
    }

  };
  # cauchy by gamma: https://betanalpha.github.io/assets/case_studies/fitting_the_cauchy.html
  if(beta_prior=="cauchy"){

    pars = c(pars,"beta")

    model_skeleton =
      gsub("begin_param",
           "begin_param\nvector[p] aux_a;\nvector<lower = 0>[p] aux_b;\n",
           model_skeleton)

    model_skeleton =
      gsub("begin_transform",
           "begin_transform\nvector[p] beta = aux_a ./ sqrt(aux_b);\n",
           model_skeleton)

    if(!is.na(intercept_scale)){
      model_skeleton =
        gsub("// likelihood",
             paste("aux_a ~ normal(0,1);\naux_b[1] ~ gamma(0.5,",intercept_scale^2,"*0.5);\naux_b ~ gamma(0.5,0.5);\n// likelihood",sep=""),
             model_skeleton)
    }else{
      model_skeleton =
        gsub("// likelihood",
             paste("aux_a ~ normal(0,1);\naux_b ~ gamma(0.5,0.5);\n// likelihood",sep=""),
             model_skeleton)
    }
  };

  # # # define small_area effects prior
  if(!is.na(small_area_prior)){
    if(small_area_prior=="fixed"){

      pars = c(pars,"gamma")

      data$Z = model.matrix(~as.factor(data$small_area_id)-1)
      colnames(data$Z) = gsub("as\\.factor\\(data\\$small_area_id\\)","",colnames(data$Z))
      for(i in 1:data$N_small_area){
        if(! as.character(i) %in% colnames(data$Z)){
          temp = data.table(x =  rep(0,data$n))
          names(temp) = as.character(i)
          data$Z = cbind(data$Z,temp)
        }
      }
      data$Z = as.matrix(data$Z)[,order(as.numeric(colnames(data$Z)))]
      data$pZ = dim(data$Z)[2]


      model_skeleton =
        gsub("\\}//end_data",
             "int<lower = 1> pZ;\nmatrix[n, pZ] Z;\nint<lower = 1> N_small_area;}//end_data",
             model_skeleton)

      model_skeleton =
        gsub("\\}//end_param",
             "vector[pZ] gamma;\n}//end_param",
             model_skeleton)

      model_skeleton =
        gsub("mu =",
             "mu = Z * gamma +",
             model_skeleton)

      model_skeleton =
        gsub("// likelihood",
             "gamma ~ normal(0,1);\nsum(gamma) ~ normal(0, 0.01 * N_small_area);// likelihood",
             model_skeleton)
    };
    if(small_area_prior=="random"){

      pars = c(pars,"gamma","sigma_gamma")

      model_skeleton =
        gsub("\\}//end_data",
             "int<lower = 1> small_area_id[n];\nint<lower = 1> N_small_area; \n}//end_data",
             model_skeleton)

      model_skeleton =
        gsub("\\}//end_param",
             "vector[N_small_area] gamma;\n real<lower = 0> sigma_gamma;\n}//end_param",
             model_skeleton)

      model_skeleton =
        gsub("mu =",
             "mu = gamma[small_area_id]*sigma_gamma +",
             model_skeleton)

      model_skeleton =
        gsub("// likelihood",
             "gamma ~ normal(0,1);\nsigma_gamma ~ normal(0,1);\n// likelihood",
             model_skeleton)
    };
    if(small_area_prior=="BYM2"){

      pars = c(pars,"gamma","sigma_gamma","lambda","phi","psi")

      model_skeleton =
        gsub("\\}//end_data",
             "int<lower = 1> small_area_id[n];
int<lower = 1> N_small_area;
int<lower = 1> N_small_area_edges;
int<lower=1, upper=N_small_area> node1_small_area[N_small_area_edges];
int<lower=1, upper=N_small_area> node2_small_area[N_small_area_edges];
real scaling_factor;
}//end_data",
             model_skeleton)

      model_skeleton =
        gsub("\\}//end_param",
             "vector[N_small_area] phi;
vector[N_small_area] psi;
real<lower = 0,upper = 1> lambda;
real<lower = 0> sigma_gamma;\n}//end_param",
             model_skeleton)

      model_skeleton =
        gsub("mu =",
             "vector[N_small_area] gamma = (sqrt(1-lambda) * phi + sqrt(lambda / scaling_factor) * psi)*sigma_gamma;\n
mu = gamma[small_area_id] +",
             model_skeleton)

      model_skeleton =
        gsub("// likelihood",
             "target += -0.5 * dot_self(psi[node1_small_area] - psi[node2_small_area]);
phi ~ normal(0,1);
sum(psi) ~ normal(0, 0.01 * N_small_area);
lambda ~ beta(0.5,0.5);
sigma_gamma ~ normal(0,1);
// likelihood",
             model_skeleton)
    };
    if(small_area_prior=="ICAR"){

      pars = c(pars,"gamma","sigma_gamma")

      model_skeleton =
        gsub("\\}//end_data",
             "int<lower = 1> small_area_id[n];
int<lower = 1> N_small_area;
int<lower = 1> N_small_area_edges;
int<lower=1, upper=N_small_area> node1_small_area[N_small_area_edges];
int<lower=1, upper=N_small_area> node2_small_area[N_small_area_edges];
real scaling_factor;
}//end_data",
             model_skeleton)

      model_skeleton =
        gsub("\\}//end_param",
             "vector[N_small_area] psi;
real<lower = 0> sigma_gamma;\n}//end_param",
             model_skeleton)

      model_skeleton =
        gsub("mu =",
             "vector[N_small_area] gamma = psi*sigma_gamma;\n
mu = gamma[small_area_id] +",
             model_skeleton)

      model_skeleton =
        gsub("// likelihood",
             "target += -0.5 * dot_self(psi[node1_small_area] - psi[node2_small_area]);
sum(psi) ~ normal(0, 0.01 * N_small_area);
sigma_gamma ~ normal(0,1);
// likelihood",
             model_skeleton)
    };
  }
  # # # define large_area effects prior
  if(!is.na(large_area_prior)){
    if(large_area_prior=="fixed"){

      pars = c(pars,"eta")

      data$Q = model.matrix(~as.factor(data$large_area_id)-1)
      colnames(data$Q) = gsub("as\\.factor\\(data\\$large_area_id\\)","",colnames(data$Q))
      data$pQ = dim(data$Q)[2]


      model_skeleton =
        gsub("\\}//end_data",
             "int<lower = 1> pQ;\nmatrix[n, pQ] Q; \n}//end_data",
             model_skeleton)

      model_skeleton =
        gsub("\\}//end_param",
             "vector[pQ] eta;\n}//end_param",
             model_skeleton)

      model_skeleton =
        gsub("mu =",
             "mu = Q * eta +",
             model_skeleton)

      model_skeleton =
        gsub("// likelihood",
             "eta ~ normal(0,1);\n// likelihood",
             model_skeleton)
    };
    if(large_area_prior=="random"){

      pars = c(pars,"eta","sigma_eta")

      model_skeleton =
        gsub("\\}//end_data",
             "int<lower = 1> large_area_id[n];\nint<lower = 1> N_large_area; \n}//end_data",
             model_skeleton)

      model_skeleton =
        gsub("\\}//end_param",
             "vector[N_large_area] eta;\n real<lower = 0> sigma_eta;\n}//end_param",
             model_skeleton)

      model_skeleton =
        gsub("mu =",
             "mu = eta[large_area_id]*sigma_eta +",
             model_skeleton)

      model_skeleton =
        gsub("// likelihood",
             "eta ~ normal(0,1);\nsigma_eta ~ normal(0,1);\n// likelihood",
             model_skeleton)
    };
  }
  # # # define mu with offset
  if(!is.null(data$pi_large_area) & any(is.na(data$pi_large_area))){warning("some missing values in pi_large_areas - cannot be used to make offset")}
  if(offset & all(is.null(data$pi_large_area))){

    if(contamination){ data$log_offset = log(sum(data$y)/(data$pi*(data$n-sum(data$y))) +1) }else{
      data$log_offset = log(((1-data$pi)/data$pi)*(mean(data$y)/(1-mean(data$y)))) }

    model_skeleton =
      gsub("\\}//end_data",
           "real log_offset;\n}//end_data",
           model_skeleton)

    model_skeleton =
      gsub("mu =",
           "mu = log_offset +",
           model_skeleton)
  };
  if(offset & all(!is.null(data$pi_large_area))){

    if(contamination){
      if(any(names(data$offset_large_area)!=data$large_area_names)){stop("Order of offsets in data doesn't match order of levels in large-area.")}
      data$log_offset = data$offset_large_area
    }else{
      if(any(names(data$offset_large_area_no.contamination)!=data$large_area_names)){stop("Order of offsets in data doesn't match order of levels in large-area.")}
      data$log_offset = data$offset_large_area_no.contamination
    }

    model_skeleton =
      gsub("\\}//end_data",
           "vector[N_large_area] log_offset;\n}//end_data",
           model_skeleton)

    model_skeleton =
      gsub("mu =",
           "mu = log_offset[large_area_id] +",
           model_skeleton)
  };
  # # # redefine the model with contamination
  if(contamination & all(is.na(data$pi_large_area))){
    data$theta = c()
    data$theta[1] = 0
    data$theta[2] = sum(data$y)/(sum(data$y) + data$pi*(data$n-sum(data$y)))

    model_skeleton =
      gsub("\\}//end_data",
           "vector[2] theta;\n}//end_data",
           model_skeleton)

    model_skeleton =
      gsub("y  ~  bernoulli_logit\\(mu\\);",
           "for (i in 1:n) {
  target += log_mix(1-inv_logit(mu[i]),
                    bernoulli_lpmf(y[i] | theta[1]),
                    bernoulli_lpmf(y[i] | theta[2]));
                    }",
           model_skeleton)

    model_skeleton =
      gsub("vector\\[n\\] y\\_gen;",
           "vector[n] y_gen;\nint r_gen[n];",
           model_skeleton)

    model_skeleton =
      gsub("log_lik\\[i\\] = bernoulli_logit_lpmf\\(y\\[i\\] \\| mu\\[i\\]\\);",
           "log_lik[i] = log_mix(1-inv_logit(mu[i]),bernoulli_lpmf(y[i] | theta[1] ),bernoulli_lpmf(y[i] | theta[2] ) );",
           model_skeleton)

    model_skeleton =
      gsub("y_gen\\[i\\] = bernoulli_logit_rng\\(mu\\[i\\]\\);",
           "r_gen[i] = bernoulli_rng(inv_logit(mu[i]));\ny_gen[i] = bernoulli_rng(theta[r_gen[i]+1]);",
           model_skeleton)

  };
  if(contamination & all(!is.null(data$pi_large_area))){
    data$theta = as.matrix(t(data.table(theta_0 = data$theta_0_large_area,
                                        theta_1 = data$theta_1_large_area)))

    model_skeleton =
      gsub("\\}//end_data",
           "matrix[2,N_large_area] theta;\n}//end_data",
           model_skeleton)

    model_skeleton =
      gsub("y  ~  bernoulli_logit\\(mu\\);",
           "for (i in 1:n) {
  target += log_mix(1-inv_logit(mu[i]),
                    bernoulli_lpmf(y[i] | theta[1,large_area_id[i]]),
                    bernoulli_lpmf(y[i] | theta[2,large_area_id[i]]));
                    }",
           model_skeleton)

    model_skeleton =
      gsub("vector\\[n\\] y\\_gen;",
           "vector[n] y_gen;\nint r_gen[n];",
           model_skeleton)

    model_skeleton =
      gsub("log_lik\\[i\\] = bernoulli_logit_lpmf\\(y\\[i\\] \\| mu\\[i\\]\\);",
           "log_lik[i] = log_mix(1-inv_logit(mu[i]),bernoulli_lpmf(y[i] | theta[1,large_area_id[i]]),bernoulli_lpmf(y[i] | theta[2,large_area_id[i]]) );",
           model_skeleton)

    model_skeleton =
      gsub("y_gen\\[i\\] = bernoulli_logit_rng\\(mu\\[i\\]\\);",
           "r_gen[i] = bernoulli_rng(inv_logit(mu[i]));\ny_gen[i] = bernoulli_rng(theta[r_gen[i]+1,large_area_id[i]]);",
           model_skeleton)

  };


  # # # avoid stan recompiling of code if possible - do not amend code string
  if(exists("model_skeleton_store")){
    if(model_skeleton == model_skeleton_store){

    }else{model_skeleton_store <<- model_skeleton}}else{model_skeleton_store <<- model_skeleton}

  fit <- stan(model_code = model_skeleton_store,
              data = data,
              iter = iter,
              warmup = warmup,
              thin = thin,
              cores = cores,
              pars = pars,
              chains = chains,
              control = control,
              verbose = verbose);


  cat(model_skeleton_store)
  # # # should we show the code at the end ?
  if(show_code == T){ cat(model_skeleton_store) }

  return(fit)
}
# # # output: a STAN model fit object , which can be analysed with traditional STAN functions.
