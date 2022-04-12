#' Prepare data
#'
#' Prepare data in format required for estimation procedure described in "Explaining Recruitment to Violent Extremism: A Bayesian Case-Control Approach.
#'
#' @param shape sf object: shapefile data.
#' @param survey  data.table data.frame, case-control data including common geographic ID.
#' @param shape_large.area_id_name string, large area name identifiers in the shapefile.
#' @param shape_large.area_id_num  integer, large area identifiers in the shapefile.
#' @param shape_small.area_id_name string, small area name identifiers in the shapefile.
#' @param shape_small.area_id_num integer, small area identifiers in the shapefile.
#' @param survey_small.area_id_num string, small area name identifiers in the survey.
#' @param survey_small.area_id_name integer, small area identifiers in the survey.
#' @param drop.incomplete.records logical, should the function return complete data? Defaults to `TRUE`.
#' @param colnames_X character vector, covariates definining the design matrix X. Must be numeric.
#' @param interactions_list list, each element is a string of the form "a*b" where a and be are the names of two variables in colnames_X.
#' @param scale_X string, takes values "1sd" or "2sd."
#' @param colname_y string, variable name for the outcome variable. Must be numeric.
#' @param contamination logical, should this offset account for contamination? Defaults to `TRUE`.
#' @param pi numeric, scalar defining the prevalence of the outcome in the population of interest.
#' @param large_area_shape logical, should the function return a large-area shapefile? Defaults to `TRUE`.
#'
#' @importFrom magrittr %>%
#' @importFrom spdep poly2nb card nb2mat droplinks n.comp.nb
#' @importFrom data.table as.data.table data.table
#' @importFrom stats complete.cases
#' @importFrom dplyr group_by summarise ungroup
#' @importFrom geodist geodist
#' @importFrom sf st_coordinates st_centroid
#' @importFrom INLA inla.qinv
#' @importFrom Matrix Diagonal sparseMatrix
#' @importFrom raster rowSums
#'
#' @importClassesFrom Matrix dsCMatrix
#' @importClassesFrom data.table data.table
#'
#' @return
#' @export
#'
#' @examples
data.prep = function(shape,
                     survey,
                     shape_large.area_id_name= NA,shape_large.area_id_num = NA,
                     shape_small.area_id_name = NA,shape_small.area_id_num = NA,
                     survey_small.area_id_num = NA,survey_small.area_id_name = NA,
                     drop.incomplete.records = NA,
                     colnames_X = NA,
                     interactions_list = NA,
                     scale_X = NA,
                     colname_y = NA,
                     contamination = T,
                     pi = NA,
                     large_area_shape = F
){

  if(all(is.na(c(shape_large.area_id_name,shape_large.area_id_num)))){warning("no large area ID (nominal or numeric) to search for in shapefile - preparing data exclusively at the small-area level")}

  if(!is.na(shape_large.area_id_name)){# factor the large area english names in shape file
    shape[[shape_large.area_id_name]] <- as.factor(shape[[shape_large.area_id_name]])
  }
  if(!is.na(shape_large.area_id_num)){# factor the large area english names in shape file
    shape[[shape_large.area_id_num]] <- as.factor(shape[[shape_large.area_id_num]])
  }



  if(all(is.na(c(shape_small.area_id_name,shape_small.area_id_num)))){stop("no small area ID (nominal or numeric) to search for in shapefile")}
  if(all(is.na(c(survey_small.area_id_name,survey_small.area_id_num)))){stop("no small area ID (nominal or numeric) to search for in survey data")}

  if(!is.na(survey_small.area_id_name) ){
    if( !survey_small.area_id_name %in% colnames(survey)){stop("the small area nominal ID provided does not exist in survey data")}
    if( !shape_small.area_id_name %in% colnames(shape)){stop("the small area nominal ID provided does not exist in shapefile")}
    # cleaning name id - if empty, give NA value, then factor
    survey[[survey_small.area_id_name]] = as.factor(ifelse(survey[[survey_small.area_id_name]]=="",NA,survey[[survey_small.area_id_name]]))
  }
  if(!is.na(survey_small.area_id_num)){
    if( !survey_small.area_id_num %in% colnames(survey)){stop("the small area numerical ID provided does not exist in survey data")}
    if( !shape_small.area_id_num %in% colnames(shape)){stop("the small area numerical ID provided does not exist in shapefile")}
    # cleaning num id- if empty, give NA value, then factor
    survey[[survey_small.area_id_num]] = as.factor(ifelse(survey[[survey_small.area_id_num]]=="",NA,survey[[survey_small.area_id_num]]))
  }


  # merge to add legitimate large-area to survey and/or ensure small-area IDs match shapefile
  index <- c(shape_small.area_id_name = shape_small.area_id_name,
             shape_small.area_id_num = shape_small.area_id_num,
             shape_large.area_id_name = shape_large.area_id_name,
             shape_large.area_id_num = shape_large.area_id_num)

  index = index[which(!is.na(index))]

  by.x = c(num = survey_small.area_id_num,name = survey_small.area_id_name)
  by.y = c(num = shape_small.area_id_num,name = shape_small.area_id_name)

  by.y = by.y[which(!is.na(by.x))]
  by.x = by.x[which(!is.na(by.x))]

  if(length(by.y)==0 | length(by.x)==0){stop("we cannot match numbers with names - if you provide a numerical/nominal id for the survey, you must also provide an equilvalent for the shapefile, and vice-versa") }


  survey = merge(survey,
                 as.data.table(shape[,index])[,1:length(index)],
                 by.x = by.x,
                 by.y =  by.y,
                 all.x = TRUE)

  # keep only relevant columns
  survey = survey[,c(by.x,index[-which(index %in% by.y)],colnames_X,colname_y),with=F]
  # and columns in X must all be numeric
  if(!all(survey[,lapply(.SD,class),.SDcols = c(colnames_X,colname_y)] == "numeric" | survey[,lapply(.SD,class),.SDcols = c(colnames_X,colname_y)] == "integer" )){stop("one of your variables in colnames_X or colname_y is not numeric - please change this in the survey data and re-run the data.prep function")}

  # if there are any interactions
  names(interactions_list) = ifelse(names(interactions_list)=="",gsub("\\*","\\_",unlist(interactions_list)),names(interactions_list))
  if(length(interactions_list)>0){
    # make sure all variables needed are in X
    if(!all(unlist(strsplit(unlist(interactions_list),split = "\\*")) %in% colnames_X)){stop("some variables needed for interactions are missing from colnames_X - please include all interaction variables in colnames_X")}
    # calculate the interactions
    int.temp = data.table()
    for(i in 1:length(interactions_list)){
      int.temp = cbind(int.temp,
                       apply(survey[,unlist(strsplit(interactions_list[[i]],split = "\\*")),with = F],1,function(x){x[1]*x[2]}))
    }
    names(int.temp) = names(interactions_list)
    # and add them to the survey
    survey = as.data.table(cbind(survey,int.temp))
  }

  # remove incomplete records
  if(drop.incomplete.records){
    if(sum(complete.cases(survey))>0){
      survey_complete = survey[complete.cases(survey),]
    }else{stop("cannot subset complete cases - a missing value exists for every record\nmaybe something went wrong with the matching - check that levels of shape and survey areas match (spelling could be a problem)")}
  }else{
    survey_complete = survey
  }

  # re-factor now that we have removed data
  if(!is.na(survey_small.area_id_name)){
    survey_complete[[ survey_small.area_id_name]] = as.factor(as.character(unlist(survey_complete[[ survey_small.area_id_name]])))
  }
  if(!is.na(survey_small.area_id_num)){
    survey_complete[[ survey_small.area_id_num]] = as.factor(as.character(unlist(survey_complete[[ survey_small.area_id_num]])))
  }

  # and extract area ids
  small_area_id = match(survey_complete[[by.x[1]]],shape[[by.y[1]]])
  small_area_names = shape[[by.y]]

  if(!all(is.na(c(shape_large.area_id_name,shape_large.area_id_num)))){
    by.y.large = c(num = shape_large.area_id_num,name = shape_large.area_id_name)
    by.y.large = by.y.large[which(!is.na(by.y.large))]
    #  large area is typically added on from shape from small-area
    large_area_id =  match(survey_complete[[by.y.large[1]]],levels(shape[[by.y.large[1]]]))
    large_area_names = levels(shape[[by.y.large[1]]])
    N_large_area = nlevels(shape[[by.y.large[1]]])
  }


  # get covariates in a matrix object
  X = survey_complete[,c(colnames_X,names(interactions_list)),with=F]
  mu_X = apply(X,2,mean,na.rm=T)
  sd_X = apply(X,2,sd,na.rm=T)
  # choose and apply standardizaton type
  if(scale_X=="1sd"){
    X = (X - t(array(mu_X,c(dim(X)[2],dim(X)[1]))))/t(array(sd_X,c(dim(X)[2],dim(X)[1])))
  }
  if(scale_X=="2sd"){
    X_dich = which(sapply(X,function(x){all(sort(unique(x))==c(0,1))}))
    X.tmp = (X - t(array(mu_X,c(dim(X)[2],dim(X)[1]))))/t(array(2*sd_X,c(dim(X)[2],dim(X)[1])))
    for(i in 1:length(X_dich)){
      X.tmp[[X_dich[i]]] = X[[X_dich[i]]]
    }
    X = X.tmp
  }
  X = cbind(intercept = 1,X)
  p = dim(X)[2]

  # get outcome
  if(!is.na(colname_y)){
    if(!all(sort(unique(survey_complete[[colname_y]])) == c(0,1)) | !is.numeric(survey_complete[[colname_y]])){stop("outcome must be a numeric, dichotomous variable")}
    y = survey_complete[[colname_y]]
  }else{stop("please enter a string for colname_y so we can identify the outcome in the surey data")}

  # number of observations
  n = length(y)
  # number of cases
  n1 = sum(y)
  # number of controls
  n0 = n-n1

  if(n1==n |n0==n){stop("no variance in outcome variable")}

  # Now we calculate the offset :

  # offset = log(P1/P0) =
  #        = log( n1/(pi*N) / n0/((1-pi)*N) =
  #        = log( n1/pi / n0/(1-pi) )

  # under contamination:
  # offset = log( (n1 + pi*n0)/(pi*N)) / ((n0-pi*n0)/((1-pi)*N)) ) =
  #          = log( (n1 + pi*n0)(1-pi)/(pi*(n0-pi*n0)) ) =
  #          = log( (n1 + pi*n0)(1-pi)/pi*n0(1-pi) ) =
  #          = log( (n1 + pi*n0)/pi*n0 ) =
  #          = log( n1/pi*n0 + 1 )

  if(contamination){
    theta_0_large_area = NA; theta_1_large_area = NA ; n1_large_area = NA; n0_large_area = NA; offset_large_area = NA
    theta_0 = NA; theta_1 = NA ;
    # calculate offset and mixture probabilities under contaminationa:
    # at the area-level
    # if(!is.na(pi_large_area)){
    #   n1_large_area = as.numeric(table(survey_complete[[by.y.large[1]]])[,"1"])
    #    n0_large_area = as.numeric(table(survey_complete[[by.y.large[1]]])[,"0"])
    #    offset_large_area = c(); for(i in 1:length(n1_large_area)){ offset_large_area = c(offset_large_area, log( n1_large_area[i]/(pi_large_area[i]*n0_large_area[i]) + 1 ) ) }
    #    theta_0_large_area = 0
    #    theta_1_large_area = n1_large_area/(n1_large_area+pi_large_area*n0_large_area)
    #  }else{warning("no large-area level prevalence was defined - will not calculate large-area level contaminated offest")}
    # at the population level
    if(!is.na(pi)){
      offset = log( n1/(pi*n0) + 1 )
      theta_0 = 0
      theta_1 = n1/(n1+pi*n0)
    }else{warning("no population-level prevalence was defined - will not calculate contaminated offest")}
    if(is.na(theta_0)){stop("could not calculate contaminated offset - likely missing prevalence pi")}

    # if(is.na(theta_0) & is.na(theta_0_large_area)){stop("could not calculate contaminated offset - likely missing prevalence (pi or pi_large_area)")}
  }else{
    # calculate offset under simple case-control (King and Zeng, 2001, 'prior correction'):
    # at the area-level
    #if(!is.na(pi_large_area)){
    #   offset_large_area = c(); for(i in 1:length(n1_large_area)){ offset_large_area = c(offset_large_area, log( (n1_large_area[i]/pi_large_area[i])/(n0_large_area[i]/(1-pi_large_area[i]))  ) ) }
    # }else{warning("no large-area level prevalence was defined - will not calculate large-area level offest")}
    # at the population level
    if(!is.na(pi)){
      offset = log( (n1/pi)/(n0/(1-pi))  )
    }else{warning("no population-level prevalence was defined - will not calculate offest")}
  }


  # now extract neighbourhood object
  # ensure you can take lon-lat centroids over weird spherical objects
  sf::sf_use_s2(FALSE)
  # get distance matrix
  lat_lon = st_coordinates(st_centroid(shape))
  colnames(lat_lon) = c("lat","lon")
  D = geodist(lat_lon,measure = "geodesic") #1m
  # fully connect the graph based on shortest distance
  nb = addnbs(sp.sample = shape,ID = shape[[by.y]],D=D);
  # extract graph information
  nb_graph = nb2graph(nb)

  if(length(grep("large",names(index)))>0){
    # get a shapefile with areas collapsed into large areas - useful to plot large-area effects
    shape_large = shape
    shape_large$large_area_id = shape_large[[index [grep("large",names(index))[1]]]]
    shape_large =
      shape_large %>%
      group_by(large_area_id) %>%
      summarise(geometry = sf::st_union(geometry)) %>%
      ungroup()
    names(shape_large)[1] = index [grep("large",names(index))[1]]
  }else{shape_large = NA}

  output = list(survey_complete = survey_complete,

                n = n,
                y = y,

                pi = pi,

                n1 = n1,
                n0 = n0,
                offset = offset,
                theta_1 = theta_1,
                theta_0 = theta_0,
                #offset_large_area = offset_large_area,
                #theta_1_large_area = theta_1_large_area,
                #theta_0_large_area = theta_0_large_area,

                X = X,
                p = p,
                mu_X = mu_X,
                sd_X = sd_X,
                scale_X = scale_X,

                small_area_id = small_area_id,
                small_area_names = small_area_names,
                large_area_id =  large_area_id,
                large_area_names = large_area_names,

                N_large_area = N_large_area,
                nb_object = nb,
                N_small_area = nb_graph$N,
                node1_small_area = nb_graph$node1,
                node2_small_area = nb_graph$node2,
                N_small_area_edges = nb_graph$N_edges,
                scaling_factor = scale_nb_components(poly2nb(shape,queen = TRUE))[1],
                nb.matrix  = nb2mat(neighbours = nb ,style = 'B',zero.policy = TRUE),

                shape_small = shape,
                shape_large = shape_large
  )

  return(output)
}
