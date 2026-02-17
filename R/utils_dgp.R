
generate_data_ORE <- function(n_raters,n_objects,
  target_icc, fixed_obj_var, rater_resid_ratio){

  # Obtain variances 
  sigma_sqr <- get_data_ORE(target_icc, fixed_obj_var, rater_resid_ratio)
  
  # Generate effects 
  obj_effects   <- rnorm(
    n_objects, mean = 0, sd = sqrt(sigma_sqr$var_object)
  )
  rater_effects <- rnorm(
    n_raters, mean =0, sd = sqrt(sigma_sqr$var_rater)
  )

  # Create data structure # fully crossed designs only
  dat <- expand.grid(
    ObjectID = 1:n_objects,
    RaterID= 1:n_raters
  )
 # Generate scores 
  data <- dat %>% tibble %>%
    mutate(
      u_i = obj_effects[ObjectID], #object effect
      v_j = rater_effects[RaterID], #rater effect
      Error = rnorm(n(), mean=0, sd = sqrt(sigma_sqr$var_residual)),
      OBJ_VAR = sigma_sqr$var_object,
      RATER_VAR = sigma_sqr$var_rater, 
      RES_VAR = sigma_sqr$var_residual, 
      #implicate grand mean such that
      # y_ij = mu + O_i + R_j + e_ij

     # Score = grand_mean + Effect_O + Effect_R + Error
    )

  return(data)
}


get_data_ORE <- function(target_icc, fixed_obj_var = 1, rater_resid_ratio) {

  #error check on ICC value 
    if(target_icc >= 1 | target_icc <=0) {
      stop("ICC must be between 0 and 1")
    }
  
  #calculate noise(combined rater and residual variance)
  # TODO: CHANGE RATIO to OBJECT/RATER VARIANCE

 # total_noise <- fixed_obj_var * ((1/target_icc) - 1) 
  
  #total_var <- (target_icc/(1-target_icc)) + 1

  #rho = sig2_obj/ (sig2_obj + sig2_rater + 1)
  #sig2_obj = (rho * (sig2_rater+1)) / (1-rho)

  var_object = (target_icc * (rater_resid_ratio + 1))/ (1-target_icc)
  
  #obtain rater variance as a function of ratio  
 # var_residual <- total_noise / (1 + rater_resid_ratio)

 # var_rater <- total_noise - var_residual

  #return variances 
  out <- list(
    var_object = var_object, 
    var_rater = rater_resid_ratio, 
    var_residual = 1)
  
  return(out)
}

