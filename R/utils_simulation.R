#' @param n_raters Number of Raters
#' @param n_objects Number of Objets
#' @param target_icc ICC [0,1]
#' @param p probability
#' @return data
#' @export
simulate_binary <- function(
  n_raters = 100,
  n_objects= 10,
  target_icc = 0.5,
  p=0.5,...){
  # set seed first
  #set.seed(SEED)
  
  # 1) Set binary data hyper-parameters 
  
  fixed_obj_var <- 1  # assume using rater-residual ratio for now
  intercept <- qnorm(p) # probit transformation 

  # Scenario A: Noise is mostly random error (Rater Variance is low)
  # Ratio 0.2 means: Rater Var is only 20% size of Residual Var
  rater_resid_ratio <- 0.2

  # 2) First, obtain (object and rater) random effects 
  # must be fully crossed.
 
  dat <- generate_data_ORE(n_raters,n_objects,
    target_icc, fixed_obj_var, rater_resid_ratio)
  
  # 3)) Generate datasets given the binomial model 
  # calculate Linear Predictor and Probabilities
  df <- dat %>% tibble () %>%
    mutate(
    # Map random effects to rows
    # u_i = object_effects[Object_ID],
    # v_j = rater_effects[Rater_ID],
    
    # Linear Predictor (probit scale)
    eta = intercept + u_i + v_j + Error, # don't add error? 
    
    # Probability scale (inverse probit)
    # prob = plogis(eta), 
    
    # Generate Binary Rating
    #Score = rbinom(n(), 1, prob) #bernouli (don't do since not stochastic?)


    Score = ifelse(eta > 0 , 1, 0), #threshold of the probit scale 
    #Seed = SEED
  )

  return(df)

}


run_one_binary <- function(n_raters, n_objects,target_icc, p,iter){
  #ones run row as stated in the params matrix
  #set set 
  #set.seed(seed)

  results <- simhelpers::repeat_and_stack(iter, expr ={
  #1) Generate Data 
    dat <- simulate_binary(n_raters, n_objects, target_icc, p)

  #2) Analyze data with error
    
  calc_icc_quietly <- purrr::quietly(purrr:::possibly(calc_icc, otherwise = NA))
  
  estimate_icc <- calc_icc_quietly(dat, subject ="ObjectID", rater = "RaterID", scores = "Score",
     k = NULL, engine = "LME")
  
  
  
    return(estimate_icc)

  }, stack = TRUE)
  
  return(results)
  
}

# Utilizing {simhelpers}
# bundle the data-generating function and icc function together

binary_sim <- simhelpers::bundle_sim(
  f_generate = simulate_binary, 
  f_analyze = calc_vardel_icc
  #seed = "SEED" #already set in params
) 

#' @param P Parameter grid
#' @param Iter Int; # of repetitions per condition
#' @return ICCs
#' @export
run_all_binary <- function(P, iter){
  res <- furrr::future_pmap(P, binary_sim, reps=iter,
      .progress = TRUE,
    .options = furrr::furrr_options(seed = NULL,
   packages = "vardel"))
  
  params$res <- res
  return(params)
}




