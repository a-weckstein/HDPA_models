###################################################################################################################################  
###################################################################################################################################  
## Illustrative code from Weckstein et al (2025) ------
# Includes sample R code for primary data-adaptive (DA) model types used for high-dimensional proxy adjustment
###################################################################################################################################  
###################################################################################################################################  

library(dplyr)
library(MatchIt)
library(glmnet)
library(survival)
'%!in%' <- function(x,y)!('%in%'(x,y))

###################################################################################################################################  
## Unpenalized logistic regression for IS and HDPS strategies -----
###################################################################################################################################
# Example code shown for illustrative purposes only. Study analyses used parallel computing processes with automated components to expedite model fitting across all included emulation studies and model types.      

## Inputs ----  
  ## @ predictors: includes vector of covariate names for inputs to unpenalized logistic regression model, for each covariate model.
        ## For IS covariate adjustment, "predictors" includes only the investigator-defined algorithms and covariate forms used for the original RCT-DUPLICATE initiative. Each trial emulation contained a distinct set of IS covariates.
        ## For the hdPS strategy, "predictors" includes the top 500 bias-ranked empirical features AND either a basic set of demographic variables (for fDA) or all IS covariates (fDA).
  
  ## @df: Deidentified patient-level dataframe for each trial_emulation with columns for patient identifier, exposure indicator, outcome status indicator, and all candidate covariate features (empirical and investigator-defined). 
    
  
   # Fit unpenalized IS or HDPS PS model for treatment (exposure)
          PS_fit <- glm(as.formula(paste("exposure_numeric ~ ", paste(c(predictors), collapse = " + "))),
                              data = df, family = binomial(link = 'logit'))
          summary(PS_fit)
          # Extract predictions (propensity scores):
            ps <- predict(PS_fit, type = "response")
            
###################################################################################################################################
#### Function for creating folds for cross-validation ------
###################################################################################################################################
# Creating folds used for cross-validation within Lasso, OAL, and elastic net models. 
# stratifyCVFoldsByYandID function created by Susan Gruber, as described in supplement of Wyss et al 2024.
# This function allows a specified variable to be distributed equally within the folds. This is particularly useful when analyzing data with rare outcome or exposure events.
  # For models fit to exposure, stratifyCVFoldsByYandID function used to evenly distribute treatment groups across folds.
  # For models fit to outcome, stratifyCVFoldsByYandID function used to evenly distribute outcome events across folds.


CustomFxn_stratifyCVFoldsByYandID <- function (V, Y, id = NULL) {
    if (is.null(id)) id <- 1:length(Y)
    case_status_by_id <- by(Y, id, sum) # this gives n.unique results, sorted by id # 
    case_ids <- names(case_status_by_id)[ case_status_by_id > 0]
    noncase_ids <- names(case_status_by_id)[ case_status_by_id == 0]
    if (V > min(length(case_ids), length(noncase_ids))) {
    stop("number of observations in minority class is less than the number of folds") }
    valSet.case_ids <- split(sample(case_ids), rep(1:V, length = length(case_ids))) 
    valSet.noncase_ids <- split(sample(noncase_ids), rep(1:V, length = length(noncase_ids))) 
    validRows <- vector("list", length = V)
    names(validRows) <- paste(seq(V))
    fold_id <- rep(NA, length(Y))
    for (v in seq(V)){
      validRows[[v]] <- which(as.character(id) %in% c(valSet.case_ids[[v]], valSet.noncase_ids[[v]]))
      fold_id[validRows[[v]]] <- v }
    return(list(validRows = validRows, fold_id = fold_id)) }


  ## example applications:
     set.seed(1234001)
      cvfolds_E <- CustomFxn_stratifyCVFoldsByYandID(V=10, Y = df$exposure_numeric)
      folds_E <- cvfolds_E$validRows
      foldid_E <- cvfolds_E$fold_id
      ### create folds, by OUTCOME:
      set.seed(1234001)
      cvfolds_O <- CustomFxn_stratifyCVFoldsByYandID(V=10, Y = df$outcome)
      folds_O <- cvfolds_O$validRows
      foldid_O <- cvfolds_O$fold_id
    
      
###################################################################################################################################
#### LASSO ------
###################################################################################################################################
# Basic LASSO (L1) regularization used to fit PS model for treatment assignment
# Function shown for illustrative purposes only. Study analyses used parallel computing processes with automated components to expedite model fitting across all included emulation studies and model types.      
      
## Function inputs ------
          ## @candidate_predictors: List, defined separately for each trial emulation and DA model strategy. List includes two components:
                ## candidate_predictors[[1]]: vector of covariate feature names to be considered in LASSO model for treatment assignment. For instance, for the primary model 2 hybrid data-adaptive (DA) approach, candidate_predictors includes vector with names of ~6000 candidate empirical features and all investigator-defined covariates in IS PS models. Candidate predictors vary across different trial emulations and sensitivitity analyses, based on prescreening criteria and emulation-specific removal of predictors with perfectly collinearity or zero cells. 
                ## candidate_predictors[[2]]: vector of adaptive penalty factors to indicate degree of shrinkage for candidate predictors during L1 regularization. Candidate empirical features generated from raw data were assigned a penalty factor of 1. Basic demographic variables (for fDA) and investigator-defined covariates (for hDA) were assigned penalty factors between 0 and 0.1, as to slow their shrinkage relative to empirical features. We used a more complicated function (not shown) to iterate across different penalty factors between 0 and 0.1 if there were convergence issues with initial model fitting.
    
    ## @df: Deidentified patient-level dataframe for each trial_emulation with columns for patient identifier, exposure indicator, outcome status indicator, and all candidate covariate features (empirical and investigator-defined). 
    
## FXN_LASSO ------
FXN_LASSO <- function(x) { 
  # restrict covariate data to candidate predictors for given DA model strategy:
  patient_data <- df %>% dplyr::select(all_of(c(candidate_predictors[[1]]))) 
  e <- df[["exposure_numeric"]] # extract binary exposure indicator
  Wmat = as.matrix(patient_data)
  sx = Matrix(Wmat, sparse=TRUE)
  
  vec_penaltyfactors <- candidate_predictors[[2]]  # extract vector of penalty factors corresponding to candidate covariate inputs

  # fitting adaptive lasso for treatment ----
  glmnet.LASSO <- NULL
  glmnet.LASSO <- cv.glmnet(x = sx, y = e,
                          family="binomial", type.measure="deviance", 
                           alpha = 1, # alpha=1 for L1; alpha = 0.5 for elastic net sensitivity analysis
                           nlambda = 100, penalty.factor = vec_penaltyfactors,
                           parallel = FALSE,  trace.it = 1,
                           standardize = TRUE, maxit=1000,
                           foldid = foldid_E, # exposure-based folds created separately with stratifyCVFoldsByYandID function
                           keep=TRUE, relax = FALSE)
  
  # extract predicted values (propensity scores) for all lambda values in glmnet.LASSO, for both out-of-fold and whole-dataset predictions
  gns <- predict(glmnet.LASSO, newx = Wmat, s=glmnet.LASSO$lambda, type = "response")
  gns_cv <- plogis(glmnet.LASSO$fit.preval[, 1:ncol(gns)])
  lpos <- which(glmnet.LASSO$lambda == glmnet.LASSO$lambda.min) 
  ps <- gns[,lpos]
  ps_cv <- gns_cv[,lpos] 
  model_coef_out = glmnet.LASSO$glmnet.fit$beta[,lpos]
  
  # Extract variable names for coefficients not shrunk to zero
  chosenVars_out = names(model_coef_out[model_coef_out !=0 & names(model_coef_out) != 'e'])
  
  # Assign the results
  ls_output_LASSO <- list()
  ls_output_LASSO[["input_vars"]] <- formulas_forOAL
  ls_output_LASSO[["chosen_vars"]] <- chosenVars_out
  ls_output_LASSO[["input_penalty"]] <- vec_penaltyfactors
  ls_output_LASSO[["ps"]] <- ps
  ls_output_LASSO[["ps_cv"]] <- ps_cv
  ls_output_LASSO[["PID"]] <- df[["patient_id"]]
  
  return(ls_output_LASSO)
}
    
    


###################################################################################################################################
#### Outcome adaptive LASSO (OAL) ------
###################################################################################################################################

## OAL algorithm originally proposed by Susan Shortreed and Ashkan Ertefaie (2017). 
    # We used a modified OAL algorithm to allow for model convergence in high-dimensional claims setting, based on code from Wyss et al 2024.
    # We further refined the modified OAL algorithm from Wyss et al 2024 to incorporate adaptive penalty factors for prioritizing investigator-defined covariates relative to empirical features.

# Function shown for illustrative purposes only. Study analyses used parallel computing processes with automated components to expedite model fitting across all included emulation studies and model types.
      
## Function inputs ------
          ## @candidate_predictors: List, defined separately for each trial emulation and DA model strategy. Each list includes two components:
                ## candidate_predictors[[1]]: vector of covariate feature names to be considered in OAL method. For instance, for the primary model 3 hybrid data-adaptive (DA) approach, candidate_predictors includes vector with names of ~6000 candidate empirical features and all investigator-defined covariates in IS PS models. Candidate predictors vary across different trial emulations and sensitivitity analyses, based on prescreening criteria and emulation-specific removal of predictors with perfectly collinearity or zero cells. 
                ## candidate_predictors[[2]]: vector of adaptive penalty factors to indicate degree of shrinkage for candidate predictors during L1 regularization. Candidate empirical features generated from raw data were assigned a penalty factor of 1. Basic demographic variables (for fDA) and investigator-defined covariates (for hDA) were assigned penalty factors between 0 and 0.1, as to slow their shrinkage relative to empirical features. We used a more complicated function (not shown) to iterate across different penalty factors between 0 and 0.1 until convergence was achieved for both OAL steps.
    
                ## @df: Deidentified patient-level dataframe for each trial_emulation with columns for patient identifier, exposure indicator, outcome status indicator, and all candidate covariate features (empirical and investigator-defined). 
                
## FXN_OAL ----------
FXN_OAL <- function(x) { 
  # restrict covariate data to candidate predictors for given DA model strategy:
  patient_data <- df %>% dplyr::select(all_of(c(candidate_predictors[[1]]))) 
  e <- df[["exposure_numeric"]] # extract binary exposure indicator
  Wmat_e = as.matrix(cbind(e=e, patient_data))
  sx_e = Matrix(Wmat_e, sparse=TRUE)
  
  Y = df[["outcome"]] # extract binary outcome indicator
  vec_penaltyfactors <- candidate_predictors[[2]]  # extract vector of penalty factors corresponding to candidate covariate inputs
  
  # OAL step1, fit outcome model with L1 regularization ----
  glmnet.OAL1 <- NULL
  glmnet.OAL1 <- cv.glmnet(x = sx_e, y = Y,
                           family="binomial", type.measure="deviance",
                           alpha = 1,  # alpha=1 for L1; alpha = 0.5 for elastic net sensitivity analysis
                           nlambda = 100,
           # "0" penalty factor for exposure indicator to prevent shrinkage; all other predictors assigned penalty factors as determined by "vec_penaltyfactor" 
                           penalty.factor = c(0, vec_penaltyfactors), 
           # standardizing candidate predictors within glmnet function to ensure shrinkage is applied evenly:
                           standardize = TRUE,  # outcome-based folds created separately with stratifyCVFoldsByYandID function
                           maxit= 1000, keep=FALSE, relax = FALSE, 
                           trace.it = 1,  parallel = FALSE, 
                           foldid = foldid_O) 
  
  # finding minimum lambda value (lambda that optimizes predictive performance for outcome)
  lambdapos_out1 <- which(glmnet.OAL1$lambda == glmnet.OAL1$lambda.min)
  lassocoef_out1 = glmnet.OAL1$glmnet.fit$beta[,lambdapos_out1]
  
  # extracting variable names for coefficients not shrunk to zero
  chosenVars_out1 = names(lassocoef_out1[lassocoef_out1 !=0 & names(lassocoef_out1) != 'e']) 

  # OAL step2, fitting adaptive lasso for treatment ----
  vec_penalty_factor_step2 <- rep(1, ncol(patient_data))
  # Assign penalty factor of zero** for variables selected by outcome lasso model (step 1)
  vec_penalty_factor_step2[names(patient_data) %in% chosenVars_out1] <- 0
        ## *Note: Iterative procedure (not shown) was applied to certain DA strategies for which initial models failed to converge, including elastic net formulations. Different penalty factors between 0 and 0.1 were tried in a stepwise fashion until achieving convergence.
      
  
  # fitting adaptive lasso for treatment ----
  Wmat = as.matrix(patient_data)
  sx = Matrix(Wmat, sparse=TRUE)
  glmnet.OAL2 <- NULL
  glmnet.OAL2 <- cv.glmnet(x = sx, y = e,
                           family="binomial", type.measure="deviance", 
                           alpha = 1, # alpha=1 for L1; alpha = 0.5 for elastic net sensitivity analysis
                           nlambda = 100, penalty.factor = vec_penalty_factor_step2,
                           parallel = FALSE,  trace.it = 1,
                           standardize = TRUE, maxit=1000,
                           foldid = foldid_E, # exposure-based folds created separately with stratifyCVFoldsByYandID function
                           keep=TRUE, relax = FALSE)
  
  # extract predicted values (propensity scores) for all lambda values in glmnet.OAL2, for both out-of-fold and whole-dataset predictions
  gns2 <- predict(glmnet.OAL2, newx = Wmat, s=glmnet.OAL2$lambda, type = "response")
  gns2_cv <- plogis(glmnet.OAL2$fit.preval[, 1:ncol(gns2)])
  lpos2 <- which(glmnet.OAL2$lambda == glmnet.OAL2$lambda.min) 
  model2_ps <- gns2[,lpos2]
  model2_ps_cv <- gns2_cv[,lpos2] 
  model2_coef_out = glmnet.OAL2$glmnet.fit$beta[,lpos2]
  
  # Extract variable names for coefficients not shrunk to zero in second L1 step
  chosenVars_out_model2 = names(model2_coef_out[model2_coef_out !=0 & names(model2_coef_out) != 'e'])
  
  
  # Assign the results
  ls_output_OAL <- list()
  ls_output_OAL[["input_vars"]] <- candidate_predictors[[1]]
  ls_output_OAL[["chosen_vars1"]] <- chosenVars_out1
  ls_output_OAL[["input_penalty1"]] <- candidate_predictors[[2]]
  ls_output_OAL[["penalty2"]] <- vec_penalty_factor_step2
  ls_output_OAL[["chosen_vars2"]] <- chosenVars_out2
  ls_output_OAL[["ps2"]] <- model2_ps
  ls_output_OAL[["ps2_cv"]] <- model2_ps_cv
  ls_output_OAL[["PID"]] <- df[["patient_id"]]
  
  return(ls_output_OAL)
}
    
###################################################################################################################################
#### Illustrative code for 1:1 propensity score matching  ------
###################################################################################################################################
# For primary analyses, all estimated PS fits were used to conduct 1:1 nearest neighbor propensity score matching, as done in the original RCT-DUPLICATE initiative.
# Code shown for illustrative purposes only. Study analyses used parallel computing processes with automated components to expedite matching across all included emulation studies and model types.

  ## Restrict input dataframe with relevant PS predictions to variables needed for PS-matching and Cox model fits. 
      # "ps" variable includes individual propensity score value for PS fits from IS, HDPS, LASSO, or OAL model strategies.
  temp_data <- data_PS_fits %>% select("exposure_numeric","patient_id",
                                       "outcome", "fup_time", 
                                       "ps")
      
      # conduct 1:1 PS-matching using "ps" variable from input dataframe
  m.out <- matchit(as.factor(exposure_numeric) ~ PS, data =  temp_data, method = "nearest", distance =  ps, 
                   caliper = .01, std.caliper=FALSE, replace = FALSE, m.order="closest")

  data_PS_match  <- match.data(m.out) # extract matched sample


###################################################################################################################################
#### Illustrative Cox model to estimate treatment effects ------
###################################################################################################################################

## Fit basic cox model within PS-matched sample to estimate hazard ratios and 95% CIs, as done in original RCT-DUPLICATE study.
## Note that follow-up time and censoring parameters are unique to each emulation study, and already accounted for within the "fup_time" and "outcome" variables. See preregistered RCT-DUPLICATE protocols for emulation-specific details.
# Example code shown for illustrative purposes only. Study analyses used parallel computing processes with automated components to expedite model fitting across all included emulation studies and model types. 

cox_model <- coxph(formula = Surv(fup_time, outcome) ~ exposure_numeric, data =  data_PS_match, robust=TRUE) 
    summary(cox_model)
    
      
      
###################################################################################################################################
###################################################################################################################################
    
    
