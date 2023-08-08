######################
# FSF Grant Power Analysis
# Ryan Larson, PhD
#####################

library(tidyverse)
library(tigris)
library(tidycensus)
library(sf)
library(sfdep)
library(spdep)
library(spatialreg)
library(broom)
library(fabricatr)
library(plm)
library(boot)


set.seed(7188)

# ZCTA Base
census_api_key("ecda17575f4d914b502c70f2bae7a5f3d253792d")

zcta <- get_acs(geography = "zcta",
                variables = "B01001_001",
                output = "wide", 
                year = 2020, #change back to 2020
                geometry = T,
                survey = "acs5") %>%
  rename(zcta = GEOID,
         pop_2020 = B01001_001E) %>%
  select(-c(NAME, B01001_001M)) %>%
  mutate(zcta = as.numeric(zcta))


#city boundaries
mpls_stp <- places(state = "MN", year = 2020, cb = F) %>%
  filter(NAME=="Minneapolis"|NAME=="St. Paul") 

#intersecting ZCTAs
zcta_msp <- zcta %>%
  st_filter(mpls_stp, .predicate = st_intersects)

#map to check
ggplot(zcta_msp)+
  geom_sf(aes(geometry=geometry))

# Aim 1 simulation

#power function
aim_1_power <- function(effect, sample_size, alpha, nsim, lambda, srmu, rho){
  
  sig_results <- c()
  
  for (i in 1:nsim) {
    
    #simulate data
    aim_1_data <- zcta_msp %>%
      mutate(sr = rnorm(sample_size, srmu, sqrt(srmu*(1-srmu)/sample_size)), 
             ps_init = rpois(sample_size, lambda = lambda),
             ps_rate_init = ifelse(pop_2020==0,
                                   NA_integer_,
                                   ps_init/pop_2020))
    
    #create spatial weights matrix
    nb <- st_contiguity(aim_1_data, queen=TRUE)
    st_weight <- st_weights(nb, style =  "W") #spatial weights
    wt <- nb2listw(nb, style = "W") #neighbor list object
    #create spatial lag variable
    sp_lag <- sfdep::st_lag(aim_1_data$ps_init, nb, st_weight)
    
    #create DV with specified effects
    aim_1_data <- aim_1_data %>% mutate(
      sp_lag = sp_lag,
      ps = ifelse(pop_2020==0,
                  NA_integer_,
                  effect*sr + rho*sp_lag + (rpois(sample_size, lambda = lambda)/pop_2020)))
    
    #Spatial AR(1) Model
    sp_ar <- lagsarlm(ps~sr, data = aim_1_data, listw=wt)
    
    #extract p-value
    sig_results[i] <- tidy(sp_ar)$p.value[3] <= alpha
  }
  
  sig_results %>%
    mean() %>%
    return()
}



#iterate simulation over varying effect sizes

power_levels <- c()

effects <- seq(from = 0, to = 1, by = .1)


for (i in 1:length(effects)) {
  power_levels[i] <- aim_1_power(effect = effects[i], sample_size = 50, 
                                        alpha = .05, nsim = 100, 
                                 lambda = 1, srmu = .5, rho = .05)
}

power_results <- tibble(effect = effects,
                        power = power_levels)

power_results

ggplot(power_results, 
       aes(x = effect, y = power)) +
  geom_line(color = 'red', size = 1.5) + 
  # horizontal line at 80% power
  geom_hline(aes(yintercept = .8), linetype = 'dashed') + 
  theme_minimal() + 
  scale_y_continuous(labels = scales::percent) + 
  labs(x = 'Linear Effect Size', y = 'Power')



# Aim 2 simulation

#panel data base

zcta_panel <- bind_rows(replicate(15, zcta_msp, simplify = F)) %>%
  arrange(zcta) %>%
  mutate(year = rep(2008:2022, 50))

#power function
aim_2_power <- function(effect, n, t, alpha, nsim, sd){
  
  sig_results <- c()
  
  for (i in 1:nsim) {
    
    #simulate data
    aim_2_data <- fabricate(
      zctas = add_level(N= n, zcta_fe = runif(N, 1, 10)),
      years = add_level(N = t, year_shock = runif(N, 1, 10), nest = FALSE),
      observations = cross_levels(
        by = join_using(zctas, years),
      ps =  zcta_fe + year_shock + rnorm(N, 0, sd),
      health = effect*ps + zcta_fe + year_shock + rnorm(N, 0, sd)
    ))
    
    #Twoway FE Panel Model
    fe <- plm(health~ps, index = c("zcta_fe", "year_shock"),
                 model = "within", effect = "twoways",
                 data = aim_2_data)
    
    #extract p-value
    sig_results[i] <- tidy(fe)$p.value[1] <= alpha
  }
  
  sig_results %>%
    mean() %>%
    return()
}

#iterate simulation over varying effect sizes

power_levels <- c()

effects <- seq(from = 0, to = 1, by = .01)


for (i in 1:length(effects)) {
  power_levels[i] <- aim_2_power(effect = effects[i], n= 50, t = 50, 
                                 alpha = .05, nsim = 100, sd = 3)
}

power_results <- tibble(effect = effects,
                        power = power_levels)

power_results

ggplot(power_results, 
       aes(x = effect, y = power)) +
  geom_line(color = 'red', size = 1.5) + 
  # horizontal line at 80% power
  geom_hline(aes(yintercept = .8), linetype = 'dashed') + 
  theme_minimal() + 
  scale_y_continuous(labels = scales::percent) + 
  labs(x = 'Linear Effect Size', y = 'Power')

# Aim 3 simulation

cmest <- function(data = NULL, model = "rb",
                  full = TRUE, casecontrol = FALSE, yrare = NULL, yprevalence = NULL,
                  estimation = "imputation", inference = "bootstrap",
                  outcome = NULL, event = NULL,
                  exposure = NULL, mediator = NULL, EMint = NULL, basec = NULL, postc = NULL,
                  yreg = NULL, mreg = NULL, wmnomreg = NULL, wmdenomreg = NULL, ereg = NULL, 
                  postcreg = NULL,
                  astar = 0, a = 1, mval = NULL, yval = NULL, basecval = NULL,
                  nboot = 200, boot.ci.type = "per", nRep = 5, multimp = FALSE, ...) {
  # function call
  cl <- match.call()
  # output list
  out <- list(call = cl)
  
  ###################################################################################################
  #################################Argument Restrictions And Warnings################################
  ###################################################################################################
  # data
  if (is.null(data)) stop("Unspecified data")
  data <- as.data.frame(data)[, c(outcome, event, exposure, mediator, basec, postc)]
  if (sum(is.na(data)) > 0 && !multimp) stop("NAs in outcome, event, exposure, mediator, basec, or postc data; delete rows with NAs in these variables from the data or set multimp = TRUE") 
  out$data <- data
  n <- nrow(data)
  # model
  if (!model %in% c("rb", "wb", "iorw", "ne", "gformula", "msm")) stop("Select model from 'rb', 'wb', 'iorw', 'ne', 'gformula', 'msm'")
  if (model == "ne") stop("The natural effect model ('ne') is temporarily unavailable")
  out$methods$model <- model
  # full
  if (!is.logical(full)) stop("full should be TRUE or FALSE")
  out$methods$full <- full
  # casecontrol, yrare, yprevalence
  if (!is.logical(casecontrol)) stop("casecontrol should be TRUE or FALSE")
  out$methods$casecontrol <- casecontrol
  if (casecontrol) {
    if (length(unique(data[, outcome])) != 2) stop("When casecontrol is TRUE, the outcome must be binary")
    if (is.null(yprevalence) && yrare != TRUE) stop("When casecontrol is TRUE, specify yprevalence or set yrare to be TRUE")
    # imputation-based estimation is biased for case-control studies without specifying yprevalence 
    if (is.null(yprevalence) && !(estimation %in% c("para", "paramfunc"))) stop("When casecontrol is TRUE, specify yprevalence or use estimation = 'paramfunc'")
    if (!is.null(yprevalence)) {
      if (!is.numeric(yprevalence)) stop("yprevalence should be numeric")
      out$methods$yprevalence <- yprevalence
    } else out$methods$yrare <- yrare
  }
  # estimation, inference, nboot
  if (estimation == "para") estimation <- "paramfunc"
  if (estimation == "impu") estimation <- "imputation"
  if (inference == "delt") inference <- "delta"
  if (inference == "boot") inference <- "bootstrap"
  if (model == "rb" && !estimation %in% c("paramfunc", "imputation")) stop("When model = 'rb', select estimation from 'paramfunc', 'imputation'")
  if (model != "rb" && !estimation == "imputation") stop("Use estimation = 'imputation'")
  if (estimation == "paramfunc") {
    if (!is.character(yreg) | length(yreg) != 1) stop(
      "When estimation = 'paramfunc', select yreg from 'linear', 'logistic', 
                       'loglinear', 'poisson', 'quasipoisson', 'negbin', 'coxph', 'aft_exp', 'aft_weibull'")
    if (!yreg %in% c("linear", "logistic", "loglinear", "poisson",
                     "quasipoisson", "negbin", "coxph", "aft_exp", "aft_weibull")) stop(
                       "When estimation = 'paramfunc', select yreg from 'linear', 'logistic', 
                       'loglinear', 'poisson', 'quasipoisson', 'negbin', 'coxph', 'aft_exp', 'aft_weibull'")
    if (length(mediator) > 1) stop("'paramfunc' only supports a single mediator")
    if (!is.character(mreg[[1]]) | length(mreg[[1]]) != 1) stop(
      "When estimation = 'paramfunc', select mreg[[1]] from 'linear', 'logistic', 'multinomial'")
    if(!mreg[[1]] %in% c("linear", "logistic", "multinomial")) stop(
      "When estimation = 'paramfunc', select mreg[[1]] from 'linear', 'logistic', 'multinomial'")
  }
  out$methods$estimation <- estimation
  if (!inference %in% c("delta", "bootstrap")) stop("Select inference from 'delta', 'bootstrap'")
  if (estimation == "imputation" && inference == "delta") stop("Use inference = 'bootstrap' when estimation = 'imputation'")
  out$methods$inference <- inference
  if (inference == "bootstrap") {
    if (!is.numeric(nboot)) stop("nboot should be numeric")
    if (!boot.ci.type %in% c("per", "bca")) stop("Select boot.ci.type from 'per', 'bca'")
    out$methods$nboot <- nboot
    out$methods$boot.ci.type <- boot.ci.type
  }
  # outcome
  if (length(outcome) == 0) stop("Unspecified outcome")
  if (length(outcome) > 1) stop("length(outcome) > 1")
  out$variables$outcome <- outcome
  # event
  if (is.character(yreg)) {
    if (yreg %in% c("coxph", "aft_exp", "aft_weibull") && !is.null(event)) out$variables$event <- event
  } else {
    if (!is.null(event)) warning("event is ignored when yreg is not character")
  }
  # exposure
  if (length(exposure) == 0) stop("Unspecified exposure")
  if (length(exposure) > 1) stop("length(exposure) > 1")
  out$variables$exposure <- exposure
  # mediator
  if (length(mediator) == 0) stop("Unspecified mediator")
  out$variables$mediator <- mediator
  # EMint
  if (model != "iorw") {
    if (!is.logical(EMint)) stop("EMint should be TRUE or FALSE")
    out$variables$EMint <- EMint
  } 
  # basec
  if (!is.null(basec)) out$variables$basec <- basec
  # postc
  if (length(postc) != 0) {
    if (!model %in% c("msm", "gformula")) stop("When postc is not empty, select model from 'msm' and 'gformula'")
    out$variables$postc <- postc
  }
  #regs
  if (model == "rb") out$reg.input <- list(yreg = yreg, mreg = mreg)
  if (model == "wb") out$reg.input <- list(yreg = yreg)
  if (model == "gformula") out$reg.input <- list(yreg = yreg, mreg = mreg)
  if (model == "msm") out$reg.input <- list(yreg = yreg, mreg = mreg, wmnomreg = wmnomreg, wmdenomreg = wmdenomreg)
  if (model == "iorw") out$reg.input <- list(yreg = yreg, ereg = ereg)
  if (model == "ne") out$reg.input <- list(yreg = yreg) 
  if (length(basec) != 0 && model %in% c("wb", "msm")) out$reg.input$ereg <- ereg
  if (length(postc) != 0 && model == "gformula") out$reg.input$postcreg <- postcreg
  # a, astar
  if (is.null(a) | is.null(astar)) stop("Unspecified a or astar")
  # mval
  if (model != "iorw") {
    if (!is.list(mval)) stop("mval should be a list")
    if (length(mval) != length(mediator)) stop("length(mval) != length(mediator)")
    for (p in 1:length(mval)) if (is.null(mval[[p]])) stop(paste0("Unspecified mval[[", p, "]]"))
  } 
  # basecval
  if (!(model == "rb" && estimation == "paramfunc" && length(basec) != 0) && !is.null(basecval)) warning("basecval is ignored")
  # nRep
  if (model == "ne") {
    if (!is.numeric(nRep)) stop("nRep should be numeric")
    out$methods$nRep <- nRep
  }
  # multimp
  out$multimp <- list(multimp = multimp)
  if (!is.logical(multimp)) stop("multimp should be TRUE or FALSE")
  # args_mice
  if (multimp) {
    args_mice <- list(...)
    args_mice$print <- FALSE
    if (!is.null(args_mice$data)) warning("args_mice$data is overwritten by data")
    args_mice$data <- data
    out$multimp$args_mice <- args_mice
  }
  if (!multimp && model %in% c("wb", "iorw", "ne", "msm")) {
    if (sum(is.na(data)) > 0) stop("Selected model doesn't support missing values; use multimp = TRUE")
  }
  
  # run regressions
  environment(regrun) <- environment()
  regs <- regrun()
  yreg <- regs$yreg
  ereg <- regs$ereg
  mreg <- regs$mreg
  wmnomreg <- regs$wmnomreg
  wmdenomreg <- regs$wmdenomreg
  postcreg <- regs$postcreg
  
  ###################################################################################################
  ############################################Estimation and Inference###############################
  ###################################################################################################
  # add a progress bar for bootstrap inference
  if (inference == "bootstrap") {
    env <- environment()
    counter <- 0
    progbar <- txtProgressBar(min = 0, max = nboot, style = 3)
  }
  # estimation and inference of causal effects
  environment(estinf) <- environment()
  out <- c(out, estinf())
  class(out) <- "cmest"
  return(out)
}


#' @describeIn cmest Print the results of \code{cmest} nicely
#' @export
print.cmest <- function(x, ...) {
  cat("Causal Mediation Analysis\n\n")
  # print regression models used
  if (!x$multimp$multimp) {
    regnames <- names(x$reg.output)
    for (name in regnames) {
      if (name == "yreg") {
        cat("# Outcome regression:\n")
        if (inherits(x$reg.output$yreg, "svyglm")) {
          x$reg.output$yreg$call <- update(x$reg.output$yreg,design = getCall(x$reg.output$yreg)$design,
                                           family = getCall(x$reg.output$yreg)$family, evaluate = FALSE)
          x$reg.output$yreg$survey.design$call <- as.call(update(summary(x$reg.output$yreg)$survey.design,
                                                                 data = getCall(summary(x$reg.output$yreg)$survey.design)$data,
                                                                 weights = getCall(summary(x$reg.output$yreg)$survey.design)$weights, 
                                                                 evaluate = FALSE))
          print(x$reg.output$yreg)
        } else {
          x$reg.output$yreg$call <- update(x$reg.output$yreg,data=getCall(x$reg.output$yreg)$data,
                                           weights=getCall(x$reg.output$yreg)$weights, evaluate = FALSE)
          print(x$reg.output$yreg)
        }
      }
      if (name == "yregTot") {
        cat("# Outcome regression for the total effect: \n")
        x$reg.output$yregTot$call <- update(x$reg.output$yregTot,data=getCall(x$reg.output$yregTot)$data,
                                            weights=getCall(x$reg.output$yregTot)$weights, evaluate = FALSE)
        print(x$reg.output$yregTot)
      }
      if (name == "yregDir") {
        cat("# Outcome regression for the direct effect: \n")
        x$reg.output$yregDir$call <- update(x$reg.output$yregDir,data=getCall(x$reg.output$yregDir)$data,
                                            weights=getCall(x$reg.output$yregDir)$weights, evaluate = FALSE)
        print(x$reg.output$yregDir)
      }
      if (name == "ereg") {
        cat("# Exposure regression for weighting: \n")
        x$reg.output$ereg$call <- update(x$reg.output$ereg,data=getCall(x$reg.output$ereg)$data,
                                         weights=getCall(x$reg.output$ereg)$weights, evaluate = FALSE)
        print(x$reg.output$ereg)
      }
      if (name == "mreg") {
        cat("# Mediator regressions: \n")
        for (i in 1:length(x$reg.output$mreg)) {
          if (inherits(x$reg.output$mreg[[i]], "svyglm")) {
            x$reg.output$mreg[[i]]$call <- eval(bquote(update(x$reg.output$mreg[[.(i)]],
                                                              design = getCall(x$reg.output$mreg[[.(i)]])$design,
                                                              family = getCall(x$reg.output$mreg[[.(i)]])$family, evaluate = FALSE)))
            x$reg.output$mreg[[i]]$survey.design$call <- eval(bquote(as.call(update(summary(x$reg.output$mreg[[.(i)]])$survey.design,
                                                                                    data = getCall(summary(x$reg.output$mreg[[.(i)]])$survey.design)$data,
                                                                                    weights = getCall(summary(x$reg.output$mreg[[.(i)]])$survey.design)$weights, 
                                                                                    evaluate = FALSE))))
            print(x$reg.output$mreg[[i]])
          } else {
            x$reg.output$mreg[[i]]$call <- eval(bquote(update(x$reg.output$mreg[[.(i)]], 
                                                              data=getCall(x$reg.output$mreg[[.(i)]])$data, 
                                                              weights=getCall(x$reg.output$mreg[[.(i)]])$weights,
                                                              evaluate = FALSE)))
            print(x$reg.output$mreg[[i]])
          }
          if (i < length(x$reg.output$mreg)) cat("\n")
        }
      }
      if (name == "wmdenomreg") {
        cat("# Mediator regressions for weighting (denominator): \n")
        for (i in 1:length(x$reg.output$wmdenomreg)) {
          x$reg.output$wmdenomreg[[i]]$call <- eval(bquote(update(x$reg.output$wmdenomreg[[i]], 
                                                                  data=getCall(x$reg.output$wmdenomreg[[.(i)]])$data, 
                                                                  weights=getCall(x$reg.output$wmdenomreg[[.(i)]])$weights,
                                                                  evaluate = FALSE)))
          print(x$reg.output$wmdenomreg[[i]])
          if (i < length(x$reg.output$wmdenomreg)) cat("\n")
        }
      }
      if (name == "wmnomreg") {
        cat("# Mediator regressions for weighting (nominator): \n")
        for (i in 1:length(x$reg.output$wmnomreg)) {
          x$reg.output$wmnomreg[[i]]$call <- eval(bquote(update(x$reg.output$wmnomreg[[i]], 
                                                                data=getCall(x$reg.output$wmnomreg[[.(i)]])$data, 
                                                                weights=getCall(x$reg.output$wmnomreg[[.(i)]])$weights,
                                                                evaluate = FALSE)))
          print(x$reg.output$wmnomreg[[i]])
          if (i < length(x$reg.output$wmnomreg)) cat("\n")
        }
      }
      if (name == "postcreg") {
        cat("# Regressions for mediator-outcome confounders affected by the exposure: \n")
        for (i in 1:length(x$reg.output$postcreg)) {
          x$reg.output$postcreg[[i]]$call <- eval(bquote(update(x$reg.output$postcreg[[i]], 
                                                                data=getCall(x$reg.output$postcreg[[.(i)]])$data, 
                                                                weights=getCall(x$reg.output$postcreg[[.(i)]])$weights,
                                                                evaluate = FALSE)))
          print(x$reg.output$postcreg[[i]])
          if (i < length(x$reg.output$postcreg)) cat("\n")
        }
      }
      cat("\n")
    }
  } else {
    for (m in 1:length(x$reg.output)){ 
      cat(paste("# Regressions with imputed dataset", m, "\n\n"))
      regnames <- names(x$reg.output[[m]])
      for (name in regnames) {
        if (name == "yreg") {
          cat("## Outcome regression: \n")
          if (inherits(x$reg.output[[m]]$yreg, "svyglm")) {
            x$reg.output[[m]]$yreg$call <- eval(bquote(update(x$reg.output[[.(m)]]$yreg,design = getCall(x$reg.output[[.(m)]]$yreg)$design,
                                                              family = getCall(x$reg.output[[.(m)]]$yreg)$family, evaluate = FALSE)))
            x$reg.output[[m]]$yreg$survey.design$call <- eval(bquote(as.call(update(summary(x$reg.output[[.(m)]]$yreg)$survey.design,
                                                                                    data = getCall(summary(x$reg.output[[.(m)]]$yreg)$survey.design)$data,
                                                                                    weights = getCall(summary(x$reg.output[[.(m)]]$yreg)$survey.design)$weights, 
                                                                                    evaluate = FALSE))))
            print(x$reg.output[[m]]$yreg)
          } else {
            x$reg.output[[m]]$yreg$call <- eval(bquote(update(x$reg.output[[.(m)]]$yreg,
                                                              data = getCall(x$reg.output[[.(m)]]$yreg)$data,
                                                              weights = getCall(x$reg.output[[.(m)]]$yreg)$weights, 
                                                              evaluate = FALSE)))
            print(x$reg.output[[m]]$yreg)
          }
        }
        if (name == "yregTot") {
          cat("## Outcome regression for the total effect: \n")
          x$reg.output[[m]]$yregTot$call <- eval(bquote(update(x$reg.output[[.(m)]]$yregTot,
                                                               data=getCall(x$reg.output[[.(m)]]$yregTot)$data,
                                                               weights=getCall(x$reg.output[[.(m)]]$yregTot)$weights, 
                                                               evaluate = FALSE)))
          print(x$reg.output[[m]]$yregTot)
        }
        if (name == "yregDir") {
          cat("## Outcome regression for the direct effect: \n")
          x$reg.output[[m]]$yregDir$call <- eval(bquote(update(x$reg.output[[.(m)]]$yregDir,
                                                               data=getCall(x$reg.output[[.(m)]]$yregDir)$data,
                                                               weights=getCall(x$reg.output[[.(m)]]$yregDir)$weights, 
                                                               evaluate = FALSE)))
          print(x$reg.output[[m]]$yregDir)
        }
        if (name == "ereg") {
          cat("## Exposure regression for weighting: \n")
          x$reg.output[[m]]$ereg$call <- eval(bquote(update(x$reg.output[[.(m)]]$ereg,
                                                            data=getCall(x$reg.output[[.(m)]]$ereg)$data,
                                                            weights=getCall(x$reg.output[[.(m)]]$ereg)$weights, 
                                                            evaluate = FALSE)))
          print(x$reg.output[[m]]$ereg)
        }
        if (name == "mreg") {
          cat("## Mediator regressions: \n")
          for (i in 1:length(x$reg.output[[m]]$mreg)) {
            if (inherits(x$reg.output[[m]]$mreg[[i]], "svyglm")) {
              x$reg.output[[m]]$mreg[[i]]$call <- eval(bquote(update(x$reg.output[[.(m)]]$mreg[[.(i)]],
                                                                     design = getCall(x$reg.output[[.(m)]]$mreg[[.(i)]])$design,
                                                                     family = getCall(x$reg.output[[.(m)]]$mreg[[.(i)]])$family, evaluate = FALSE)))
              x$reg.output[[m]]$mreg[[i]]$survey.design$call <- eval(bquote(as.call(update(summary(x$reg.output[[.(m)]]$mreg[[.(i)]])$survey.design,
                                                                                           data = getCall(summary(x$reg.output[[.(m)]]$mreg[[.(i)]])$survey.design)$data,
                                                                                           weights = getCall(summary(x$reg.output[[.(m)]]$mreg[[.(i)]])$survey.design)$weights, 
                                                                                           evaluate = FALSE))))
              print(x$reg.output[[m]]$mreg[[i]])
            } else {
              x$reg.output[[m]]$mreg[[i]]$call <- eval(bquote(update(x$reg.output[[.(m)]]$mreg[[.(i)]], 
                                                                     data=getCall(x$reg.output[[.(m)]]$mreg[[.(i)]])$data, 
                                                                     weights=getCall(x$reg.output[[.(m)]]$mreg[[.(i)]])$weights, 
                                                                     evaluate = FALSE)))
              print(x$reg.output[[m]]$mreg[[i]])
            }
            if (i < length(x$reg.output[[m]]$mreg)) cat("\n")
          }
        }
        if (name == "wmdenomreg") {
          if (!is.null(x$reg.output[[m]]$wmdenomreg)) {
            cat("## Mediator regressions for weighting (denominator): \n")
            for (i in 1:length(x$reg.output[[m]]$wmdenomreg)) {
              x$reg.output[[m]]$wmdenomreg[[i]]$call <- eval(bquote(update(x$reg.output[[.(m)]]$wmdenomreg[[i]], 
                                                                           data=getCall(x$reg.output[[.(m)]]$wmdenomreg[[.(i)]])$data, 
                                                                           weights=getCall(x$reg.output[[.(m)]]$wmdenomreg[[.(i)]])$weights, 
                                                                           evaluate = FALSE)))
              print(x$reg.output[[m]]$wmdenomreg[[i]])
              if (i < length(x$reg.output[[m]]$wmdenomreg)) cat("\n")
            }
          }
        }
        if (name == "wmnomreg") {
          if (!is.null(x$reg.output[[m]]$wmnomreg)) {
            cat("## Mediator regressions for weighting (nominator): \n")
            for (i in 1:length(x$reg.output[[m]]$wmnomreg)) {
              x$reg.output[[m]]$wmnomreg[[i]]$call <- eval(bquote(update(x$reg.output[[.(m)]]$wmnomreg[[i]], 
                                                                         data=getCall(x$reg.output[[.(m)]]$wmnomreg[[.(i)]])$data, 
                                                                         weights=getCall(x$reg.output[[.(m)]]$wmnomreg[[.(i)]])$weights, 
                                                                         evaluate = FALSE)))
              print(x$reg.output[[m]]$wmnomreg[[i]])
              if (i < length(x$reg.output[[m]]$wmnomreg)) cat("\n")
            }
          }
        }
        if (name == "postcreg") {
          if (!is.null(x$reg.output[[m]]$postcreg)) {
            cat("## Regressions for mediator-outcome confounders affected by the exposure: \n")
            for (i in 1:length(x$reg.output[[m]]$postcreg)) {
              x$reg.output[[m]]$postcreg[[i]]$call <- eval(bquote(update(x$reg.output[[.(m)]]$postcreg[[i]], 
                                                                         data=getCall(x$reg.output[[.(m)]]$postcreg[[.(i)]])$data, 
                                                                         weights=getCall(x$reg.output[[.(m)]]$postcreg[[.(i)]])$weights, 
                                                                         evaluate = FALSE)))
              print(x$reg.output[[m]]$postcreg[[i]])
              if (i < length(x$reg.output[[m]]$postcreg)) cat("\n")
            }
          }
        }
        cat("\n")
      }
    }
  }
  
  # scale and legend
  full <- x$methods$full
  model <- x$methods$model
  EMint <- x$variables$EMint
  if (model == "iorw") {
    if (x$multimp$multimp) yreg_mid <- x$reg.output[[1]]$yregTot
    if (!x$multimp$multimp) yreg_mid <- x$reg.output$yregTot
  } else {
    if (x$multimp$multimp) yreg_mid <- x$reg.output[[1]]$yreg
    if (!x$multimp$multimp) yreg_mid <- x$reg.output$yreg
  }
  if (inherits(yreg_mid, "rcreg") | inherits(yreg_mid, "simexreg")) yreg_mid <- yreg_mid$NAIVEreg
  is_lm <- inherits(yreg_mid, "lm")
  is_glm <- inherits(yreg_mid, "glm")
  is_svyglm <- inherits(yreg_mid, "svyglm")
  is_gam <- inherits(yreg_mid, "gam")
  if (is_lm | is_glm) family_yreg <- family(yreg_mid)
  is_multinom <- inherits(yreg_mid, "multinom")
  is_svymultinom <- inherits(yreg_mid, "svymultinom")
  is_polr <- inherits(yreg_mid, "polr")
  is_survreg <- inherits(yreg_mid, "survreg")
  is_coxph <- inherits(yreg_mid, "coxph")
  if ((is_lm | is_glm) && (family_yreg$family %in% c("gaussian", "inverse.gaussian", "quasi", "Gamma"))) {
    scale <- "mean difference scale"
    if (model == "iorw") {
      if (full) legend <- "(te: total effect; pnde: pure natural direct effect; tnie: total natural indirect effect; pm: proportion mediated)"   
      if (!full) legend <- "(te: total effect; pnde: pure natural direct effect; tnie: total natural indirect effect)"
    } else if (length(x$variables$postc) != 0) {
      if (full) {
        if (EMint) legend <- "(cde: controlled direct effect; rpnde: randomized analogue of pure natural direct effect; rtnde: randomized analogue of total natural direct effect; rpnie: randomized analogue of pure natural indirect effect; rtnie: randomized analogue of total natural indirect effect; te: total effect; rintref: randomized analogue of reference interaction; rintmed: randomized analogue of mediated interaction; cde(prop): proportion cde; rintref(prop): proportion rintref; rintmed(prop): proportion rintmed; rpnie(prop): proportion rpnie; rpm: randomized analogue of overall proportion mediated; rint: randomized analogue of overall proportion attributable to interaction; rpe: randomized analogue of overall proportion eliminated)"
        if (!EMint) legend <- "(cde: controlled direct effect; rpnde: randomized analogue of pure natural direct effect; rtnde: randomized analogue of total natural direct effect; rpnie: randomized analogue of pure natural indirect effect; rtnie: randomized analogue of total natural indirect effect; te: total effect; rpm: randomized analogue of overall proportion mediated)"
      } else legend <- "(cde: controlled direct effect; rpnde: randomized analogue of pure natural direct effect; rtnde: randomized analogue of total natural direct effect; rpnie: randomized analogue of pure natural indirect effect; rtnie: randomized analogue of total natural indirect effect; te: total effect)"
    } else {
      if (full) {
        if (EMint) legend <- "(cde: controlled direct effect; pnde: pure natural direct effect; tnde: total natural direct effect; pnie: pure natural indirect effect; tnie: total natural indirect effect; te: total effect; intref: reference interaction; intmed: mediated interaction; cde(prop): proportion cde; intref(prop): proportion intref; intmed(prop): proportion intmed; pnie(prop): proportion pnie; pm: overall proportion mediated; int: overall proportion attributable to interaction; pe: overall proportion eliminated)"
        if (!EMint) legend <- "(cde: controlled direct effect; pnde: pure natural direct effect; tnde: total natural direct effect; pnie: pure natural indirect effect; tnie: total natural indirect effect; te: total effect; pm: overall proportion mediated)"
      } else legend <- "(cde: controlled direct effect; pnde: pure natural direct effect; tnde: total natural direct effect; pnie: pure natural indirect effect; tnie: total natural indirect effect; te: total effect)"
    }
  } else if ((is_lm | is_glm) && (family_yreg$family %in% c("poisson", "quasipoisson", "ziplss") |
                                  startsWith(family_yreg$family, "Negative Binomial") |
                                  startsWith(family_yreg$family, "Zero inflated Poisson"))) {
    scale <- "rate ratio scale"
    if (model == "iorw") {
      if (full) legend <- "(Rte: total effect rate ratio; Rpnde: pure natural direct effect rate ratio; Rtnie: total natural indirect effect rate ratio; pm: proportion mediated)"   
      if (!full) legend <- "(Rte: total effect rate ratio; Rpnde: pure natural direct effect rate ratio; Rtnie: total natural indirect effect rate ratio)" 
    } else if (length(x$variables$postc) != 0) {
      if (full) {
        if (EMint) legend <- "(Rcde: controlled direct effect rate ratio; rRpnde: randomized analogue of pure natural direct effect rate ratio; rRtnde: randomized analogue of total natural direct effect rate ratio; rRpnie: randomized analogue of pure natural indirect effect rate ratio; rRtnie: randomized analogue of total natural indirect effect rate ratio; Rte: total effect rate ratio; ERcde: excess relative rate due to controlled direct effect; rERintref: randomized analogue of excess relative rate due to reference interaction; rERintmed: randomized analogue of excess relative rate due to mediated interaction; rERpnie: randomized analogue of excess relative rate due to pure natural indirect effect; ERcde(prop): proportion ERcde; rERintref(prop): proportion rERintref; rERintmed(prop): proportion rERintmed; rERpnie(prop): proportion rERpnie; rpm: randomized analogue of overall proportion mediated; rint: randomized analogue of overall proportion attributable to interaction; rpe: randomized analogue of overall proportion eliminated)"
        if (!EMint) legend <- "(Rcde: controlled direct effect rate ratio; rRpnde: randomized analogue of pure natural direct effect rate ratio; rRtnde: randomized analogue of total natural direct effect rate ratio; rRpnie: randomized analogue of pure natural indirect effect rate ratio; rRtnie: randomized analogue of total natural indirect effect rate ratio; Rte: total effect rate ratio; rpm: randomized analogue of overall proportion mediated)"
      } else legend <- "(Rcde: controlled direct effect rate ratio; rRpnde: randomized analogue of pure natural direct effect rate ratio; rRtnde: randomized analogue of total natural direct effect rate ratio; rRpnie: randomized analogue of pure natural indirect effect rate ratio; rRtnie: randomized analogue of total natural indirect effect rate ratio; Rte: total effect rate ratio)"
    } else {
      if (full) {
        if (EMint) legend <- "(Rcde: controlled direct effect rate ratio; Rpnde: pure natural direct effect rate ratio; Rtnde: total natural direct effect rate ratio; Rpnie: pure natural indirect effect rate ratio; Rtnie: total natural indirect effect rate ratio; Rte: total effect rate ratio; ERcde: excess relative rate due to controlled direct effect; ERintref: excess relative rate due to reference interaction; ERintmed: excess relative rate due to mediated interaction; ERpnie: excess relative rate due to pure natural indirect effect; ERcde(prop): proportion ERcde; ERintref(prop): proportion ERintref; ERintmed(prop): proportion ERintmed; ERpnie(prop): proportion ERpnie; pm: overall proportion mediated; int: overall proportion attributable to interaction; pe: overall proportion eliminated)"
        if (!EMint) legend <- "(Rcde: controlled direct effect rate ratio; Rpnde: pure natural direct effect rate ratio; Rtnde: total natural direct effect rate ratio; Rpnie: pure natural indirect effect rate ratio; Rtnie: total natural indirect effect rate ratio; Rte: total effect rate ratio; pm: overall proportion mediated)"
      } else legend <- "(Rcde: controlled direct effect rate ratio; Rpnde: pure natural direct effect rate ratio; Rtnde: total natural direct effect rate ratio; Rpnie: pure natural indirect effect rate ratio; Rtnie: total natural indirect effect rate ratio; Rte: total effect rate ratio)"
    }
  } else if (((is_lm | is_glm) && (family_yreg$family %in% c("binomial", "quasibinomial", "multinom") |
                                   startsWith(family_yreg$family, "Ordered Categorical"))) |
             is_multinom | is_polr) {
    
    if (is_glm && family_yreg$family %in% c("binomial", "quasibinomial") && yreg_mid$family$link == "logit") {
      scale <- "odds ratio scale"
      if (model == "iorw") {
        if (full) legend <- "(Rte: total effect odds ratio; Rpnde: pure natural direct effect odds ratio; Rtnie: total natural indirect effect odds ratio; pm: proportion mediated)"   
        if (!full) legend <- "(Rte: total effect odds ratio; Rpnde: pure natural direct effect odds ratio; Rtnie: total natural indirect effect odds ratio)" 
      } else if (length(x$variables$postc) != 0) {
        if (full) {
          if (EMint) legend <- "(Rcde: controlled direct effect odds ratio; rRpnde: randomized analogue of pure natural direct effect odds ratio; rRtnde: randomized analogue of total natural direct effect odds ratio; rRpnie: randomized analogue of pure natural indirect effect odds ratio; rRtnie: randomized analogue of total natural indirect effect odds ratio; Rte: total effect odds ratio; ERcde: excess relative risk due to controlled direct effect; rERintref: randomized analogue of excess relative risk due to reference interaction; rERintmed: randomized analogue of excess relative risk due to mediated interaction; rERpnie: randomized analogue of excess relative risk due to pure natural indirect effect; ERcde(prop): proportion ERcde; rERintref(prop): proportion rERintref; rERintmed(prop): proportion rERintmed; rERpnie(prop): proportion rERpnie; rpm: randomized analogue of overall proportion mediated; rint: randomized analogue of overall proportion attributable to interaction; rpe: randomized analogue of overall proportion eliminated)"
          if (!EMint) legend <- "(Rcde: controlled direct effect odds ratio; rRpnde: randomized analogue of pure natural direct effect odds ratio; rRtnde: randomized analogue of total natural direct effect odds ratio; rRpnie: randomized analogue of pure natural indirect effect odds ratio; rRtnie: randomized analogue of total natural indirect effect odds ratio; Rte: total effect odds ratio; rpm: randomized analogue of overall proportion mediated)"
        } else legend <- "(Rcde: controlled direct effect odds ratio; rRpnde: randomized analogue of pure natural direct effect odds ratio; rRtnde: randomized analogue of total natural direct effect odds ratio; rRpnie: randomized analogue of pure natural indirect effect odds ratio; rRtnie: randomized analogue of total natural indirect effect odds ratio; Rte: total effect odds ratio)"
      } else {
        if (full) {
          if (EMint) legend <- "(Rcde: controlled direct effect odds ratio; Rpnde: pure natural direct effect odds ratio; Rtnde: total natural direct effect odds ratio; Rpnie: pure natural indirect effect odds ratio; Rtnie: total natural indirect effect odds ratio; Rte: total effect odds ratio; ERcde: excess relative risk due to controlled direct effect; ERintref: excess relative risk due to reference interaction; ERintmed: excess relative risk due to mediated interaction; ERpnie: excess relative risk due to pure natural indirect effect; ERcde(prop): proportion ERcde; ERintref(prop): proportion ERintref; ERintmed(prop): proportion ERintmed; ERpnie(prop): proportion ERpnie; pm: overall proportion mediated; int: overall proportion attributable to interaction; pe: overall proportion eliminated)"
          if (!EMint) legend <- "(Rcde: controlled direct effect odds ratio; Rpnde: pure natural direct effect odds ratio; Rtnde: total natural direct effect odds ratio; Rpnie: pure natural indirect effect odds ratio; Rtnie: total natural indirect effect odds ratio; Rte: total effect odds ratio; pm: overall proportion mediated)"
        } else legend <- "(Rcde: controlled direct effect odds ratio; Rpnde: pure natural direct effect odds ratio; Rtnde: total natural direct effect odds ratio; Rpnie: pure natural indirect effect odds ratio; Rtnie: total natural indirect effect odds ratio; Rte: total effect odds ratio)"
      }
    } else {
      scale <- "risk ratio scale"
      if (model == "iorw") {
        if (full) legend <- "(Rte: total effect risk ratio; Rpnde: pure natural direct effect risk ratio; Rtnie: total natural indirect effect risk ratio; pm: proportion mediated)"   
        if (!full) legend <- "(Rte: total effect risk ratio; Rpnde: pure natural direct effect risk ratio; Rtnie: total natural indirect effect risk ratio)" 
      } else if (length(x$variables$postc) != 0) {
        if (full) {
          if (EMint) legend <- "(Rcde: controlled direct effect risk ratio; rRpnde: randomized analogue of pure natural direct effect risk ratio; rRtnde: randomized analogue of total natural direct effect risk ratio; rRpnie: randomized analogue of pure natural indirect effect risk ratio; rRtnie: randomized analogue of total natural indirect effect risk ratio; Rte: total effect risk ratio; ERcde: excess relative risk due to controlled direct effect; rERintref: randomized analogue of excess relative risk due to reference interaction; rERintmed: randomized analogue of excess relative risk due to mediated interaction; rERpnie: randomized analogue of excess relative risk due to pure natural indirect effect; ERcde(prop): proportion ERcde; rERintref(prop): proportion rERintref; rERintmed(prop): proportion rERintmed; rERpnie(prop): proportion rERpnie; rpm: randomized analogue of overall proportion mediated; rint: randomized analogue of overall proportion attributable to interaction; rpe: randomized analogue of overall proportion eliminated)"
          if (!EMint) legend <- "(Rcde: controlled direct effect risk ratio; rRpnde: randomized analogue of pure natural direct effect risk ratio; rRtnde: randomized analogue of total natural direct effect risk ratio; rRpnie: randomized analogue of pure natural indirect effect risk ratio; rRtnie: randomized analogue of total natural indirect effect risk ratio; Rte: total effect risk ratio; rpm: randomized analogue of overall proportion mediated)"
        } else legend <- "(Rcde: controlled direct effect risk ratio; rRpnde: randomized analogue of pure natural direct effect risk ratio; rRtnde: randomized analogue of total natural direct effect risk ratio; rRpnie: randomized analogue of pure natural indirect effect risk ratio; rRtnie: randomized analogue of total natural indirect effect risk ratio; Rte: total effect risk ratio)"
      } else {
        if (full) {
          if (EMint) legend <- "(Rcde: controlled direct effect risk ratio; Rpnde: pure natural direct effect risk ratio; Rtnde: total natural direct effect risk ratio; Rpnie: pure natural indirect effect risk ratio; Rtnie: total natural indirect effect risk ratio; Rte: total effect risk ratio; ERcde: excess relative risk due to controlled direct effect; ERintref: excess relative risk due to reference interaction; ERintmed: excess relative risk due to mediated interaction; ERpnie: excess relative risk due to pure natural indirect effect; ERcde(prop): proportion ERcde; ERintref(prop): proportion ERintref; ERintmed(prop): proportion ERintmed; ERpnie(prop): proportion ERpnie; pm: overall proportion mediated; int: overall proportion attributable to interaction; pe: overall proportion eliminated)"
          if (!EMint) legend <- "(Rcde: controlled direct effect risk ratio; Rpnde: pure natural direct effect risk ratio; Rtnde: total natural direct effect risk ratio; Rpnie: pure natural indirect effect risk ratio; Rtnie: total natural indirect effect risk ratio; Rte: total effect risk ratio; pm: overall proportion mediated)"
        } else legend <- "(Rcde: controlled direct effect risk ratio; Rpnde: pure natural direct effect risk ratio; Rtnde: total natural direct effect risk ratio; Rpnie: pure natural indirect effect risk ratio; Rtnie: total natural indirect effect risk ratio; Rte: total effect risk ratio)"
      }
    }
    
  } else if (is_coxph) {
    scale <- "hazard ratio scale"
    if (model == "iorw") {
      if (full) legend <- "(Rte: total effect hazard ratio; Rpnde: pure natural direct effect hazard ratio; Rtnie: total natural indirect effect hazard ratio; pm: proportion mediated)"   
      if (!full) legend <- "(Rte: total effect hazard ratio; Rpnde: pure natural direct effect hazard ratio; Rtnie: total natural indirect effect hazard ratio)" 
    }  else if (length(x$variables$postc) != 0) {
      if (full) {
        if (EMint) legend <- "(Rcde: controlled direct effect hazard ratio; rRpnde: randomized analogue of pure natural direct effect hazard ratio; rRtnde: randomized analogue of total natural direct effect hazard ratio; rRpnie: randomized analogue of pure natural indirect effect hazard ratio; rRtnie: randomized analogue of total natural indirect effect hazard ratio; Rte: total effect hazard ratio; ERcde: excess relative hazard due to controlled direct effect; rERintref: randomized analogue of excess relative hazard due to reference interaction; rERintmed: randomized analogue of excess relative hazard due to mediated interaction; rERpnie: randomized analogue of excess relative hazard due to pure natural indirect effect; ERcde(prop): proportion ERcde; rERintref(prop): proportion rERintref; rERintmed(prop): proportion rERintmed; rERpnie(prop): proportion rERpnie; rpm: randomized analogue of overall proportion mediated; rint: randomized analogue of overall proportion attributable to interaction; rpe: randomized analogue of overall proportion eliminated)"
        if (!EMint) legend <- "(Rcde: controlled direct effect hazard ratio; rRpnde: randomized analogue of pure natural direct effect hazard ratio; rRtnde: randomized analogue of total natural direct effect hazard ratio; rRpnie: randomized analogue of pure natural indirect effect hazard ratio; rRtnie: randomized analogue of total natural indirect effect hazard ratio; Rte: total effect hazard ratio; rpm: randomized analogue of overall proportion mediated)"
      } else legend <- "(Rcde: controlled direct effect hazard ratio; rRpnde: randomized analogue of pure natural direct effect hazard ratio; rRtnde: randomized analogue of total natural direct effect hazard ratio; rRpnie: randomized analogue of pure natural indirect effect hazard ratio; rRtnie: randomized analogue of total natural indirect effect hazard ratio; Rte: total effect hazard ratio)"
    } else {
      if (full) {
        if (EMint) legend <- "(Rcde: controlled direct effect hazard ratio; Rpnde: pure natural direct effect hazard ratio; Rtnde: total natural direct effect hazard ratio; Rpnie: pure natural indirect effect hazard ratio; Rtnie: total natural indirect effect hazard ratio; Rte: total effect hazard ratio; ERcde: excess relative hazard due to controlled direct effect; ERintref: excess relative hazard due to reference interaction; ERintmed: excess relative hazard due to mediated interaction; ERpnie: excess relative hazard due to pure natural indirect effect; ERcde(prop): proportion ERcde; ERintref(prop): proportion ERintref; ERintmed(prop): proportion ERintmed; ERpnie(prop): proportion ERpnie; pm: overall proportion mediated; int: overall proportion attributable to interaction; pe: overall proportion eliminated)"
        if (!EMint) legend <- "(Rcde: controlled direct effect hazard ratio; Rpnde: pure natural direct effect hazard ratio; Rtnde: total natural direct effect hazard ratio; Rpnie: pure natural indirect effect hazard ratio; Rtnie: total natural indirect effect hazard ratio; Rte: total effect hazard ratio; pm: overall proportion mediated)"
      } else legend <- "(Rcde: controlled direct effect hazard ratio; Rpnde: pure natural direct effect hazard ratio; Rtnde: total natural direct effect hazard ratio; Rpnie: pure natural indirect effect hazard ratio; Rtnie: total natural indirect effect hazard ratio; Rte: total effect hazard ratio)"
    }
  } else if (is_survreg) {
    scale <- "mean survival scale"
    if (model == "iorw") {
      if (full) legend <- "(Rte: total effect mean survival ratio; Rpnde: pure natural direct effect mean survival ratio; Rtnie: total natural indirect effect mean survival ratio; pm: proportion mediated)"   
      if (!full) legend <- "(Rte: total effect mean survival ratio; Rpnde: pure natural direct effect mean survival ratio; Rtnie: total natural indirect effect mean survival ratio)" 
    } else if (length(x$variables$postc) != 0) {
      if (full) {
        if (EMint) legend <- "(Rcde: controlled direct effect mean survival ratio; rRpnde: randomized analogue of pure natural direct effect mean survival ratio; rRtnde: randomized analogue of total natural direct effect mean survival ratio; rRpnie: randomized analogue of pure natural indirect effect mean survival ratio; rRtnie: randomized analogue of total natural indirect effect mean survival ratio; Rte: total effect mean survival ratio; ERcde: excess mean survival ratio due to controlled direct effect; rERintref: randomized analogue of excess mean survival ratio due to reference interaction; rERintmed: randomized analogue of excess mean survival ratio due to mediated interaction; rERpnie: randomized analogue of excess mean survival ratio due to pure natural indirect effect; ERcde(prop): proportion ERcde; rERintref(prop): proportion rERintref; rERintmed(prop): proportion rERintmed; rERpnie(prop): proportion rERpnie; rpm: randomized analogue of overall proportion mediated; rint: randomized analogue of overall proportion attributable to interaction; rpe: randomized analogue of overall proportion eliminated)"
        if (!EMint) legend <- "(Rcde: controlled direct effect mean survival ratio; rRpnde: randomized analogue of pure natural direct effect mean survival ratio; rRtnde: randomized analogue of total natural direct effect mean survival ratio; rRpnie: randomized analogue of pure natural indirect effect mean survival ratio; rRtnie: randomized analogue of total natural indirect effect mean survival ratio; Rte: total effect mean survival ratio; rpm: randomized analogue of overall proportion mediated)"
      } else legend <- "(Rcde: controlled direct effect mean survival ratio; rRpnde: randomized analogue of pure natural direct effect mean survival ratio; rRtnde: randomized analogue of total natural direct effect mean survival ratio; rRpnie: randomized analogue of pure natural indirect effect mean survival ratio; rRtnie: randomized analogue of total natural indirect effect mean survival ratio; Rte: total effect mean survival ratio)"
    } else {
      if (full) {
        if (EMint) legend <- "(Rcde: controlled direct effect mean survival ratio; Rpnde: pure natural direct effect mean survival ratio; Rtnde: total natural direct effect mean survival ratio; Rpnie: pure natural indirect effect mean survival ratio; Rtnie: total natural indirect effect mean survival ratio; Rte: total effect mean survival ratio; ERcde: excess mean survival ratio due to controlled direct effect; ERintref: excess mean survival ratio due to reference interaction; ERintmed: excess mean survival ratio due to mediated interaction; ERpnie: excess mean survival ratio due to pure natural indirect effect; ERcde(prop): proportion ERcde; ERintref(prop): proportion ERintref; ERintmed(prop): proportion ERintmed; ERpnie(prop): proportion ERpnie; pm: overall proportion mediated; int: overall proportion attributable to interaction; pe: overall proportion eliminated)"
        if (!EMint) legend <- "(Rcde: controlled direct effect mean survival ratio; Rpnde: pure natural direct effect mean survival ratio; Rtnde: total natural direct effect mean survival ratio; Rpnie: pure natural indirect effect mean survival ratio; Rtnie: total natural indirect effect mean survival ratio; Rte: total effect mean survival ratio; pm: overall proportion mediated)"
      } else legend <- "(Rcde: controlled direct effect mean survival ratio; Rpnde: pure natural direct effect mean survival ratio; Rtnde: total natural direct effect mean survival ratio; Rpnie: pure natural indirect effect mean survival ratio; Rtnie: total natural indirect effect mean survival ratio; Rte: total effect mean survival ratio)"
    }
  }
  
  # print causal mediation analysis results
  if (x$methods$model == "rb") model_str <- "regression-based approach"
  if (x$methods$model == "wb") model_str <- "weighting-based approach"
  if (x$methods$model == "ne") model_str <- "natural effect model"
  if (x$methods$model == "iorw") model_str <- "inverse odds ratio weighting approach"
  if (x$methods$model == "msm") model_str <- "marginal structural model"
  if (x$methods$model == "gformula") model_str <- "g-formula approach"
  if (x$multimp$multimp) model_str <- paste(model_str, "with multiple imputation")
  if (x$methods$estimation == "paramfunc") est_str <- "Closed-form parameter function estimation"
  if (x$methods$estimation == "imputation") est_str <- "Direct counterfactual imputation estimation"
  if (x$methods$inference == "delta") inf_str <- "delta method standard errors, confidence intervals and p-values"
  if (x$methods$inference == "bootstrap") {
    if (x$methods$boot.ci.type == "per") inf_str <- "bootstrap standard errors, percentile confidence intervals and p-values"
    if (x$methods$boot.ci.type == "bca") inf_str <- "bootstrap standard errors, bias-corrected and accelerated confidence intervals and p-values"
  }
  if (x$methods$casecontrol) cat(paste("# Effect decomposition on the", scale, "for a case control study via the "))
  if (!(x$methods$casecontrol)) cat(paste("# Effect decomposition on the", scale, "via the "))
  cat(model_str)
  cat("\n \n")
  cat(est_str)
  cat(paste(" with \n", inf_str, "\n \n"))
  print(x$effect.pe)
  cat("\n")
  cat(legend)
  cat("\n\nRelevant variable values: \n")
  print(x$ref)
}


#' @describeIn cmest Summarize the results of \code{cmest} nicely
#' @export
summary.cmest <- function(object, ...) {
  # summarize regressions
  # print regression models used
  out <- object
  out$reg.output.summary <- out$reg.output
  if (!object$multimp$multimp) {
    regnames <- names(object$reg.output)
    for (name in regnames) {
      if (name == "yreg") out$reg.output.summary$yreg <- summary(object$reg.output$yreg)
      if (name == "yregTot") out$reg.output.summary$yregTot <- summary(object$reg.output$yregTot)
      if (name == "yregDir") out$reg.output.summary$yregDir <- summary(object$reg.output$yregDir)
      if (name == "ereg") out$reg.output.summary$ereg <- summary(object$reg.output$ereg)
      if (name == "mreg") out$reg.output.summary$mreg <- lapply(1:length(object$reg.output$mreg), function(i) 
        summary(object$reg.output$mreg[[i]]))
      if (name == "wmnomreg") out$reg.output.summary$wmnomreg <- lapply(1:length(object$reg.output$wmnomreg), function(i) 
        summary(object$reg.output$wmnomreg[[i]]))
      if (name == "wmdenomreg") out$reg.output.summary$wmdenomreg <- lapply(1:length(object$reg.output$wmdenomreg), function(i) 
        summary(object$reg.output$wmdenomreg[[i]]))
      if (name == "postcreg") out$reg.output.summary$postcreg <- lapply(1:length(object$reg.output$postcreg), function(i) 
        summary(object$reg.output$postcreg[[i]]))
    }
  } else {
    for (m in 1:length(object$reg.output)){ 
      regnames <- names(object$reg.output[[m]])
      for (name in regnames) {
        if (name == "yreg") out$reg.output.summary[[m]]$yreg <- summary(object$reg.output[[m]]$yreg)
        if (name == "yregTot") out$reg.output.summary[[m]]$yregTot <- summary(object$reg.output[[m]]$yregTot)
        if (name == "yregDir") out$reg.output.summary[[m]]$yregDir <- summary(object$reg.output[[m]]$yregDir)
        if (name == "ereg") out$reg.output.summary[[m]]$ereg <- summary(object$reg.output[[m]]$ereg)
        if (name == "mreg") out$reg.output.summary[[m]]$mreg <- lapply(1:length(object$reg.output[[m]]$mreg), function(i) 
          summary(object$reg.output[[m]]$mreg[[i]]))
        if (name == "wmnomreg") out$reg.output.summary[[m]]$wmnomreg <- lapply(1:length(object$reg.output[[m]]$wmnomreg), function(i) 
          summary(object$reg.output[[m]]$wmnomreg[[i]]))
        if (name == "wmdenomreg") out$reg.output.summary[[m]]$wmdenomreg <- lapply(1:length(object$reg.output[[m]]$wmdenomreg), function(i) 
          summary(object$reg.output[[m]]$wmdenomreg[[i]]))
        if (name == "postcreg") out$reg.output.summary[[m]]$postcreg <- lapply(1:length(object$reg.output[[m]]$postcreg), function(i) 
          summary(object$reg.output[[m]]$postcreg[[i]]))
      }
    }
  }
  
  # summarize causal mediation analysis results
  summarydf <- data.frame(object$effect.pe, object$effect.se, object$effect.ci.low, 
                          object$effect.ci.high, object$effect.pval)
  colnames(summarydf) <- c("Estimate", "Std.error", "95% CIL", "95% CIU", "P.val")
  out$summarydf <- summarydf
  class(out) <- c("summary.cmest")
  return(out)
}


#' @describeIn cmest Print the summary of \code{cmest} nicely
#' @export
print.summary.cmest <- function(x, digits = 4, ...) {
  cat("Causal Mediation Analysis\n\n")
  # print summary of regression models used
  if (!x$multimp$multimp) {
    regnames <- names(x$reg.output.summary)
    for (name in regnames) {
      if (name == "yreg") {
        cat("# Outcome regression:\n")
        if (inherits(x$reg.output$yreg, "svyglm")) {
          x$reg.output.summary$yreg$call <- update(x$reg.output$yreg,design = getCall(x$reg.output$yreg)$design,
                                                   family = getCall(x$reg.output$yreg)$family, evaluate = FALSE)
          x$reg.output.summary$yreg$survey.design$call <- as.call(update(x$reg.output.summary$yreg$survey.design,
                                                                         data = getCall(x$reg.output.summary$yreg$survey.design)$data,
                                                                         weights = getCall(x$reg.output.summary$yreg$survey.design)$weights, 
                                                                         evaluate = FALSE)) 
          print(x$reg.output.summary$yreg)
        } else {
          x$reg.output.summary$yreg$call <- update(x$reg.output$yreg,data=getCall(x$reg.output$yreg)$data,
                                                   weights=getCall(x$reg.output$yreg)$weights, evaluate = FALSE)
          print(x$reg.output.summary$yreg)
        }
      }
      if (name == "yregTot") {
        cat("# Outcome regression for the total effect: \n")
        x$reg.output.summary$yregTot$call <- update(x$reg.output$yregTot,data=getCall(x$reg.output$yregTot)$data,
                                                    weights=getCall(x$reg.output$yregTot)$weights, evaluate = FALSE)
        print(x$reg.output.summary$yregTot)
      }
      if (name == "yregDir") {
        cat("# Outcome regression for the direct effect: \n")
        x$reg.output.summary$yregDir$call <- update(x$reg.output$yregDir,data=getCall(x$reg.output$yregDir)$data,
                                                    weights=getCall(x$reg.output$yregDir)$weights, evaluate = FALSE)
        print(x$reg.output.summary$yregDir)
      }
      if (name == "ereg") {
        cat("# Exposure regression for weighting: \n")
        x$reg.output.summary$ereg$call <- update(x$reg.output$ereg,data=getCall(x$reg.output$ereg)$data,
                                                 weights=getCall(x$reg.output$ereg)$weights, evaluate = FALSE)
        print(x$reg.output.summary$ereg)
      }
      if (name == "mreg") {
        cat("# Mediator regressions: \n")
        for (i in 1:length(x$reg.output.summary$mreg)) {
          if (inherits(x$reg.output$mreg[[i]], "svyglm")) {
            x$reg.output.summary$mreg[[i]]$call <- eval(bquote(update(x$reg.output$mreg[[.(i)]],design = getCall(x$reg.output$mreg[[.(i)]])$design,
                                                                      family = getCall(x$reg.output$mreg[[.(i)]])$family, evaluate = FALSE)))
            x$reg.output.summary$mreg[[i]]$survey.design$call <- eval(bquote(as.call(update(x$reg.output.summary$mreg[[.(i)]]$survey.design,
                                                                                            data = getCall(x$reg.output.summary$mreg[[.(i)]]$survey.design)$data,
                                                                                            weights = getCall(x$reg.output.summary$mreg[[.(i)]]$survey.design)$weights, 
                                                                                            evaluate = FALSE))))
            print(x$reg.output.summary$mreg[[i]])
          } else {
            x$reg.output.summary$mreg[[i]]$call <- eval(bquote(update(x$reg.output$mreg[[.(i)]], 
                                                                      data=getCall(x$reg.output$mreg[[.(i)]])$data, 
                                                                      weights=getCall(x$reg.output$mreg[[.(i)]])$weights,
                                                                      evaluate = FALSE)))
            print(x$reg.output.summary$mreg[[i]])
          }
          if (i < length(x$reg.output$mreg)) cat("\n")
        }
      }
      if (name == "wmdenomreg") {
        cat("# Mediator regressions for weighting (denominator): \n")
        for (i in 1:length(x$reg.output.summary$wmdenomreg)) {
          x$reg.output.summary$wmdenomreg[[i]]$call <- eval(bquote(update(x$reg.output$wmdenomreg[[i]], 
                                                                          data=getCall(x$reg.output$wmdenomreg[[.(i)]])$data, 
                                                                          weights=getCall(x$reg.output$wmdenomreg[[.(i)]])$weights,
                                                                          evaluate = FALSE)))
          print(x$reg.output.summary$wmdenomreg[[i]])
          if (i < length(x$reg.output.summary$wmdenomreg)) cat("\n")
        }
      }
      if (name == "wmnomreg") {
        cat("# Mediator regressions for weighting (nominator): \n")
        for (i in 1:length(x$reg.output.summary$wmnomreg)) {
          x$reg.output.summary$wmnomreg[[i]]$call <- eval(bquote(update(x$reg.output$wmnomreg[[i]], 
                                                                        data=getCall(x$reg.output$wmnomreg[[.(i)]])$data, 
                                                                        weights=getCall(x$reg.output$wmnomreg[[.(i)]])$weights,
                                                                        evaluate = FALSE)))
          print(x$reg.output.summary$wmnomreg[[i]])
          if (i < length(x$reg.output.summary$wmnomreg)) cat("\n")
        }
      }
      if (name == "postcreg") {
        cat("# Regressions for mediator-outcome confounders affected by the exposure: \n")
        for (i in 1:length(x$reg.output.summary$postcreg)) {
          x$reg.output.summary$postcreg[[i]]$call <- eval(bquote(update(x$reg.output$postcreg[[i]], 
                                                                        data=getCall(x$reg.output$postcreg[[.(i)]])$data, 
                                                                        weights=getCall(x$reg.output$postcreg[[.(i)]])$weights,
                                                                        evaluate = FALSE)))
          print(x$reg.output.summary$postcreg[[i]])
          if (i < length(x$reg.output.summary$postcreg)) cat("\n")
        }
      }
      cat("\n")
    }
  } else {
    for (m in 1:length(x$reg.output.summary)){ 
      cat(paste("# Regressions with imputed dataset", m, "\n\n"))
      regnames <- names(x$reg.output.summary[[m]])
      for (name in regnames) {
        if (name == "yreg") {
          cat("## Outcome regression: \n")
          if (inherits(x$reg.output[[m]]$yreg, "svyglm")) {
            x$reg.output.summary[[m]]$yreg$call <- eval(bquote(update(x$reg.output[[.(m)]]$yreg,design = getCall(x$reg.output[[.(m)]]$yreg)$design,
                                                                      family = getCall(x$reg.output[[.(m)]]$yreg)$family, evaluate = FALSE)))
            x$reg.output.summary[[m]]$yreg$survey.design$call <- eval(bquote(as.call(update(x$reg.output.summary[[.(m)]]$yreg$survey.design,
                                                                                            data = getCall(x$reg.output.summary[[.(m)]]$yreg$survey.design)$data,
                                                                                            weights = getCall(x$reg.output.summary[[.(m)]]$yreg$survey.design)$weights, 
                                                                                            evaluate = FALSE))))
            print(x$reg.output.summary[[m]]$yreg)
          } else {
            x$reg.output.summary[[m]]$yreg$call <- eval(bquote(update(x$reg.output[[.(m)]]$yreg,
                                                                      data = getCall(x$reg.output[[.(m)]]$yreg)$data,
                                                                      weights = getCall(x$reg.output[[.(m)]]$yreg)$weights, 
                                                                      evaluate = FALSE)))
            print(x$reg.output.summary[[m]]$yreg)
          }
        }
        if (name == "yregTot") {
          cat("## Outcome regression for the total effect: \n")
          x$reg.output.summary[[m]]$yregTot$call <- eval(bquote(update(x$reg.output[[.(m)]]$yregTot,
                                                                       data=getCall(x$reg.output[[.(m)]]$yregTot)$data,
                                                                       weights=getCall(x$reg.output[[.(m)]]$yregTot)$weights, 
                                                                       evaluate = FALSE)))
          print(x$reg.output.summary[[m]]$yregTot)
        }
        if (name == "yregDir") {
          cat("## Outcome regression for the direct effect: \n")
          x$reg.output.summary[[m]]$yregDir$call <- eval(bquote(update(x$reg.output[[.(m)]]$yregDir,
                                                                       data=getCall(x$reg.output[[.(m)]]$yregDir)$data,
                                                                       weights=getCall(x$reg.output[[.(m)]]$yregDir)$weights, 
                                                                       evaluate = FALSE)))
          print(x$reg.output.summary[[m]]$yregDir)
        }
        if (name == "ereg") {
          cat("## Exposure regression for weighting: \n")
          x$reg.output.summary[[m]]$ereg$call <- eval(bquote(update(x$reg.output[[.(m)]]$ereg,
                                                                    data=getCall(x$reg.output[[.(m)]]$ereg)$data,
                                                                    weights=getCall(x$reg.output[[.(m)]]$ereg)$weights, 
                                                                    evaluate = FALSE)))
          print(x$reg.output.summary[[m]]$ereg)
        }
        if (name == "mreg") {
          cat("## Mediator regressions: \n")
          for (i in 1:length(x$reg.output.summary[[m]]$mreg)) {
            if (inherits(x$reg.output[[m]]$mreg[[i]], "svyglm")) {
              x$reg.output.summary[[m]]$mreg[[i]]$call <- eval(bquote(update(x$reg.output[[.(m)]]$mreg[[.(i)]],
                                                                             design = getCall(x$reg.output[[.(m)]]$mreg[[.(i)]])$design,
                                                                             family = getCall(x$reg.output[[.(m)]]$mreg[[.(i)]])$family, evaluate = FALSE)))
              x$reg.output.summary[[m]]$mreg[[i]]$survey.design$call <- eval(bquote(as.call(update(x$reg.output.summary[[.(m)]]$mreg[[.(i)]]$survey.design,
                                                                                                   data = getCall(x$reg.output.summary[[.(m)]]$mreg[[.(i)]]$survey.design)$data,
                                                                                                   weights = getCall(x$reg.output.summary[[.(m)]]$mreg[[.(i)]]$survey.design)$weights, 
                                                                                                   evaluate = FALSE))))
              print(x$reg.output.summary[[m]]$mreg[[i]])
            } else {
              x$reg.output.summary[[m]]$mreg[[i]]$call <- eval(bquote(update(x$reg.output[[.(m)]]$mreg[[.(i)]], 
                                                                             data=getCall(x$reg.output[[.(m)]]$mreg[[.(i)]])$data, 
                                                                             weights=getCall(x$reg.output[[.(m)]]$mreg[[.(i)]])$weights, 
                                                                             evaluate = FALSE)))
              print(x$reg.output.summary[[m]]$mreg[[i]])
            }
            if (i < length(x$reg.output.summary[[m]]$mreg)) cat("\n")
          }
        }
        if (name == "wmdenomreg") {
          if (!is.null(x$reg.output.summary[[m]]$wmdenomreg)) {
            cat("## Mediator regressions for weighting (denominator): \n")
            for (i in 1:length(x$reg.output.summary[[m]]$wmdenomreg)) {
              x$reg.output.summary[[m]]$wmdenomreg[[i]]$call <- eval(bquote(update(x$reg.output[[.(m)]]$wmdenomreg[[i]], 
                                                                                   data=getCall(x$reg.output[[.(m)]]$wmdenomreg[[.(i)]])$data, 
                                                                                   weights=getCall(x$reg.output[[.(m)]]$wmdenomreg[[.(i)]])$weights, 
                                                                                   evaluate = FALSE)))
              print(x$reg.output.summary[[m]]$wmdenomreg[[i]])
              if (i < length(x$reg.output.summary[[m]]$wmdenomreg)) cat("\n")
            }
          }
        }
        if (name == "wmnomreg") {
          if (!is.null(x$reg.output.summary[[m]]$wmnomreg)) {
            cat("## Mediator regressions for weighting (nominator): \n")
            for (i in 1:length(x$reg.output.summary[[m]]$wmnomreg)) {
              x$reg.output.summary[[m]]$wmnomreg[[i]]$call <- eval(bquote(update(x$reg.output[[.(m)]]$wmnomreg[[i]], 
                                                                                 data=getCall(x$reg.output[[.(m)]]$wmnomreg[[.(i)]])$data, 
                                                                                 weights=getCall(x$reg.output[[.(m)]]$wmnomreg[[.(i)]])$weights, 
                                                                                 evaluate = FALSE)))
              print(x$reg.output.summary[[m]]$wmnomreg[[i]])
              if (i < length(x$reg.output.summary[[m]]$wmnomreg)) cat("\n")
            }
          }
        }
        if (name == "postcreg") {
          if (!is.null(x$reg.output.summary[[m]]$postcreg)) {
            cat("## Regressions for mediator-outcome confounders affected by the exposure: \n")
            for (i in 1:length(x$reg.output.summary[[m]]$postcreg)) {
              x$reg.output.summary[[m]]$postcreg[[i]]$call <- eval(bquote(update(x$reg.output[[.(m)]]$postcreg[[i]], 
                                                                                 data=getCall(x$reg.output[[.(m)]]$postcreg[[.(i)]])$data, 
                                                                                 weights=getCall(x$reg.output[[.(m)]]$postcreg[[.(i)]])$weights, 
                                                                                 evaluate = FALSE)))
              print(x$reg.output.summary[[m]]$postcreg[[i]])
              if (i < length(x$reg.output.summary[[m]]$postcreg)) cat("\n")
            }
          }
        }
        cat("\n")
      }
    }
  }
  
  # scale and legend
  full <- x$methods$full
  model <- x$methods$model
  EMint <- x$variables$EMint
  if (model == "iorw") {
    if(x$multimp$multimp) yreg_mid <- x$reg.output[[1]]$yregTot
    if(!x$multimp$multimp) yreg_mid <- x$reg.output$yregTot
  } else {
    if(x$multimp$multimp) yreg_mid <- x$reg.output[[1]]$yreg
    if(!x$multimp$multimp) yreg_mid <- x$reg.output$yreg
  }
  if (inherits(yreg_mid, "rcreg") | inherits(yreg_mid, "simexreg")) yreg_mid <- yreg_mid$NAIVEreg
  is_lm <- inherits(yreg_mid, "lm")
  is_glm <- inherits(yreg_mid, "glm")
  is_svyglm <- inherits(yreg_mid, "svyglm")
  is_gam <- inherits(yreg_mid, "gam")
  if (is_lm | is_glm) family_yreg <- family(yreg_mid)
  is_multinom <- inherits(yreg_mid, "multinom")
  is_svymultinom <- inherits(yreg_mid, "svymultinom")
  is_polr <- inherits(yreg_mid, "polr")
  is_survreg <- inherits(yreg_mid, "survreg")
  is_coxph <- inherits(yreg_mid, "coxph")
  if ((is_lm | is_glm) && (family_yreg$family %in% c("gaussian", "inverse.gaussian", "quasi", "Gamma"))) {
    scale <- "mean difference scale"
    if (model == "iorw") {
      if (full) legend <- "(te: total effect; pnde: pure natural direct effect; tnie: total natural indirect effect; pm: proportion mediated)"   
      if (!full) legend <- "(te: total effect; pnde: pure natural direct effect; tnie: total natural indirect effect)"
    } else if (length(x$variables$postc) != 0) {
      if (full) {
        if (EMint) legend <- "(cde: controlled direct effect; rpnde: randomized analogue of pure natural direct effect; rtnde: randomized analogue of total natural direct effect; rpnie: randomized analogue of pure natural indirect effect; rtnie: randomized analogue of total natural indirect effect; te: total effect; rintref: randomized analogue of reference interaction; rintmed: randomized analogue of mediated interaction; cde(prop): proportion cde; rintref(prop): proportion rintref; rintmed(prop): proportion rintmed; rpnie(prop): proportion rpnie; rpm: randomized analogue of overall proportion mediated; rint: randomized analogue of overall proportion attributable to interaction; rpe: randomized analogue of overall proportion eliminated)"
        if (!EMint) legend <- "(cde: controlled direct effect; rpnde: randomized analogue of pure natural direct effect; rtnde: randomized analogue of total natural direct effect; rpnie: randomized analogue of pure natural indirect effect; rtnie: randomized analogue of total natural indirect effect; te: total effect; rpm: randomized analogue of overall proportion mediated)"
      } else legend <- "(cde: controlled direct effect; rpnde: randomized analogue of pure natural direct effect; rtnde: randomized analogue of total natural direct effect; rpnie: randomized analogue of pure natural indirect effect; rtnie: randomized analogue of total natural indirect effect; te: total effect)"
    } else {
      if (full) {
        if (EMint) legend <- "(cde: controlled direct effect; pnde: pure natural direct effect; tnde: total natural direct effect; pnie: pure natural indirect effect; tnie: total natural indirect effect; te: total effect; intref: reference interaction; intmed: mediated interaction; cde(prop): proportion cde; intref(prop): proportion intref; intmed(prop): proportion intmed; pnie(prop): proportion pnie; pm: overall proportion mediated; int: overall proportion attributable to interaction; pe: overall proportion eliminated)"
        if (!EMint) legend <- "(cde: controlled direct effect; pnde: pure natural direct effect; tnde: total natural direct effect; pnie: pure natural indirect effect; tnie: total natural indirect effect; te: total effect; pm: overall proportion mediated)"
      } else legend <- "(cde: controlled direct effect; pnde: pure natural direct effect; tnde: total natural direct effect; pnie: pure natural indirect effect; tnie: total natural indirect effect; te: total effect)"
    }
  } else if ((is_lm | is_glm) && (family_yreg$family %in% c("poisson", "quasipoisson", "ziplss") |
                                  startsWith(family_yreg$family, "Negative Binomial") |
                                  startsWith(family_yreg$family, "Zero inflated Poisson"))) {
    scale <- "rate ratio scale"
    if (model == "iorw") {
      if (full) legend <- "(Rte: total effect rate ratio; Rpnde: pure natural direct effect rate ratio; Rtnie: total natural indirect effect rate ratio; pm: proportion mediated)"   
      if (!full) legend <- "(Rte: total effect rate ratio; Rpnde: pure natural direct effect rate ratio; Rtnie: total natural indirect effect rate ratio)" 
    } else if (length(x$variables$postc) != 0) {
      if (full) {
        if (EMint) legend <- "(Rcde: controlled direct effect rate ratio; rRpnde: randomized analogue of pure natural direct effect rate ratio; rRtnde: randomized analogue of total natural direct effect rate ratio; rRpnie: randomized analogue of pure natural indirect effect rate ratio; rRtnie: randomized analogue of total natural indirect effect rate ratio; Rte: total effect rate ratio; ERcde: excess relative rate due to controlled direct effect; rERintref: randomized analogue of excess relative rate due to reference interaction; rERintmed: randomized analogue of excess relative rate due to mediated interaction; rERpnie: randomized analogue of excess relative rate due to pure natural indirect effect; ERcde(prop): proportion ERcde; rERintref(prop): proportion rERintref; rERintmed(prop): proportion rERintmed; rERpnie(prop): proportion rERpnie; rpm: randomized analogue of overall proportion mediated; rint: randomized analogue of overall proportion attributable to interaction; rpe: randomized analogue of overall proportion eliminated)"
        if (!EMint) legend <- "(Rcde: controlled direct effect rate ratio; rRpnde: randomized analogue of pure natural direct effect rate ratio; rRtnde: randomized analogue of total natural direct effect rate ratio; rRpnie: randomized analogue of pure natural indirect effect rate ratio; rRtnie: randomized analogue of total natural indirect effect rate ratio; Rte: total effect rate ratio; rpm: randomized analogue of overall proportion mediated)"
      } else legend <- "(Rcde: controlled direct effect rate ratio; rRpnde: randomized analogue of pure natural direct effect rate ratio; rRtnde: randomized analogue of total natural direct effect rate ratio; rRpnie: randomized analogue of pure natural indirect effect rate ratio; rRtnie: randomized analogue of total natural indirect effect rate ratio; Rte: total effect rate ratio)"
    } else {
      if (full) {
        if (EMint) legend <- "(Rcde: controlled direct effect rate ratio; Rpnde: pure natural direct effect rate ratio; Rtnde: total natural direct effect rate ratio; Rpnie: pure natural indirect effect rate ratio; Rtnie: total natural indirect effect rate ratio; Rte: total effect rate ratio; ERcde: excess relative rate due to controlled direct effect; ERintref: excess relative rate due to reference interaction; ERintmed: excess relative rate due to mediated interaction; ERpnie: excess relative rate due to pure natural indirect effect; ERcde(prop): proportion ERcde; ERintref(prop): proportion ERintref; ERintmed(prop): proportion ERintmed; ERpnie(prop): proportion ERpnie; pm: overall proportion mediated; int: overall proportion attributable to interaction; pe: overall proportion eliminated)"
        if (!EMint) legend <- "(Rcde: controlled direct effect rate ratio; Rpnde: pure natural direct effect rate ratio; Rtnde: total natural direct effect rate ratio; Rpnie: pure natural indirect effect rate ratio; Rtnie: total natural indirect effect rate ratio; Rte: total effect rate ratio; pm: overall proportion mediated)"
      } else legend <- "(Rcde: controlled direct effect rate ratio; Rpnde: pure natural direct effect rate ratio; Rtnde: total natural direct effect rate ratio; Rpnie: pure natural indirect effect rate ratio; Rtnie: total natural indirect effect rate ratio; Rte: total effect rate ratio)"
    }
  } else if (((is_lm | is_glm) && (family_yreg$family %in% c("binomial", "quasibinomial", "multinom") |
                                   startsWith(family_yreg$family, "Ordered Categorical"))) |
             is_multinom | is_polr) {
    
    if (is_glm && family_yreg$family %in% c("binomial", "quasibinomial") && yreg_mid$family$link == "logit") {
      scale <- "odds ratio scale"
      if (model == "iorw") {
        if (full) legend <- "(Rte: total effect odds ratio; Rpnde: pure natural direct effect odds ratio; Rtnie: total natural indirect effect odds ratio; pm: proportion mediated)"   
        if (!full) legend <- "(Rte: total effect odds ratio; Rpnde: pure natural direct effect odds ratio; Rtnie: total natural indirect effect odds ratio)" 
      } else if (length(x$variables$postc) != 0) {
        if (full) {
          if (EMint) legend <- "(Rcde: controlled direct effect odds ratio; rRpnde: randomized analogue of pure natural direct effect odds ratio; rRtnde: randomized analogue of total natural direct effect odds ratio; rRpnie: randomized analogue of pure natural indirect effect odds ratio; rRtnie: randomized analogue of total natural indirect effect odds ratio; Rte: total effect odds ratio; ERcde: excess relative risk due to controlled direct effect; rERintref: randomized analogue of excess relative risk due to reference interaction; rERintmed: randomized analogue of excess relative risk due to mediated interaction; rERpnie: randomized analogue of excess relative risk due to pure natural indirect effect; ERcde(prop): proportion ERcde; rERintref(prop): proportion rERintref; rERintmed(prop): proportion rERintmed; rERpnie(prop): proportion rERpnie; rpm: randomized analogue of overall proportion mediated; rint: randomized analogue of overall proportion attributable to interaction; rpe: randomized analogue of overall proportion eliminated)"
          if (!EMint) legend <- "(Rcde: controlled direct effect odds ratio; rRpnde: randomized analogue of pure natural direct effect odds ratio; rRtnde: randomized analogue of total natural direct effect odds ratio; rRpnie: randomized analogue of pure natural indirect effect odds ratio; rRtnie: randomized analogue of total natural indirect effect odds ratio; Rte: total effect odds ratio; rpm: randomized analogue of overall proportion mediated)"
        } else legend <- "(Rcde: controlled direct effect odds ratio; rRpnde: randomized analogue of pure natural direct effect odds ratio; rRtnde: randomized analogue of total natural direct effect odds ratio; rRpnie: randomized analogue of pure natural indirect effect odds ratio; rRtnie: randomized analogue of total natural indirect effect odds ratio; Rte: total effect odds ratio)"
      } else {
        if (full) {
          if (EMint) legend <- "(Rcde: controlled direct effect odds ratio; Rpnde: pure natural direct effect odds ratio; Rtnde: total natural direct effect odds ratio; Rpnie: pure natural indirect effect odds ratio; Rtnie: total natural indirect effect odds ratio; Rte: total effect odds ratio; ERcde: excess relative risk due to controlled direct effect; ERintref: excess relative risk due to reference interaction; ERintmed: excess relative risk due to mediated interaction; ERpnie: excess relative risk due to pure natural indirect effect; ERcde(prop): proportion ERcde; ERintref(prop): proportion ERintref; ERintmed(prop): proportion ERintmed; ERpnie(prop): proportion ERpnie; pm: overall proportion mediated; int: overall proportion attributable to interaction; pe: overall proportion eliminated)"
          if (!EMint) legend <- "(Rcde: controlled direct effect odds ratio; Rpnde: pure natural direct effect odds ratio; Rtnde: total natural direct effect odds ratio; Rpnie: pure natural indirect effect odds ratio; Rtnie: total natural indirect effect odds ratio; Rte: total effect odds ratio; pm: overall proportion mediated)"
        } else legend <- "(Rcde: controlled direct effect odds ratio; Rpnde: pure natural direct effect odds ratio; Rtnde: total natural direct effect odds ratio; Rpnie: pure natural indirect effect odds ratio; Rtnie: total natural indirect effect odds ratio; Rte: total effect odds ratio)"
      }
    } else {
      scale <- "risk ratio scale"
      if (model == "iorw") {
        if (full) legend <- "(Rte: total effect risk ratio; Rpnde: pure natural direct effect risk ratio; Rtnie: total natural indirect effect risk ratio; pm: proportion mediated)"   
        if (!full) legend <- "(Rte: total effect risk ratio; Rpnde: pure natural direct effect risk ratio; Rtnie: total natural indirect effect risk ratio)" 
      } else if (length(x$variables$postc) != 0) {
        if (full) {
          if (EMint) legend <- "(Rcde: controlled direct effect risk ratio; rRpnde: randomized analogue of pure natural direct effect risk ratio; rRtnde: randomized analogue of total natural direct effect risk ratio; rRpnie: randomized analogue of pure natural indirect effect risk ratio; rRtnie: randomized analogue of total natural indirect effect risk ratio; Rte: total effect risk ratio; ERcde: excess relative risk due to controlled direct effect; rERintref: randomized analogue of excess relative risk due to reference interaction; rERintmed: randomized analogue of excess relative risk due to mediated interaction; rERpnie: randomized analogue of excess relative risk due to pure natural indirect effect; ERcde(prop): proportion ERcde; rERintref(prop): proportion rERintref; rERintmed(prop): proportion rERintmed; rERpnie(prop): proportion rERpnie; rpm: randomized analogue of overall proportion mediated; rint: randomized analogue of overall proportion attributable to interaction; rpe: randomized analogue of overall proportion eliminated)"
          if (!EMint) legend <- "(Rcde: controlled direct effect risk ratio; rRpnde: randomized analogue of pure natural direct effect risk ratio; rRtnde: randomized analogue of total natural direct effect risk ratio; rRpnie: randomized analogue of pure natural indirect effect risk ratio; rRtnie: randomized analogue of total natural indirect effect risk ratio; Rte: total effect risk ratio; rpm: randomized analogue of overall proportion mediated)"
        } else legend <- "(Rcde: controlled direct effect risk ratio; rRpnde: randomized analogue of pure natural direct effect risk ratio; rRtnde: randomized analogue of total natural direct effect risk ratio; rRpnie: randomized analogue of pure natural indirect effect risk ratio; rRtnie: randomized analogue of total natural indirect effect risk ratio; Rte: total effect risk ratio)"
      } else {
        if (full) {
          if (EMint) legend <- "(Rcde: controlled direct effect risk ratio; Rpnde: pure natural direct effect risk ratio; Rtnde: total natural direct effect risk ratio; Rpnie: pure natural indirect effect risk ratio; Rtnie: total natural indirect effect risk ratio; Rte: total effect risk ratio; ERcde: excess relative risk due to controlled direct effect; ERintref: excess relative risk due to reference interaction; ERintmed: excess relative risk due to mediated interaction; ERpnie: excess relative risk due to pure natural indirect effect; ERcde(prop): proportion ERcde; ERintref(prop): proportion ERintref; ERintmed(prop): proportion ERintmed; ERpnie(prop): proportion ERpnie; pm: overall proportion mediated; int: overall proportion attributable to interaction; pe: overall proportion eliminated)"
          if (!EMint) legend <- "(Rcde: controlled direct effect risk ratio; Rpnde: pure natural direct effect risk ratio; Rtnde: total natural direct effect risk ratio; Rpnie: pure natural indirect effect risk ratio; Rtnie: total natural indirect effect risk ratio; Rte: total effect risk ratio; pm: overall proportion mediated)"
        } else legend <- "(Rcde: controlled direct effect risk ratio; Rpnde: pure natural direct effect risk ratio; Rtnde: total natural direct effect risk ratio; Rpnie: pure natural indirect effect risk ratio; Rtnie: total natural indirect effect risk ratio; Rte: total effect risk ratio)"
      }
    }
    
  } else if (is_coxph) {
    scale <- "hazard ratio scale"
    if (model == "iorw") {
      if (full) legend <- "(Rte: total effect hazard ratio; Rpnde: pure natural direct effect hazard ratio; Rtnie: total natural indirect effect hazard ratio; pm: proportion mediated)"   
      if (!full) legend <- "(Rte: total effect hazard ratio; Rpnde: pure natural direct effect hazard ratio; Rtnie: total natural indirect effect hazard ratio)" 
    }  else if (length(x$variables$postc) != 0) {
      if (full) {
        if (EMint) legend <- "(Rcde: controlled direct effect hazard ratio; rRpnde: randomized analogue of pure natural direct effect hazard ratio; rRtnde: randomized analogue of total natural direct effect hazard ratio; rRpnie: randomized analogue of pure natural indirect effect hazard ratio; rRtnie: randomized analogue of total natural indirect effect hazard ratio; Rte: total effect hazard ratio; ERcde: excess relative hazard due to controlled direct effect; rERintref: randomized analogue of excess relative hazard due to reference interaction; rERintmed: randomized analogue of excess relative hazard due to mediated interaction; rERpnie: randomized analogue of excess relative hazard due to pure natural indirect effect; ERcde(prop): proportion ERcde; rERintref(prop): proportion rERintref; rERintmed(prop): proportion rERintmed; rERpnie(prop): proportion rERpnie; rpm: randomized analogue of overall proportion mediated; rint: randomized analogue of overall proportion attributable to interaction; rpe: randomized analogue of overall proportion eliminated)"
        if (!EMint) legend <- "(Rcde: controlled direct effect hazard ratio; rRpnde: randomized analogue of pure natural direct effect hazard ratio; rRtnde: randomized analogue of total natural direct effect hazard ratio; rRpnie: randomized analogue of pure natural indirect effect hazard ratio; rRtnie: randomized analogue of total natural indirect effect hazard ratio; Rte: total effect hazard ratio; rpm: randomized analogue of overall proportion mediated)"
      } else legend <- "(Rcde: controlled direct effect hazard ratio; rRpnde: randomized analogue of pure natural direct effect hazard ratio; rRtnde: randomized analogue of total natural direct effect hazard ratio; rRpnie: randomized analogue of pure natural indirect effect hazard ratio; rRtnie: randomized analogue of total natural indirect effect hazard ratio; Rte: total effect hazard ratio)"
    } else {
      if (full) {
        if (EMint) legend <- "(Rcde: controlled direct effect hazard ratio; Rpnde: pure natural direct effect hazard ratio; Rtnde: total natural direct effect hazard ratio; Rpnie: pure natural indirect effect hazard ratio; Rtnie: total natural indirect effect hazard ratio; Rte: total effect hazard ratio; ERcde: excess relative hazard due to controlled direct effect; ERintref: excess relative hazard due to reference interaction; ERintmed: excess relative hazard due to mediated interaction; ERpnie: excess relative hazard due to pure natural indirect effect; ERcde(prop): proportion ERcde; ERintref(prop): proportion ERintref; ERintmed(prop): proportion ERintmed; ERpnie(prop): proportion ERpnie; pm: overall proportion mediated; int: overall proportion attributable to interaction; pe: overall proportion eliminated)"
        if (!EMint) legend <- "(Rcde: controlled direct effect hazard ratio; Rpnde: pure natural direct effect hazard ratio; Rtnde: total natural direct effect hazard ratio; Rpnie: pure natural indirect effect hazard ratio; Rtnie: total natural indirect effect hazard ratio; Rte: total effect hazard ratio; pm: overall proportion mediated)"
      } else legend <- "(Rcde: controlled direct effect hazard ratio; Rpnde: pure natural direct effect hazard ratio; Rtnde: total natural direct effect hazard ratio; Rpnie: pure natural indirect effect hazard ratio; Rtnie: total natural indirect effect hazard ratio; Rte: total effect hazard ratio)"
    }
  } else if (is_survreg) {
    scale <- "mean survival scale"
    if (model == "iorw") {
      if (full) legend <- "(Rte: total effect mean survival ratio; Rpnde: pure natural direct effect mean survival ratio; Rtnie: total natural indirect effect mean survival ratio; pm: proportion mediated)"   
      if (!full) legend <- "(Rte: total effect mean survival ratio; Rpnde: pure natural direct effect mean survival ratio; Rtnie: total natural indirect effect mean survival ratio)" 
    } else if (length(x$variables$postc) != 0) {
      if (full) {
        if (EMint) legend <- "(Rcde: controlled direct effect mean survival ratio; rRpnde: randomized analogue of pure natural direct effect mean survival ratio; rRtnde: randomized analogue of total natural direct effect mean survival ratio; rRpnie: randomized analogue of pure natural indirect effect mean survival ratio; rRtnie: randomized analogue of total natural indirect effect mean survival ratio; Rte: total effect mean survival ratio; ERcde: excess mean survival ratio due to controlled direct effect; rERintref: randomized analogue of excess mean survival ratio due to reference interaction; rERintmed: randomized analogue of excess mean survival ratio due to mediated interaction; rERpnie: randomized analogue of excess mean survival ratio due to pure natural indirect effect; ERcde(prop): proportion ERcde; rERintref(prop): proportion rERintref; rERintmed(prop): proportion rERintmed; rERpnie(prop): proportion rERpnie; rpm: randomized analogue of overall proportion mediated; rint: randomized analogue of overall proportion attributable to interaction; rpe: randomized analogue of overall proportion eliminated)"
        if (!EMint) legend <- "(Rcde: controlled direct effect mean survival ratio; rRpnde: randomized analogue of pure natural direct effect mean survival ratio; rRtnde: randomized analogue of total natural direct effect mean survival ratio; rRpnie: randomized analogue of pure natural indirect effect mean survival ratio; rRtnie: randomized analogue of total natural indirect effect mean survival ratio; Rte: total effect mean survival ratio; rpm: randomized analogue of overall proportion mediated)"
      } else legend <- "(Rcde: controlled direct effect mean survival ratio; rRpnde: randomized analogue of pure natural direct effect mean survival ratio; rRtnde: randomized analogue of total natural direct effect mean survival ratio; rRpnie: randomized analogue of pure natural indirect effect mean survival ratio; rRtnie: randomized analogue of total natural indirect effect mean survival ratio; Rte: total effect mean survival ratio)"
    } else {
      if (full) {
        if (EMint) legend <- "(Rcde: controlled direct effect mean survival ratio; Rpnde: pure natural direct effect mean survival ratio; Rtnde: total natural direct effect mean survival ratio; Rpnie: pure natural indirect effect mean survival ratio; Rtnie: total natural indirect effect mean survival ratio; Rte: total effect mean survival ratio; ERcde: excess mean survival ratio due to controlled direct effect; ERintref: excess mean survival ratio due to reference interaction; ERintmed: excess mean survival ratio due to mediated interaction; ERpnie: excess mean survival ratio due to pure natural indirect effect; ERcde(prop): proportion ERcde; ERintref(prop): proportion ERintref; ERintmed(prop): proportion ERintmed; ERpnie(prop): proportion ERpnie; pm: overall proportion mediated; int: overall proportion attributable to interaction; pe: overall proportion eliminated)"
        if (!EMint) legend <- "(Rcde: controlled direct effect mean survival ratio; Rpnde: pure natural direct effect mean survival ratio; Rtnde: total natural direct effect mean survival ratio; Rpnie: pure natural indirect effect mean survival ratio; Rtnie: total natural indirect effect mean survival ratio; Rte: total effect mean survival ratio; pm: overall proportion mediated)"
      } else legend <- "(Rcde: controlled direct effect mean survival ratio; Rpnde: pure natural direct effect mean survival ratio; Rtnde: total natural direct effect mean survival ratio; Rpnie: pure natural indirect effect mean survival ratio; Rtnie: total natural indirect effect mean survival ratio; Rte: total effect mean survival ratio)"
    }
  }
  
  # print summary of causal mediation analysis results
  if (x$methods$model == "rb") model_str <- "regression-based approach"
  if (x$methods$model == "wb") model_str <- "weighting-based approach"
  if (x$methods$model == "ne") model_str <- "natural effect model"
  if (x$methods$model == "iorw") model_str <- "inverse odds ratio weighting approach"
  if (x$methods$model == "msm") model_str <- "marginal structural model"
  if (x$methods$model == "gformula") model_str <- "g-formula approach"
  if (x$multimp$multimp) model_str <- paste(model_str, "with multiple imputation")
  if (x$methods$estimation == "paramfunc") est_str <- "Closed-form parameter function estimation"
  if (x$methods$estimation == "imputation") est_str <- "Direct counterfactual imputation estimation"
  if (x$methods$inference == "delta") inf_str <- "delta method standard errors, confidence intervals and p-values"
  if (x$methods$inference == "bootstrap") {
    if (x$methods$boot.ci.type == "per") inf_str <- "bootstrap standard errors, percentile confidence intervals and p-values"
    if (x$methods$boot.ci.type == "bca") inf_str <- "bootstrap standard errors, bias-corrected and accelerated confidence intervals and p-values"
  }
  if (x$methods$casecontrol) cat(paste("# Effect decomposition on the", scale, "for a case control study via the "))
  if (!(x$methods$casecontrol)) cat(paste("# Effect decomposition on the", scale, "via the "))
  cat(model_str)
  cat("\n \n")
  cat(est_str)
  cat(paste(" with \n", inf_str, "\n \n"))
  printCoefmat(x$summarydf, digits = digits, has.Pvalue = TRUE)
  cat("\n")
  cat(legend)
  cat("\n\nRelevant variable values: \n")
  print(x$ref)
}


#' Plotting Point Estimates and Confidence Intervals of Causal Effects
#' 
#' \code{ggcmest} is used to plot results of \code{cmest} nicely with plotting functions
#' in the \link{ggplot2} package. Additional layers can be added to this plot using other 
#' plotting functions in the \link{ggplot2} package.
#' 
#' @param x an object of class \code{cmest}
#' @param errorbar.width width of errorbars for confidence intervals. Default is \code{0.3}.
#' @param errorbar.size size of errorbars for confidence intervals. Default is \code{0.3}.
#' @param errorbar.colour colour of errorbars for confidence intervals. Default is \code{black}.
#' @param point.size size of points for point estimates. Default is \code{1}.
#' @param point.colour colour of points for point estimates. Default is \code{blue}.
#' @param refline a logical value. If \code{true}, include a reference line at 
#' \code{y = 0} when effects are on the difference scale and include a reference line at 
#' \code{y = 1} when effects are on the ratio scale. Default is \code{TRUE}.
#' @param refline.colour colour of the reference line. Default is \code{red}.
#' @param refline.size size of the reference line. Default is \code{0.3}.
#' 
#' @seealso \code{\link{cmest}}, \code{\link{ggplot2}}.
#' 
#' @examples
#' 
#' library(CMAverse)
#' library(ggplot2)
#' 
#' x <- cmest(data = cma2020, model = "rb", outcome = "contY", 
#' exposure = "A", mediator = "M2", basec = c("C1", "C2"), 
#' EMint = TRUE, mreg = list("multinomial"), yreg = "linear", 
#' astar = 0, a = 1, mval = list("M2_0"), estimation = "paramfunc", 
#' inference = "delta")
#' 
#' ggcmest(x) +
#' theme(axis.text.x = element_text(angle = 45))
#' 
#' ggcmest(x) +
#' coord_flip(xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")
#' 
#' @export
#' 
ggcmest <- function(x, errorbar.width = 0.3, errorbar.size = 0.3, errorbar.colour = "black",
                    point.size = 1, point.colour = "blue", 
                    refline = TRUE, refline.colour = "red", refline.size = 0.3) {
  # create a data frame for results of cmest
  effect_df <- data.frame(Effect = factor(names(x$effect.pe), levels = names(x$effect.pe)),
                          PE = x$effect.pe, CIlower = x$effect.ci.low,
                          CIupper = x$effect.ci.high)
  # reference line
  if (refline) {
    if (!x$multimp$multimp) {
      if ((inherits(x$reg.output$yreg, "lm") | inherits(x$reg.output$yreg, "glm")) &&
          (family(x$reg.output$yreg)$family %in% c("gaussian","Gamma","inverse.gaussian","quasi"))) {
        ref <- 0
      } else {
        if (x$methods$full) ref <- c(0, 1)
        if (!x$methods$full) ref <- 1
      }
    } else {
      if ((inherits(x$reg.output[[1]]$yreg, "lm") | inherits(x$reg.output[[1]]$yreg, "glm")) &&
          (family(x$reg.output[[1]]$yreg)$family %in% c("gaussian","Gamma","inverse.gaussian","quasi"))) {
        ref <- 0
      } else {
        if (x$methods$full) ref <- c(0, 1)
        if (!x$methods$full) ref <- 1
      }
    }
  } else ref <- NULL
  # plot
  ggplot() +
    geom_errorbar(aes(x = Effect, ymin = CIlower, ymax = CIupper),
                  width = errorbar.width, size = errorbar.size, colour = errorbar.colour,
                  data = effect_df) +
    geom_point(aes(x = Effect, y = PE),
               size = point.size, colour = point.colour, data = effect_df) +
    ylab("Point Estimate and 95% CI") +
    geom_hline(yintercept = ref, color = refline.colour, size = refline.size)
  
}

regrun <- function() {
  # Y: outcome; M: mediator; A: exposure; C: basec; L: postc
  # the variable used to calculate weights is required to be categorical
  
  # Exposure Regression For Weighting
  # for wb and msm, the exposure regression is required for calculating weights if basec is not empty, w_{a,i}=P(A=A_i)/P(A=A_i|C_i)
  # for iorw, the exposure regression is required for calculating weights, w_{a,i}=P(A=0|M_i,C_i)/P(A=A_i|M_i,C_i)
  if ((model %in% c("wb", "msm") && length(basec) > 0) | model == "iorw") {
    if (is.null(ereg)) stop("ereg is required when model is 'wb' or 'msm' with length(basec) > 0 and when model is 'iorw'")
    if (is.character(ereg)) {
      # fit glm with family = poisson() rather than family = binomial("log") for "loglinear"
      if (ereg == "loglinear" && length(unique(na.omit(data[, exposure]))) != 2) stop("When ereg is 'loglinear', exposure should be binary")
      if (!ereg %in% c("logistic", "loglinear", "multinomial", "ordinal")) stop("Select character ereg from 'logistic', 'loglinear', 'multinomial', 'ordinal'")
      exposure_formula <- switch((model == "iorw") + 1, "1" = paste0(exposure, "~", paste0(basec, collapse = "+")),
                                 "2" = paste0(exposure, "~", paste0(c(mediator, basec), collapse = "+")))
      switch(ereg,
             logistic = ereg <- eval(bquote(glm(.(as.formula(exposure_formula)), family = binomial(), data = .(data)))),
             loglinear = ereg <- eval(bquote(glm(.(as.formula(exposure_formula)), family = poisson(), data = .(data)))),
             multinomial = ereg <- eval(bquote(nnet::multinom(.(as.formula(exposure_formula)), data = .(data), trace = FALSE))),
             ordinal = ereg <- eval(bquote(MASS::polr(.(as.formula(exposure_formula)), data = .(data)))))
    } else if (!(identical(class(ereg), "lm") | identical(class(ereg), c("glm", "lm")) |
                 identical(class(ereg), c("negbin", "glm", "lm")) | identical(class(ereg), c("multinom", "nnet")) |
                 identical(class(ereg), c("gam", "glm", "lm")) | identical(class(ereg), "polr"))) {
      stop("Fit ereg by lm, glm, glm.nb, gam, multinom, polr")
    } 
  } else {
    if (!is.null(ereg)) warning("ereg is ignored when model is 'wb' or 'msm' with length(basec) = 0 or model is 'rb', 'ne' or 'gformula'")
    ereg <- NULL
  }
  
  # Mediator Regression For Weighting
  # for msm, a mediator regression for weighting is required for each mediator
  # w_{m_p,i}=P(M_p=M_{p,i}|A=A_i,M_1=M_{1,i},...,M_{p-1}=M_{p-1,i})/P(M_p=M_{p,i}|A=A_i,C=C_i,L=L_i,M_1=M_{1,i},...,M_{p-1}=M_{p-1,i})
  if (model == "msm") {
    if (!is.list(wmnomreg)) stop("wmnomreg should be a list")
    if (length(wmnomreg) != length(mediator)) stop("length(wmnomreg) != length(mediator)")
    for (p in 1:length(wmnomreg)) {
      if (is.null(wmnomreg[[p]])) stop(paste0("Unspecified wmnomreg[[", p, "]]"))
      if (is.character(wmnomreg[[p]])) {
        if (wmnomreg[[p]] == "loglinear" && length(unique(na.omit(data[, mediator[p]]))) != 2) stop(paste0("When wmnomreg[[", p, "]] is 'loglinear', mediator[[", p, "]] should be binary"))
        if (!wmnomreg[[p]] %in% c("logistic", "loglinear", "multinomial", "ordinal")) stop(paste0("Select character wmnomreg[[", p, "]] from 'logistic', 'loglinear', 'multinomial', 'ordinal'"))
        wmnomreg_formula <- paste0(mediator[p], "~", paste(c(exposure, mediator[0:(p-1)]), collapse = "+"))
        # regression for the nominator of w_{m_p,i}
        switch(wmnomreg[[p]],
               logistic = wmnomreg[[p]] <- eval(bquote(glm(.(as.formula(wmnomreg_formula)), family = binomial(), data = .(data)))),
               loglinear = wmnomreg[[p]] <- eval(bquote(glm(.(as.formula(wmnomreg_formula)), family = poisson(), data = .(data)))),
               multinomial = wmnomreg[[p]] <- eval(bquote(nnet::multinom(.(as.formula(wmnomreg_formula)), data = .(data), trace = FALSE))),
               ordinal = wmnomreg[[p]] <- eval(bquote(MASS::polr(.(as.formula(wmnomreg_formula)), data = .(data)))))
      } else if (!(identical(class(wmnomreg[[p]]), "lm") | identical(class(wmnomreg[[p]]), c("glm", "lm")) |
                   identical(class(wmnomreg[[p]]), c("negbin", "glm", "lm")) | identical(class(wmnomreg[[p]]), c("multinom", "nnet")) |
                   identical(class(wmnomreg[[p]]), c("gam", "glm", "lm")) | identical(class(wmnomreg[[p]]), "polr"))) {
        stop(paste0("Fit wmnomreg[[", p, "]] by lm, glm, glm.nb, gam, multinom, polr"))
      } 
    }
    if (!is.list(wmdenomreg)) stop("wmdenomreg should be a list")
    if (length(wmdenomreg) != length(mediator)) stop("length(wmdenomreg) != length(mediator)")
    for (p in 1:length(wmdenomreg)) {
      if (is.null(wmdenomreg[[p]])) stop(paste0("Unspecified wmdenomreg[[", p, "]]"))
      if (is.character(wmdenomreg[[p]])) {
        if (wmdenomreg[[p]] == "loglinear" && length(unique(na.omit(data[, mediator[p]]))) != 2) stop(paste0("When wmdenomreg[[", p, "]] is 'loglinear', mediator[[", p, "]] should be binary"))
        if (!wmdenomreg[[p]] %in% c("logistic", "loglinear", "multinomial", "ordinal")) stop(paste0("Select character wmdenomreg[[", p, "]] from 'logistic', 'loglinear', 'multinomial', 'ordinal'"))
        wmdenomreg_formula <- paste0(mediator[p], "~", paste(c(exposure, mediator[0:(p-1)], basec, postc), collapse = "+"))
        # regression for the denominator of w_{m_p,i}
        switch(wmdenomreg[[p]],
               logistic = wmdenomreg[[p]] <- eval(bquote(glm(.(as.formula(wmdenomreg_formula)), family = binomial(), data = .(data)))),
               loglinear = wmdenomreg[[p]] <- eval(bquote(glm(.(as.formula(wmdenomreg_formula)), family = poisson(), data = .(data)))),
               multinomial = wmdenomreg[[p]] <- eval(bquote(nnet::multinom(.(as.formula(wmdenomreg_formula)), data = .(data), trace = FALSE))),
               ordinal = wmdenomreg[[p]] <- eval(bquote(MASS::polr(.(as.formula(wmdenomreg_formula)), data = .(data)))))
      } else if (!(identical(class(wmdenomreg[[p]]), "lm") | identical(class(wmdenomreg[[p]]), c("glm", "lm")) |
                   identical(class(wmdenomreg[[p]]), c("negbin", "glm", "lm")) | identical(class(wmdenomreg[[p]]), c("multinom", "nnet")) |
                   identical(class(wmdenomreg[[p]]), c("gam", "glm", "lm")) | identical(class(wmdenomreg[[p]]), "polr"))) {
        stop(paste0("Fit wmdenomreg[[", p, "]] by lm, glm, glm.nb, gam, multinom, polr"))
      } 
    }
  } else {
    if (!is.null(wmnomreg)) warning("wmnomreg is ignored when model is not 'msm'")
    if (!is.null(wmdenomreg)) warning("wmdenomreg is ignored when model is not 'msm'")
    wmnomreg <- wmdenomreg <- NULL
  }
  
  # Mediator Regression
  # for rb, msm and gformula, a mediator regression is required for each mediator
  if (model %in% c("rb", "msm", "gformula")) {
    if (!is.list(mreg)) stop("mreg should be a list")
    if (length(mreg) != length(mediator)) stop("length(mreg) != length(mediator)")
    for (p in 1:length(mreg)) {
      if (is.null(mreg[[p]])) stop(paste0("Unspecified mreg[[", p, "]]"))
      if (is.character(mreg[[p]])) {
        if (mreg[[p]] == "loglinear" && length(unique(na.omit(data[, mediator[p]]))) != 2) stop(paste0("When mreg[[", p, "]] is 'loglinear', mediator[[", p, "]] should be binary"))
        if (!mreg[[p]] %in% c("linear", "logistic", "loglinear", "poisson", "quasipoisson",
                              "negbin", "multinomial", "ordinal")) stop(
                                paste0("Select character mreg[[", p, "]] from 'linear', 'logistic',
                                       'loglinear', 'poisson', 'quasipoisson', 'negbin', 'multinomial', 'ordinal'"))
        # for rb, regress each mediator on A, C
        # for msm, regress each mediator on A
        # for gformula, regress each mediator on A, C, L
        switch(model,
               rb = mediator_formula <- paste0(mediator[p], "~", paste(c(exposure, basec), collapse = "+")),
               msm = mediator_formula <- paste0(mediator[p], "~", exposure),
               gformula = mediator_formula <- paste0(mediator[p], "~", paste(c(exposure, basec, postc), collapse = "+")))
        switch(mreg[[p]],
               linear = mreg[[p]] <- eval(bquote(glm(.(as.formula(mediator_formula)), family = gaussian(), data = .(data)))),
               logistic = mreg[[p]] <- eval(bquote(glm(.(as.formula(mediator_formula)), family = binomial(), data = .(data)))),
               loglinear = mreg[[p]] <- eval(bquote(glm(.(as.formula(mediator_formula)), family = poisson(), data = .(data)))),
               poisson = mreg[[p]]  <- eval(bquote(glm(.(as.formula(mediator_formula)), family = poisson(), data = .(data)))),
               quasipoisson = mreg[[p]] <- eval(bquote(glm(.(as.formula(mediator_formula)), family = quasipoisson(), data = .(data)))),
               negbin = mreg[[p]] <- eval(bquote(MASS::glm.nb(.(as.formula(mediator_formula)), data = .(data)))),
               multinomial = mreg[[p]] <- eval(bquote(nnet::multinom(.(as.formula(mediator_formula)), data = .(data), trace = FALSE))),
               ordinal = mreg[[p]] <- eval(bquote(MASS::polr(.(as.formula(mediator_formula)), data = .(data)))))
      } else if (!(identical(class(mreg[[p]]), "lm") | identical(class(mreg[[p]]), c("glm", "lm")) |
                   identical(class(mreg[[p]]), c("negbin", "glm", "lm")) | identical(class(mreg[[p]]), c("multinom", "nnet")) |
                   identical(class(mreg[[p]]), c("gam", "glm", "lm")) | identical(class(mreg[[p]]), "polr"))) {
        stop(paste0("Fit mreg[[", p, "]] by lm, glm, glm.nb, gam, multinom, polr"))
      } 
    }
  } else {
    if (!is.null(mreg)) warning("mreg is ignored when model is 'wb', 'iorw' or 'ne'")
    mreg <- NULL
  }
  
  # postc Regression
  # for gformula, a regression is required for each L
  if (model == "gformula" && length(postc) > 0) {
    if (!is.list(postcreg)) stop("postcreg should be a list")
    if (length(postcreg) != length(postc)) stop("length(postcreg) != length(postc)")
    for (p in 1:length(postcreg)) {
      if (is.null(postcreg[[p]])) stop(paste0("Unspecified postcreg[[", p, "]]"))
      if (is.character(postcreg[[p]])) {
        if (postcreg[[p]] == "loglinear" && length(unique(na.omit(data[, postc[p]]))) != 2) stop(paste0("When postcreg[[", p, "]] is 'loglinear', postc[[", p, "]] should be binary"))
        if (!postcreg[[p]] %in% c("linear", "logistic", "loglinear", "poisson", "quasipoisson",
                                  "negbin", "multinomial", "ordinal"))  stop(
                                    paste0("Select character postcreg[[", p, "]] from 'linear', 'logistic',
                                           'loglinear', 'poisson', 'quasipoisson', 'negbin', 'multinomial', 'ordinal'"))
        # regress each L on A, C
        postc_formula <- paste0(postc[p], "~", paste(c(exposure, basec), collapse = "+"))
        switch(postcreg[[p]],
               linear = postcreg[[p]] <- eval(bquote(glm(.(as.formula(postc_formula)), family = gaussian(), data = .(data)))),
               logistic = postcreg[[p]] <- eval(bquote(glm(.(as.formula(postc_formula)), family = binomial(), data = .(data)))),
               loglinear = postcreg[[p]] <- eval(bquote(glm(.(as.formula(postc_formula)), family = poisson(), data = .(data)))),
               poisson = postcreg[[p]]  <- eval(bquote(glm(.(as.formula(postc_formula)), family = poisson(), data = .(data)))),
               quasipoisson = postcreg[[p]] <- eval(bquote(glm(.(as.formula(postc_formula)), family = quasipoisson(), data = .(data)))),
               negbin = postcreg[[p]] <- eval(bquote(MASS::glm.nb(.(as.formula(postc_formula)), data = .(data)))),
               multinomial = postcreg[[p]] <- eval(bquote(nnet::multinom(.(as.formula(postc_formula)), data = .(data), trace = FALSE))),
               ordinal = postcreg[[p]] <- eval(bquote(MASS::polr(.(as.formula(postc_formula)), data = .(data)))))
      } else if (!(identical(class(postcreg[[p]]), "lm") | identical(class(postcreg[[p]]), c("glm", "lm")) |
                   identical(class(postcreg[[p]]), c("negbin", "glm", "lm")) | identical(class(postcreg[[p]]), c("multinom", "nnet")) |
                   identical(class(postcreg[[p]]), c("gam", "glm", "lm")) | identical(class(postcreg[[p]]), "polr"))) {
        stop(paste0("Fit postcreg[[", p, "]] by lm, glm, glm.nb, gam, multinom, polr"))
      } 
    }
  } else {
    if (!is.null(postcreg)) warning("postcreg is ignored when model is not 'gformula' and when length(postc) = 0")
    postcreg <- NULL
  }
  
  # Outcome Regression
  if (is.null(yreg)) stop("yreg is required")
  if (is.character(yreg)) {
    if (yreg == "loglinear" && length(unique(na.omit(data[, outcome]))) != 2) stop("When yreg is 'loglinear', outcome should be binary")
    if (!yreg %in% c("linear", "logistic", "loglinear", "poisson", "quasipoisson",
                     "negbin", "multinomial", "ordinal", "coxph", "aft_exp",
                     "aft_weibull")) stop(
                       paste0("Select character yreg from 'linear', 'logistic',
                              'loglinear', 'poisson', 'quasipoisson', 'negbin', 'multinomial', 'ordinal',
                              'coxph', 'aft_exp', 'aft_weibull'"))
    if (estimation == "paramfunc" && yreg %in% c("logistic", "coxph")) warning("When estimation is 'paramfunc' and yreg is 'logistic' or 'coxph', the outcome must be rare; ignore this warning if the outcome is rare")
    if (model != "iorw") {
      out$variables$EMint <- EMint
      int.terms <- switch(EMint + 1, "1" = NULL, "2" = paste(exposure, mediator, sep = "*"))
    }
    # for rb, wb and ne, regress Y on A, M and C
    # for iorw, regress Y on A and C
    # for msm, regress Y on A and M
    # for gformula, regress Y on A, M, C and L
    switch(model,
           rb = outcome_formula <- paste0(outcome, "~", paste(c(exposure, mediator, int.terms, basec), collapse = "+")),
           wb = outcome_formula <- paste0(outcome, "~", paste(c(exposure, mediator, int.terms, basec), collapse = "+")),
           ne = outcome_formula <- paste0(outcome, "~", paste(c(exposure, mediator, int.terms, basec), collapse = "+")),
           iorw = outcome_formula <- paste0(outcome, "~", paste(c(exposure, basec), collapse = "+")),
           msm = outcome_formula <- paste0(outcome, "~", paste(c(exposure, mediator, int.terms), collapse = "+")),
           gformula = outcome_formula <- paste0(outcome, "~", paste(c(exposure, mediator, int.terms, basec, postc), collapse = "+")))
    if (yreg %in% c("coxph","aft_exp","aft_weibull")) {
      if (!is.null(event)) {
        outcome_formula <- paste(paste0("Surv(", outcome, ", ", event, ")"),
                                 strsplit(outcome_formula, split = "~")[[1]][2], sep = " ~ ")
      } else outcome_formula <- paste(paste0("Surv(", outcome, ")"),
                                      strsplit(outcome_formula, split = "~")[[1]][2], sep = " ~ ")
    }
    switch(yreg,
           linear = yreg <- eval(bquote(glm(formula = .(as.formula(outcome_formula)),
                                            family = gaussian(), data = .(data)))),
           logistic = yreg <- eval(bquote(glm(formula = .(as.formula(outcome_formula)),
                                              family = binomial(), data = .(data)))),
           loglinear = yreg <- eval(bquote(glm(formula = .(as.formula(outcome_formula)),
                                               family = poisson(), data = .(data)))),
           poisson = yreg  <- eval(bquote(glm(formula = .(as.formula(outcome_formula)),
                                              family = poisson(), data = .(data)))),
           quasipoisson = yreg <- eval(bquote(glm(formula = .(as.formula(outcome_formula)),
                                                  family = quasipoisson(), data = .(data)))),
           negbin = yreg <- eval(bquote(MASS::glm.nb(formula = .(as.formula(outcome_formula)),
                                                     data = .(data)))),
           multinomial = yreg <- eval(bquote(nnet::multinom(formula = .(as.formula(outcome_formula)),
                                                            data = .(data), trace = FALSE))),
           ordinal = yreg <- eval(bquote(MASS::polr(formula = .(as.formula(outcome_formula)),
                                                    data = .(data)))),
           coxph = yreg <- eval(bquote(survival::coxph(formula = .(as.formula(outcome_formula)),
                                                       data = .(data)))),
           aft_exp = yreg <- eval(bquote(survival::survreg(formula = .(as.formula(outcome_formula)),
                                                           dist = "exponential", data = .(data)))),
           aft_weibull = yreg <- eval(bquote(survival::survreg(formula = .(as.formula(outcome_formula)),
                                                               dist = "weibull", data = .(data)))))
  } else if (!(identical(class(yreg), "lm") | identical(class(yreg), c("glm", "lm")) |
               identical(class(yreg), c("negbin", "glm", "lm")) | identical(class(yreg), c("multinom", "nnet")) |
               identical(class(yreg), c("gam", "glm", "lm")) | identical(class(yreg), "polr") |
               identical(class(yreg), "coxph") | identical(class(yreg), "survreg"))) {
    stop("Fit yreg by lm, glm, glm.nb, gam, multinom, polr, coxph, survreg")
  } else {
    if (model %in% c("rb", "gformula")) warning("When model is 'rb' or 'gformula', make sure there is no mediator-mediator interaction in yreg; ignore this warning if there isn't")
  }
  
  # for delta method inference, use survey regressions for yreg and mreg when weights are applied
  if (inference == "delta" && casecontrol && !is.null(yprevalence)) {
    yreg <- suppressWarnings(eval(bquote(survey::svyglm(formula = .(formula(yreg)), family = .(family(yreg)),
                                                        design = survey::svydesign(~ 1, data = .(data))))))
  }
  if (inference == "delta" && casecontrol && !is.null(yprevalence)) {
    if (inherits(mreg[[1]], "glm")) mreg[[1]] <- suppressWarnings(eval(bquote(survey::svyglm(formula = .(formula(mreg[[1]])), 
                                                                                             family = .(family(mreg[[1]])),
                                                                                             design = survey::svydesign(~ 1, data = .(data))))))
    if (inherits(mreg[[1]], "multinom")) mreg[[1]] <- suppressWarnings(eval(bquote(svymultinom(formula = .(formula(mreg[[1]])), 
                                                                                               data = .(data)))))
  }
  out <- list(yreg = yreg, ereg = ereg, mreg = mreg, wmnomreg = wmnomreg, 
              wmdenomreg = wmdenomreg, postcreg = postcreg)
  return(out)
}

estinf <- function() {
  # restrict data types of variables
  allvar <- c(outcome, exposure, mediator, postc, basec)
  for (i in 1:length(allvar))
    if (!(is.numeric(data[, allvar[i]]) | is.factor(data[, allvar[i]]) |
          is.character(data[, allvar[i]]))) stop(paste0("The variable ", allvar[i], " should be numeric, factor or character"))
  # output list
  out <- list()
  # obtain calls, weights, classes, and families of regs
  for (reg_name in c("yreg", "ereg", "mreg", "wmnomreg", "wmdenomreg", "postcreg")) {
    reg <- get(reg_name)
    if (!is.null(reg)) {
      if (reg_name %in% c("yreg", "ereg")) {
        assign(paste0("call_", reg_name), getCall(reg))
        assign("reg_mid", switch((inherits(reg, "rcreg") | inherits(reg, "simexreg")) + 1, "1" = reg, "2" = reg$NAIVEreg))
        assign(paste0("is_lm_", reg_name), inherits(reg_mid, "lm"))
        assign(paste0("is_glm_", reg_name), inherits(reg_mid, "glm"))
        assign(paste0("is_svyglm_", reg_name), inherits(reg_mid, "svyglm"))
        assign(paste0("is_gam_", reg_name), inherits(reg_mid, "gam"))
        if (get(paste0("is_lm_", reg_name)) | get(paste0("is_glm_", reg_name))) assign(paste0("family_", reg_name), family(reg_mid))
        assign(paste0("is_multinom_", reg_name), inherits(reg_mid, "multinom"))
        assign(paste0("is_svymultinom_", reg_name), inherits(reg_mid, "svymultinom"))
        assign(paste0("is_polr_", reg_name), inherits(reg_mid, "polr"))
        if (reg_name == "yreg") {
          assign(paste0("is_survreg_", reg_name), inherits(reg_mid, "survreg"))
          assign(paste0("is_coxph_", reg_name), inherits(reg_mid, "coxph"))
        }
        if (get(paste0("is_svyglm_", reg_name))) {assign(paste0("weights_", reg_name), get(reg_name)$data$.survey.prob.weights)
        } else assign(paste0("weights_", reg_name), model.frame(get(reg_name))$'(weights)')  
      } else {
        assign(paste0("call_", reg_name), lapply(1:length(reg), function(x) getCall(reg[[x]])))
        assign("reg_mid", lapply(1:length(reg), function(x)
          switch((inherits(reg[[x]], "rcreg") | inherits(reg[[x]], "simexreg")) + 1, "1" = reg[[x]], "2" = reg[[x]]$NAIVEreg)))
        assign(paste0("is_lm_", reg_name), sapply(1:length(reg_mid), function(x) inherits(reg_mid[[x]], "lm")))
        assign(paste0("is_glm_", reg_name), sapply(1:length(reg_mid), function(x) inherits(reg_mid[[x]], "glm")))
        assign(paste0("is_svyglm_", reg_name), sapply(1:length(reg_mid), function(x) inherits(reg_mid[[x]], "svyglm")))
        assign(paste0("is_gam_", reg_name), sapply(1:length(reg_mid), function(x) inherits(reg_mid[[x]], "gam")))
        assign(paste0("family_", reg_name), lapply(1:length(reg_mid), function(x)
          if (get(paste0("is_lm_", reg_name))[x] | get(paste0("is_glm_", reg_name))[x]) family(reg_mid[[x]])))
        assign(paste0("is_multinom_", reg_name), sapply(1:length(reg_mid), function(x) inherits(reg_mid[[x]], "multinom")))
        assign(paste0("is_svymultinom_", reg_name), sapply(1:length(reg_mid), function(x) inherits(reg_mid[[x]], "svymultinom")))
        assign(paste0("is_polr_", reg_name), sapply(1:length(reg_mid), function(x) inherits(reg_mid[[x]], "polr")))
        assign(paste0("weights_", reg_name), lapply(1:length(reg_mid), function(x) {
          if (get(paste0("is_svyglm_", reg_name))[x]) { get(reg_name)[[x]]$data$.survey.prob.weights
          } else model.frame(get(reg_name)[[x]])$'(weights)'
        }))
      }
      rm(reg_mid)
    }
  }
  
  # restrict formulas, classes and families of regression objects
  if (!(((is_lm_yreg | is_glm_yreg) && 
         (family_yreg$family %in% 
          c("gaussian", "inverse.gaussian", "quasi", "poisson", "quasipoisson", 
            "Gamma", "binomial", "quasibinomial", "multinom", "ziplss") |
          startsWith(family_yreg$family, "Negative Binomial") |
          startsWith(family_yreg$family, "Zero inflated Poisson") |
          startsWith(family_yreg$family, "Ordered Categorical"))) |
        is_multinom_yreg | is_polr_yreg | is_survreg_yreg | is_coxph_yreg |
        inference == "delta")) stop("Unsupported yreg")
  yreg_formula <- formula(yreg)
  d_var <- unique(all.vars(yreg_formula[[2]]))
  ind_var <- unique(all.vars(yreg_formula[[3]]))
  if (model %in% c("rb", "wb", "ne") && (!(outcome %in% d_var) | !all(ind_var %in% c(exposure, mediator, basec)))) stop(
    "For yreg, please regress outcome on variables in c(exposure, mediator, basec) when model is rb or wb or ne")
  if (model == "iorw" && (!(outcome %in% d_var) | !all(ind_var %in% c(exposure, basec)))) stop(
    "For yreg, please regress outcome on variables in c(exposure, basec) when model is iorw")
  if (model == "msm" && (!(outcome %in% d_var) | !all(ind_var %in% c(exposure, mediator)))) stop(
    "For yreg, please regress outcome on variables in c(exposure, mediator) when model is msm")
  if (model == "gformula" && (!(outcome %in% d_var) | !all(ind_var %in% c(exposure, mediator, basec, postc)))) stop(
    "For yreg, please regress outcome on variables in c(exposure, mediator, basec, postc) when model is gformula")
  rm(yreg_formula, d_var, ind_var)
  if (!is.null(ereg)) {
    if (!(((is_lm_ereg | is_glm_ereg) && 
           (family_ereg$family %in% c("binomial", "quasibinomial", "multinom") |
            startsWith(family_ereg$family, "Ordered Categorical"))) |
          is_multinom_ereg | is_polr_ereg)) stop("Unsupported ereg")
    ereg_formula <- formula(ereg)
    d_var <- unique(all.vars(ereg_formula[[2]]))
    ind_var <- unique(all.vars(ereg_formula[[3]]))
    if (model != "iorw" && ((exposure != d_var) | !all(ind_var %in% basec))) stop("For ereg, please regress the exposure on variables in basec when model is wb or msm")
    if (model == "iorw" && ((exposure != d_var) | !all(mediator %in% ind_var) | 
                            !all(ind_var %in% c(mediator, basec)))) stop("For ereg, please regress the exposure on variables in basec and all mediators when model is iorw")
    rm(ereg_formula, d_var, ind_var)
  }
  if (!is.null(mreg) && inference == "bootstrap") {
    for (p in 1:length(mreg)) {
      if (!(((is_lm_mreg[[p]] | is_glm_mreg[[p]]) && 
             (family_mreg[[p]]$family %in% 
              c("gaussian", "inverse.gaussian", "poisson", "quasipoisson", 
                "Gamma", "binomial", "multinom") |
              startsWith(family_mreg[[p]]$family, "Negative Binomial") |
              startsWith(family_mreg[[p]]$family, "Ordered Categorical"))) |
            is_multinom_mreg[[p]] | is_polr_mreg[[p]])) stop(paste0("Unsupported mreg[[", p, "]]"))
      mreg_formula <- formula(mreg[[p]])
      d_var <- unique(all.vars(mreg_formula[[2]]))
      ind_var <- unique(all.vars(mreg_formula[[3]]))
      if (model == "rb" && ((mediator[p] != d_var) | !all(ind_var %in% c(exposure, basec)))) stop(
        paste0("For mreg[[", p, "]], please regress mediator[", p, "] on variables in c(exposure, basec) when model is rb"))
      if (model == "msm" && ((mediator[p] != d_var) | !all(ind_var %in% c(exposure)))) stop(
        paste0("For mreg[[", p, "]], please regress mediator[", p, "] on exposure when model is msm"))
      if (model == "gformula" && ((mediator[p] != d_var) | !all(ind_var %in% c(exposure, basec, postc)))) stop(
        paste0("For mreg[[", p, "]], please regress mediator[", p, "] on variables in c(exposure, basec, postc) when model is gformula"))
      rm(mreg_formula, d_var, ind_var)
    }
  }
  if (!is.null(wmnomreg)) {
    for (p in 1:length(wmnomreg)) {
      if (!((((is_lm_wmnomreg[[p]] | is_glm_wmnomreg[[p]]) && 
              (family_wmnomreg[[p]]$family %in% 
               c("binomial", "quasibinomial", "multinom") |
               startsWith(family_wmnomreg[[p]]$family, "Ordered Categorical"))) |
             is_multinom_wmnomreg[[p]] | is_polr_wmnomreg[[p]]))) stop(paste0("Unsupported wmnomreg[[", p, "]]"))
      wmnomreg_formula <- formula(wmnomreg[[p]])
      d_var <- unique(all.vars(wmnomreg_formula[[2]]))
      ind_var <- unique(all.vars(wmnomreg_formula[[3]]))
      if ((mediator[p] != d_var) | !all(ind_var %in% c(exposure, mediator[0:(p-1)]))) stop(
        paste0("For wmnomreg[[", p, "]], please regress mediator[", p, "] on variables in c(exposure, mediator[0:", p-1, "])"))
      rm(wmnomreg_formula, d_var, ind_var)
    }
  }
  if (!is.null(wmdenomreg)) {
    for (p in 1:length(wmdenomreg)) {
      if (!((((is_lm_wmdenomreg[[p]] | is_glm_wmdenomreg[[p]]) && 
              (family_wmdenomreg[[p]]$family %in% 
               c("binomial", "quasibinomial", "multinom") |
               startsWith(family_wmdenomreg[[p]]$family, "Ordered Categorical"))) |
             is_multinom_wmdenomreg[[p]] | is_polr_wmdenomreg[[p]]))) stop(paste0("Unsupported wmdenomreg[[", p, "]]"))
      wmdenomreg_formula <- formula(wmdenomreg[[p]])
      d_var <- unique(all.vars(wmdenomreg_formula[[2]]))
      ind_var <- unique(all.vars(wmdenomreg_formula[[3]]))
      if ((mediator[p] != d_var) | !all(ind_var %in% c(exposure, mediator[0:(p-1)], basec, postc))) stop(
        paste0("For wmdenomreg[[", p, "]], please regress mediator[", p, "] on variables in c(exposure, mediator[0:", p-1, "], basec, postc)"))
      rm(wmdenomreg_formula, d_var, ind_var)
    }
  }
  if (!is.null(postcreg)) {
    for (p in 1:length(postcreg)) {
      if (!(((is_lm_postcreg[[p]] | is_glm_postcreg[[p]]) && 
             (family_postcreg[[p]]$family %in% 
              c("gaussian", "inverse.gaussian", "poisson", "quasipoisson", 
                "Gamma", "binomial", "multinom") |
              startsWith(family_postcreg[[p]]$family, "Negative Binomial") |
              startsWith(family_postcreg[[p]]$family, "Ordered Categorical"))) |
            is_multinom_postcreg[[p]] | is_polr_postcreg[[p]])) stop(paste0("Unsupported postcreg[[", p, "]]"))
      postcreg_formula <- formula(postcreg[[p]])
      d_var <- unique(all.vars(postcreg_formula[[2]]))
      ind_var <- unique(all.vars(postcreg_formula[[3]]))
      if ((postc[p] != d_var) | !all(ind_var %in% c(exposure, basec))) stop(
        paste0("For postcreg[[", p, "]], please regress postc[", p, "] on variables in c(exposure, basec)"))
      rm(postcreg_formula, d_var, ind_var)
    }
  }
  
  # reference values for the exposure
  if (is.factor(data[, exposure]) | is.character(data[, exposure])) {
    a_lev <- levels(droplevels(as.factor(data[, exposure])))
    if (!a %in% a_lev) {
      a <- a_lev[length(a_lev)]
      warning(paste0("a is not a value of the exposure; ", a, " is used"))
    }
    if (!astar %in% a_lev) {
      astar <- a_lev[1]
      warning(paste0("astar is not a value of the exposure; ", astar, " is used"))
    }
  }
  out$ref$a <- a
  out$ref$astar <- astar
  
  # yval: the reference level for a categorical outcome
  if ((is_glm_yreg && (family_yreg$family %in% c("binomial", "quasibinomial", "multinom") |
                       startsWith(family_yreg$family, "Ordered Categorical"))) |
      is_multinom_yreg | is_polr_yreg) {
    y_lev <- levels(droplevels(as.factor(data[, outcome])))
    # if yval is not provided or yval provided is not a value of the outcome, use the last level of the outcome
    if (is.null(yval)) {
      yval <- y_lev[length(y_lev)]
      warning(paste0("yval is not specified; ", yval, " is used"))
    }
    if (!yval %in% y_lev) {
      yval <- y_lev[length(y_lev)]
      warning(paste0("yval is not a value of the outcome; ", yval, " is used"))
    }
    out$ref$yval <- yval
  }
  
  # reference values for the mediators
  if (model != "iorw") {
    if (model %in% c("wb", "ne")) {
      for (i in 1:length(mediator)) {
        if (is.factor(data[, mediator[i]]) | is.character(data[, mediator[i]])) {
          m_lev <- levels(droplevels(as.factor(data[, mediator[i]])))
          if (!mval[[i]] %in% m_lev) {
            mval[[i]] <- m_lev[length(m_lev)]
            warning(paste0("mval[[", i, "]] is not a value of mediator[", i, "]; ", mval[[i]], " is used"))
          }
        }
      }
    } else {
      for (i in 1:length(mediator)) {
        if ((is_glm_mreg[i] && ((family_mreg[[i]]$family %in% c("binomial", "multinom")) |
                                startsWith(family_mreg[[i]]$family, "Ordered Categorical")))|
            is_multinom_mreg[i] | is_polr_mreg[i]) {
          m_lev <- levels(droplevels(as.factor(data[, mediator[i]])))
          if (is.numeric(data[, mediator[i]])) m_lev <- as.numeric(m_lev)
          if (!mval[[i]] %in% m_lev) {
            mval[[i]] <- m_lev[length(m_lev)]
            warning(paste0("mval[[", i, "]] is not a value of mediator[", i, "]; ", mval[[i]], " is used"))
          }
        }
      }
    }
    out$ref$mval <- mval
  }
  
  # get the level of the case and the level of the control
  if (casecontrol) {
    y_lev <- levels(droplevels(as.factor(data[, outcome])))
    if (length(y_lev) != 2) stop("outcome with more than 2 levels")
    y_control <- y_lev[1]
    y_case <- y_lev[2]
  }
  
  if (model == "rb") {
    ###################################################################################################
    ############################################Regression-based Approach##############################
    ###################################################################################################
    # closed-form parameter function estimation
    if (estimation == "paramfunc") {
      # create a list of covariate values to calculate conditional causal effects
      if (length(basec) != 0) {
        if (!is.null(basecval)) {
          if (!is.list(basecval)) stop("basecval should be a list")
          if (length(basecval) != length(basec)) stop("length(basecval) != length(basec)")
        }
        if (is.null(basecval)) basecval <- rep(list(NULL), length(basec))
        # if NULL, set basecval[[i]] to be the mean value of basec[i]
        for (i in 1:length(basec)) {
          if (is.factor(data[, basec[i]]) | is.character(data[, basec[i]])) {
            c_lev <- levels(droplevels(as.factor(data[, basec[i]])))
            if (is.null(basecval[[i]])) {
              c_data <- data[, basec[i], drop = FALSE]
              c_data[, basec[i]] <- factor(c_data[, basec[i]], levels = c_lev)
              # set basecval[[i]] to be the mean values of dummy variables
              basecval[[i]] <- unname(colMeans(as.matrix(model.matrix(as.formula(paste0("~", basec[i])),
                                                                      data = model.frame(~., data = c_data, 
                                                                                         na.action = na.pass))[, -1]), 
                                               na.rm = TRUE))
              rm(c_data)
              # dummy values of basecval[[i]]
            } else basecval[[i]] <- as.numeric(c_lev == basecval[[i]])[-1]
            rm(c_lev)
          } else if (is.numeric(data[, basec[i]])) {
            if (is.null(basecval[[i]])) {
              # set basecval[[i]] to be the mean value of basec[i]
              basecval[[i]] <- mean(data[, basec[i]], na.rm = TRUE)
            } else basecval[[i]] <- basecval[[i]]
          } 
        }
        out$ref$basecval <- basecval
      }
    }
    
    # estimation and inference
    environment(est.rb) <- environment()
    if (!multimp) {
      # point estimates of causal effects
      est <- est.rb(data = data, indices = NULL, outReg = TRUE, full = full)
      effect.pe <- est$est
      n_effect <- length(effect.pe)
      out$reg.output <- est$reg.output
      if (inference == "bootstrap") {
        # bootstrap results
        boots <- boot(data = data, statistic = est.rb, R = nboot, outReg = FALSE, full = full)
        # bootstrap CIs
        environment(boot.ci) <- environment()
        effect.ci <- boot.ci(boots = boots)
        effect.ci.low <- effect.ci[, 1]
        effect.ci.high <- effect.ci[, 2]
        # bootstrap p-values
        effect.pval <- sapply(1:n_effect, function(x) boot.pval(boots = boots$t[, x], pe = effect.pe[x]))
      } else if (inference == "delta") {
        yreg <- est$reg.output$yreg
        mreg <- est$reg.output$mreg[[1]]
        # standard errors by the delta method
        environment(inf.delta) <- environment()
        effect.se <- inf.delta(data = data, yreg = yreg, mreg = mreg)
        # critical value
        z0 <- qnorm(0.975)
        z <- effect.pe/effect.se
        # delta method CIs
        effect.ci.low <- effect.pe - z0 * effect.se
        effect.ci.high <- effect.pe + z0 * effect.se
        # delta method p-values
        effect.pval <- 2 * (1 - pnorm(abs(z)))
      }
    } else {
      # imputed data sets
      data_imp <- complete(do.call(mice, args_mice), action = "all")
      m <- length(data_imp)
      # estimate causal effects for each imputed data set
      est_imp <- lapply(1:m, function(x)
        est.rb(data = data_imp[[x]], indices = NULL, outReg = TRUE, full = full))
      est_imp_df <- do.call(rbind, lapply(1:m, function(x) est_imp[[x]]$est))
      effect.pe <- colMeans(est_imp_df)
      n_effect <- length(effect.pe)
      out$reg.output <- lapply(1:m, function(x) est_imp[[x]]$reg.output)
      
      if (inference == "bootstrap") {
        boot.step <- function(data = NULL, indices = NULL) {
          data_boot <- data[indices, ]
          args_mice$data <- data_boot
          data_imp <- complete(do.call(mice, args_mice), action = "all")
          curVal <- get("counter", envir = env)
          assign("counter", curVal + 1, envir = env)
          setTxtProgressBar(get("progbar", envir = env), curVal + 1)
          return(colMeans(do.call(rbind, lapply(1:m, function(x)
            est.rb(data = data_imp[[x]], outReg = FALSE, full = full)))))
        }
        environment(boot.step) <- environment()
        # bootstrap results
        boots <- boot(data = data, statistic = boot.step, R = nboot)
        # bootstrap CIs
        environment(boot.ci) <- environment()
        effect.ci <- boot.ci(boots = boots)
        effect.ci.low <- effect.ci[, 1]
        effect.ci.high <- effect.ci[, 2]
        # bootstrap p-values
        effect.pval <- sapply(1:n_effect, function(x) boot.pval(boots = boots$t[, x], pe = effect.pe[x]))
      } else if (inference == "delta") {
        environment(inf.delta) <- environment()
        # standard errors by the delta method
        se_imp <- do.call(rbind, lapply(1:m, function(x)
          inf.delta(data = data_imp[[x]], yreg = est_imp[[x]]$reg.output$yreg, mreg = est_imp[[x]]$reg.output$mreg[[1]])))
        # pool the results by Rubin's rule
        var_within <- colMeans(se_imp ^ 2)
        var_between <- colSums((est_imp_df - t(replicate(m, effect.pe)))^2)/(m - 1)
        effect.se <- sqrt(var_within + var_between * (m + 1) / m)
        z0 <- qnorm(0.975)
        z <- effect.pe/effect.se
        effect.ci.low <- effect.pe - z0 * effect.se
        effect.ci.high <- effect.pe + z0 * effect.se
        effect.pval <- 2 * (1 - pnorm(abs(z)))
      }
    }
    
    if ((is_lm_yreg | is_glm_yreg) &&
        (family_yreg$family %in% c("gaussian", "inverse.gaussian", "Gamma", "quasi"))) {
      # standard errors by bootstrapping
      if (inference == "bootstrap") effect.se <- sapply(1:n_effect, function(x) sd(boots$t[, x], na.rm = TRUE))
      # effect names
      if (full) {
        if (EMint) effect_name <- c("cde", "pnde", "tnde", "pnie", "tnie", "te", 
                                    "intref", "intmed", "cde(prop)", "intref(prop)", "intmed(prop)", "pnie(prop)",
                                    "pm", "int", "pe")
        if (!EMint) effect_name <- c("cde", "pnde", "tnde", "pnie", "tnie", "te", "pm")
      }
      if (!full) effect_name <- c("cde", "pnde", "tnde", "pnie", "tnie", "te")
    } else {
      # transform standard errors of effects in log scale
      if (inference == "bootstrap") effect.se <- sapply(1:n_effect, function(x)
        ifelse(x <= 6, sd(exp(boots$t[, x])), sd(boots$t[, x]))) #
      if (inference == "delta") effect.se[1:6] <- sapply(1:6, function(x)
        deltamethod(as.formula("~exp(x1)"), effect.pe[x], effect.se[x]^2))
      # transform effects in log ratio scale into effects in ratio scale
      effect.pe[1:6] <- exp(effect.pe[1:6])
      effect.ci.low[1:6] <- exp(effect.ci.low[1:6])
      effect.ci.high[1:6] <- exp(effect.ci.high[1:6])
      # effect names
      if (full) {
        if (EMint) effect_name <- c("Rcde", "Rpnde", "Rtnde", "Rpnie", "Rtnie", "Rte", 
                                    "ERcde", "ERintref", "ERintmed", "ERpnie",
                                    "ERcde(prop)", "ERintref(prop)", "ERintmed(prop)", "ERpnie(prop)",
                                    "pm", "int", "pe")
        if (!EMint) effect_name <- c("Rcde", "Rpnde", "Rtnde", "Rpnie", "Rtnie", "Rte", "pm")
      }
      if (!full) effect_name <- c("Rcde", "Rpnde", "Rtnde", "Rpnie", "Rtnie", "Rte")
    }
    
    names(effect.pe) <- names(effect.se) <- names(effect.ci.low) <- names(effect.ci.high) <-
      names(effect.pval) <- effect_name
    out$effect.pe <- effect.pe
    out$effect.se <- effect.se
    out$effect.ci.low <- effect.ci.low
    out$effect.ci.high <- effect.ci.high
    out$effect.pval <- effect.pval
    
  } else if (model == "gformula") {
    ###################################################################################################
    #############################################G-formula Approach####################################
    ###################################################################################################
    environment(est.gformula) <- environment()
    if (!multimp) {
      # point estimates of causal effects
      est <- est.gformula(data = data, indices = NULL, outReg = TRUE, full = full)
      effect.pe <- est$est
      n_effect <- length(effect.pe)
      out$reg.output <- est$reg.output
      # bootstrap results
      boots <- boot(data = data, statistic = est.gformula, R = nboot, outReg = FALSE, full = full)
      # bootstrap CIs
      environment(boot.ci) <- environment()
      effect.ci <- boot.ci(boots = boots)
      effect.ci.low <- effect.ci[, 1]
      effect.ci.high <- effect.ci[, 2]
      # bootstrap p-values
      effect.pval <- sapply(1:n_effect, function(x) boot.pval(boots = boots$t[, x], pe = effect.pe[x]))
    } else {
      # imputed data sets
      data_imp <- complete(do.call(mice, args_mice), action = "all")
      m <- length(data_imp)
      # estimate causal effects for each imputed data set
      est_imp <- lapply(1:m, function(x)
        est.gformula(data = data_imp[[x]], indices = NULL, outReg = TRUE, full = full))
      est_imp_df <- do.call(rbind, lapply(1:m, function(x) est_imp[[x]]$est))
      effect.pe <- colMeans(est_imp_df)
      n_effect <- length(effect.pe)
      out$reg.output <- lapply(1:m, function(x) est_imp[[x]]$reg.output)
      
      boot.step <- function(data = NULL, indices = NULL) {
        data_boot <- data[indices, ]
        args_mice$data <- data_boot
        data_imp <- complete(do.call(mice, args_mice), action = "all")
        curVal <- get("counter", envir = env)
        assign("counter", curVal + 1, envir = env)
        setTxtProgressBar(get("progbar", envir = env), curVal + 1)
        return(colMeans(do.call(rbind, lapply(1:m, function(x)
          est.gformula(data = data_imp[[x]], outReg = FALSE, full = full)))))
      }
      environment(boot.step) <- environment()
      # bootstrap results
      boots <- boot(data = data, statistic = boot.step, R = nboot)
      # bootstrap CIs
      environment(boot.ci) <- environment()
      effect.ci <- boot.ci(boots = boots)
      effect.ci.low <- effect.ci[, 1]
      effect.ci.high <- effect.ci[, 2]
      # bootstrap p-values
      effect.pval <- sapply(1:n_effect, function(x) boot.pval(boots = boots$t[, x], pe = effect.pe[x]))
    }
    
    if ((is_lm_yreg | is_glm_yreg) &&
        (family_yreg$family %in% c("gaussian", "inverse.gaussian", "Gamma", "quasi"))) {
      # standard errors by bootstrapping
      effect.se <- sapply(1:n_effect, function(x) sd(boots$t[, x]))
      # effect names
      if (length(postc) == 0) {
        if (full) {
          if (EMint) effect_name <-
              c("cde", "pnde", "tnde", "pnie", "tnie", "te", "intref", "intmed", "cde(prop)", 
                "intref(prop)", "intmed(prop)", "pnie(prop)", "pm", "int", "pe")
          if (!EMint) effect_name <- c("cde", "pnde", "tnde", "pnie", "tnie", "te", "pm")
        } else effect_name <- c("cde", "pnde", "tnde", "pnie", "tnie", "te")
      } else {
        if (full) {
          if (EMint) effect_name <-
              c("cde", "rpnde", "rtnde", "rpnie", "rtnie", "te", "rintref", "rintmed", "cde(prop)", 
                "rintref(prop)", "rintmed(prop)", "rpnie(prop)", "rpm", "rint", "rpe")
          if (!EMint) effect_name <- c("cde", "rpnde", "rtnde", "rpnie", "rtnie", "te", "pm")
        } else effect_name <- c("cde", "rpnde", "rtnde", "rpnie", "rtnie", "te")
      }
    } else {
      # transform standard errors of effects on the log scale
      effect.se <- sapply(1:n_effect, function(x) ifelse(x <= 6, sd(exp(boots$t[, x])), sd(boots$t[, x])))
      # transform effects on the log ratio scale into effects on the ratio scale
      effect.pe[1:6] <- exp(effect.pe[1:6])
      effect.ci.low[1:6] <- exp(effect.ci.low[1:6])
      effect.ci.high[1:6] <- exp(effect.ci.high[1:6])
      # effect names
      if (length(postc) == 0) {
        if (full) {
          if (EMint) effect_name <-
              c("Rcde", "Rpnde", "Rtnde", "Rpnie", "Rtnie", "Rte", "ERcde", "ERintref", "ERintmed", "ERpnie",
                "ERcde(prop)", "ERintref(prop)", "ERintmed(prop)", "ERpnie(prop)", "pm", "int", "pe")
          if (!EMint) effect_name <- c("Rcde", "Rpnde", "Rtnde", "Rpnie", "Rtnie", "Rte", "pm")
        } else effect_name <- c("Rcde", "Rpnde", "Rtnde", "Rpnie", "Rtnie", "Rte")
      } else {
        if (full) {
          if (EMint) effect_name <-
              c("Rcde", "rRpnde", "rRtnde", "rRpnie", "rRtnie", "Rte", "ERcde", "rERintref", "rERintmed", 
                "rERpnie", "ERcde(prop)", "rERintref(prop)", "rERintmed(prop)", "rERpnie(prop)", "rpm", "rint", "rpe")
          if (!EMint) effect_name <- c("Rcde", "rRpnde", "rRtnde", "rRpnie", "rRtnie", "Rte", "pm")
        } else effect_name <- c("Rcde", "rRpnde", "rRtnde", "rRpnie", "rRtnie", "Rte")
      }
    }
    names(effect.pe) <- names(effect.se) <- names(effect.ci.low) <- names(effect.ci.high) <-
      names(effect.pval) <- effect_name
    out$effect.pe <- effect.pe
    out$effect.se <- effect.se
    out$effect.ci.low <- effect.ci.low
    out$effect.ci.high <- effect.ci.high
    out$effect.pval <- effect.pval
    
  } else if (model == "wb") {
    ###################################################################################################
    ############################################Weighting-based Approach##############################
    ###################################################################################################
    if (length(basec) != 0 && (!((is_glm_ereg && (family_ereg$family %in% c("binomial", "quasibinomial", "multinom") |
                                                  startsWith(family_ereg$family, "Ordered Categorical"))) |
                                 is_multinom_ereg | is_polr_ereg))) stop(
                                   "model = 'wb' only supports categorical exposure when length(basec) != 0")
    if (is_survreg_yreg | is_coxph_yreg) stop("model = 'wb' doesn't support survival outcomes")
    
    environment(est.wb) <- environment()
    if (!multimp) {
      # point estimates of causal effects
      est <- est.wb(data = data, indices = NULL, outReg = TRUE, full = full)
      effect.pe <- est$est
      n_effect <- length(effect.pe)
      out$reg.output <- est$reg.output
      # bootstrap results
      boots <- boot(data = data, statistic = est.wb, R = nboot, outReg = FALSE, full = full)
      # bootstrap CIs
      environment(boot.ci) <- environment()
      effect.ci <- boot.ci(boots = boots)
      effect.ci.low <- effect.ci[, 1]
      effect.ci.high <- effect.ci[, 2]
      # bootstrap p-values
      effect.pval <- sapply(1:n_effect, function(x) boot.pval(boots = boots$t[, x], pe = effect.pe[x]))
    } else {
      # imputed data sets
      data_imp <- complete(do.call(mice, args_mice), action = "all")
      m <- length(data_imp)
      # estimate causal effects for each imputed data set
      est_imp <- lapply(1:m, function(x)
        est.wb(data = data_imp[[x]], indices = NULL, outReg = TRUE, full = full))
      est_imp_df <- do.call(rbind, lapply(1:m, function(x) est_imp[[x]]$est))
      effect.pe <- colMeans(est_imp_df)
      n_effect <- length(effect.pe)
      out$reg.output <- lapply(1:m, function(x) est_imp[[x]]$reg.output)
      boot.step <- function(data = NULL, indices = NULL) {
        data_boot <- data[indices, ]
        args_mice$data <- data_boot
        data_imp <- complete(do.call(mice, args_mice), action = "all")
        curVal <- get("counter", envir = env)
        assign("counter", curVal + 1, envir = env)
        setTxtProgressBar(get("progbar", envir = env), curVal + 1)
        return(colMeans(do.call(rbind, lapply(1:m, function(x)
          est.wb(data = data_imp[[x]], outReg = FALSE, full = full)))))
      }
      environment(boot.step) <- environment()
      # bootstrap results
      boots <- boot(data = data, statistic = boot.step, R = nboot)
      # bootstrap CIs
      environment(boot.ci) <- environment()
      effect.ci <- boot.ci(boots = boots)
      effect.ci.low <- effect.ci[, 1]
      effect.ci.high <- effect.ci[, 2]
      # bootstrap p-values
      effect.pval <- sapply(1:n_effect, function(x) boot.pval(boots = boots$t[, x], pe = effect.pe[x]))
    }
    
    if ((is_lm_yreg | is_glm_yreg) &&
        (family_yreg$family %in% c("gaussian", "inverse.gaussian", "Gamma", "quasi"))) {
      # standard errors by bootstrapping
      effect.se <- sapply(1:n_effect, function(x) sd(boots$t[, x]))
      # effect names
      if (full) {
        if (EMint) effect_name <- c("cde", "pnde", "tnde", "pnie", "tnie", "te", 
                                    "intref", "intmed", "cde(prop)", "intref(prop)", "intmed(prop)", "pnie(prop)",
                                    "pm", "int", "pe")
        if (!EMint) effect_name <- c("cde", "pnde", "tnde", "pnie", "tnie", "te", "pm")
      }
      if (!full) effect_name <- c("cde", "pnde", "tnde", "pnie", "tnie", "te")
    } else {
      # transform standard errors of effects on the log scale
      effect.se <- sapply(1:n_effect, function(x) ifelse(x <= 6, sd(exp(boots$t[, x])), sd(boots$t[, x])))
      # transform effects on the log ratio scale into effects on the ratio scale
      effect.pe[1:6] <- exp(effect.pe[1:6])
      effect.ci.low[1:6] <- exp(effect.ci.low[1:6])
      effect.ci.high[1:6] <- exp(effect.ci.high[1:6])
      # effect names
      if (full) {
        if (EMint) effect_name <- c("Rcde", "Rpnde", "Rtnde", "Rpnie", "Rtnie", "Rte", 
                                    "ERcde", "ERintref", "ERintmed", "ERpnie",
                                    "ERcde(prop)", "ERintref(prop)", "ERintmed(prop)", "ERpnie(prop)",
                                    "pm", "int", "pe")
        if (!EMint) effect_name <- c("Rcde", "Rpnde", "Rtnde", "Rpnie", "Rtnie", "Rte", "pm")
      }
      if (!full) effect_name <- c("Rcde", "Rpnde", "Rtnde", "Rpnie", "Rtnie", "Rte")
    }
    
    names(effect.pe) <- names(effect.se) <- names(effect.ci.low) <- names(effect.ci.high) <-
      names(effect.pval) <- effect_name
    out$effect.pe <- effect.pe
    out$effect.se <- effect.se
    out$effect.ci.low <- effect.ci.low
    out$effect.ci.high <- effect.ci.high
    out$effect.pval <- effect.pval
    
  } else if (model == "iorw") {
    ###################################################################################################
    ###################################Inverse Odds Ratio Weighting Approach###########################
    ###################################################################################################
    if (!((is_glm_ereg && (family_ereg$family %in% c("binomial", "quasibinomial", "multinom") |
                           startsWith(family_ereg$family, "Ordered Categorical"))) |
          is_multinom_ereg | is_polr_ereg)) stop("model = 'iorw' only supports categorical exposure")
    
    environment(est.iorw) <- environment()
    if (!multimp) {
      # point estimates of causal effects
      est <- est.iorw(data = data, indices = NULL, outReg = TRUE, full = full)
      effect.pe <- est$est
      n_effect <- length(effect.pe)
      out$reg.output <- est$reg.output
      # bootstrap results
      boots <- boot(data = data, statistic = est.iorw, R = nboot, outReg = FALSE, full = full)
      # bootstrap CIs
      environment(boot.ci) <- environment()
      effect.ci <- boot.ci(boots = boots)
      effect.ci.low <- effect.ci[, 1]
      effect.ci.high <- effect.ci[, 2]
      # bootstrap p-values
      effect.pval <- sapply(1:n_effect, function(x) boot.pval(boots = boots$t[, x], pe = effect.pe[x]))
    } else {
      # imputed data sets
      data_imp <- complete(do.call(mice, args_mice), action = "all")
      m <- length(data_imp)
      # estimate causal effects for each imputed data set
      est_imp <- lapply(1:m, function(x)
        est.iorw(data = data_imp[[x]], indices = NULL, outReg = TRUE, full = full))
      est_imp_df <- do.call(rbind, lapply(1:m, function(x) est_imp[[x]]$est))
      effect.pe <- colMeans(est_imp_df)
      n_effect <- length(effect.pe)
      out$reg.output <- lapply(1:m, function(x) est_imp[[x]]$reg.output)
      boot.step <- function(data = NULL, indices = NULL) {
        data_boot <- data[indices, ]
        args_mice$data <- data_boot
        data_imp <- complete(do.call(mice, args_mice), action = "all")
        curVal <- get("counter", envir = env)
        assign("counter", curVal + 1, envir = env)
        setTxtProgressBar(get("progbar", envir = env), curVal + 1)
        return(colMeans(do.call(rbind, lapply(1:m, function(x)
          est.iorw(data = data_imp[[x]], outReg = FALSE, full = full)))))
      }
      environment(boot.step) <- environment()
      # bootstrap results
      boots <- boot(data = data, statistic = boot.step, R = nboot)
      # bootstrap CIs
      environment(boot.ci) <- environment()
      effect.ci <- boot.ci(boots = boots)
      effect.ci.low <- effect.ci[, 1]
      effect.ci.high <- effect.ci[, 2]
      # bootstrap p-values
      effect.pval <- sapply(1:n_effect, function(x) boot.pval(boots = boots$t[, x], pe = effect.pe[x]))
    }
    
    if ((is_lm_yreg | is_glm_yreg) &&
        (family_yreg$family %in% c("gaussian", "inverse.gaussian", "Gamma", "quasi"))) {
      # standard errors by bootstrapping
      effect.se <- sapply(1:n_effect, function(x) sd(boots$t[, x]))
      # effect names
      if (full) effect_name <- c("te", "pnde", "tnie", "pm")
      if (!full) effect_name <- c("te", "pnde", "tnie")
    } else {
      # transform standard errors of effects on the log scale
      effect.se <- sapply(1:n_effect, function(x) ifelse(x <= 3, sd(exp(boots$t[, x])), sd(boots$t[, x])))
      # transform effects on thelog ratio scale into effects on the ratio scale
      effect.pe[1:3] <- exp(effect.pe[1:3])
      effect.ci.low[1:3] <- exp(effect.ci.low[1:3])
      effect.ci.high[1:3] <- exp(effect.ci.high[1:3])
      # effect names
      if (full) effect_name <- c("Rte", "Rpnde", "Rtnie", "pm")
      if (!full) effect_name <- c("Rte", "Rpnde", "Rtnie")
    }
    
    names(effect.pe) <- names(effect.se) <- names(effect.ci.low) <- names(effect.ci.high) <-
      names(effect.pval) <- effect_name
    out$effect.pe <- effect.pe
    out$effect.se <- effect.se
    out$effect.ci.low <- effect.ci.low
    out$effect.ci.high <- effect.ci.high
    out$effect.pval <- effect.pval
    
  } else if (model == "msm") {
    ###################################################################################################
    #########################################Marginal Structural Model#################################
    ###################################################################################################
    if (length(basec) != 0 && (!((is_glm_ereg && (family_ereg$family %in% c("binomial", "quasibinomial", "multinom") |
                                                  startsWith(family_ereg$family, "Ordered Categorical"))) |
                                 is_multinom_ereg | is_polr_ereg))) stop(
                                   "model = 'msm' only supports categorical exposure when length(basec) != 0")
    for (p in 1:length(mediator)) {
      if (!((is_glm_mreg[p] && (family_mreg[[p]]$family %in% c("binomial", "multinom") |
                                startsWith(family_mreg[[p]]$family, "Ordered Categorical"))) |
            is_multinom_mreg[p] | is_polr_mreg[p])) stop(
              "model = 'msm' only supports categorical mediators")
      if (!((is_glm_wmnomreg[p] && (family_wmnomreg[[p]]$family %in% c("binomial", "quasibinomial", "multinom") |
                                    startsWith(family_wmnomreg[[p]]$family, "Ordered Categorical"))) |
            is_multinom_wmnomreg[p] | is_polr_wmnomreg[p])) stop(
              "model = 'msm' only supports categorical mediators")
      if (!((is_glm_wmdenomreg[p] && (family_wmdenomreg[[p]]$family %in% c("binomial", "quasibinomial", "multinom") |
                                      startsWith(family_wmdenomreg[[p]]$family, "Ordered Categorical"))) |
            is_multinom_wmdenomreg[p] | is_polr_wmdenomreg[p])) stop(
              "model = 'msm' only supports categorical mediators")
    }
    
    environment(est.msm) <- environment()
    if (!multimp) {
      # point estimates of causal effects
      est <- est.msm(data = data, indices = NULL, outReg = TRUE, full = full)
      effect.pe <- est$est
      n_effect <- length(effect.pe)
      out$reg.output <- est$reg.output
      # bootstrap results
      boots <- boot(data = data, statistic = est.msm, R = nboot, outReg = FALSE, full = full)
      # bootstrap CIs
      environment(boot.ci) <- environment()
      effect.ci <- boot.ci(boots = boots)
      effect.ci.low <- effect.ci[, 1]
      effect.ci.high <- effect.ci[, 2]
      # bootstrap p-values
      effect.pval <- sapply(1:n_effect, function(x) boot.pval(boots = boots$t[, x], pe = effect.pe[x]))
    } else {
      # imputed data sets
      data_imp <- complete(do.call(mice, args_mice), action = "all")
      m <- length(data_imp)
      # estimate causal effects for each imputed data set
      est_imp <- lapply(1:m, function(x)
        est.msm(data = data_imp[[x]], indices = NULL, outReg = TRUE, full = full))
      est_imp_df <- do.call(rbind, lapply(1:m, function(x) est_imp[[x]]$est))
      effect.pe <- colMeans(est_imp_df)
      n_effect <- length(effect.pe)
      out$reg.output <- lapply(1:m, function(x) est_imp[[x]]$reg.output)
      
      boot.step <- function(data = NULL, indices = NULL) {
        data_boot <- data[indices, ]
        args_mice$data <- data_boot
        data_imp <- complete(do.call(mice, args_mice), action = "all")
        curVal <- get("counter", envir = env)
        assign("counter", curVal + 1, envir = env)
        setTxtProgressBar(get("progbar", envir = env), curVal + 1)
        return(colMeans(do.call(rbind, lapply(1:m, function(x)
          est.msm(data = data_imp[[x]], outReg = FALSE, full = full)))))
      }
      environment(boot.step) <- environment()
      # bootstrap results
      boots <- boot(data = data, statistic = boot.step, R = nboot)
      # bootstrap CIs
      environment(boot.ci) <- environment()
      effect.ci <- boot.ci(boots = boots)
      effect.ci.low <- effect.ci[, 1]
      effect.ci.high <- effect.ci[, 2]
      # bootstrap p-values
      effect.pval <- sapply(1:n_effect, function(x) boot.pval(boots = boots$t[, x], pe = effect.pe[x]))
    }
    
    if ((is_lm_yreg | is_glm_yreg) &&
        (family_yreg$family %in% c("gaussian", "inverse.gaussian", "Gamma", "quasi"))) {
      # standard errors by bootstrapping
      effect.se <- sapply(1:n_effect, function(x) sd(boots$t[, x]))
      # effect names
      if (length(postc) == 0) {
        if (full) {
          if (EMint) effect_name <-
              c("cde", "pnde", "tnde", "pnie", "tnie", "te", "intref", "intmed", "cde(prop)", 
                "intref(prop)", "intmed(prop)", "pnie(prop)", "pm", "int", "pe")
          if (!EMint) effect_name <- c("cde", "pnde", "tnde", "pnie", "tnie", "te", "pm")
        } else effect_name <- c("cde", "pnde", "tnde", "pnie", "tnie", "te")
      } else {
        if (full) {
          if (EMint) effect_name <-
              c("cde", "rpnde", "rtnde", "rpnie", "rtnie", "te", "rintref", "rintmed", "cde(prop)", 
                "rintref(prop)", "rintmed(prop)", "rpnie(prop)", "rpm", "rint", "rpe")
          if (!EMint) effect_name <- c("cde", "rpnde", "rtnde", "rpnie", "rtnie", "te", "pm")
        } else effect_name <- c("cde", "rpnde", "rtnde", "rpnie", "rtnie", "te")
      }
    } else {
      # transform standard errors of effects on the log scale
      effect.se <- sapply(1:n_effect, function(x) ifelse(x <= 6, sd(exp(boots$t[, x])), sd(boots$t[, x])))
      # transform effects on the log ratio scale into effects on the ratio scale
      effect.pe[1:6] <- exp(effect.pe[1:6])
      effect.ci.low[1:6] <- exp(effect.ci.low[1:6])
      effect.ci.high[1:6] <- exp(effect.ci.high[1:6])
      # effect names
      if (length(postc) == 0) {
        if (full) {
          if (EMint) effect_name <-
              c("Rcde", "Rpnde", "Rtnde", "Rpnie", "Rtnie", "Rte", "ERcde", "ERintref", "ERintmed", "ERpnie",
                "ERcde(prop)", "ERintref(prop)", "ERintmed(prop)", "ERpnie(prop)", "pm", "int", "pe")
          if (!EMint) effect_name <- c("Rcde", "Rpnde", "Rtnde", "Rpnie", "Rtnie", "Rte", "pm")
        } else effect_name <- c("Rcde", "Rpnde", "Rtnde", "Rpnie", "Rtnie", "Rte")
      } else {
        if (full) {
          if (EMint) effect_name <-
              c("Rcde", "rRpnde", "rRtnde", "rRpnie", "rRtnie", "Rte", "ERcde", "rERintref", "rERintmed", 
                "rERpnie", "ERcde(prop)", "rERintref(prop)", "rERintmed(prop)", "rERpnie(prop)", "rpm", "rint", "rpe")
          if (!EMint) effect_name <- c("Rcde", "rRpnde", "rRtnde", "rRpnie", "rRtnie", "Rte", "pm")
        } else effect_name <- c("Rcde", "rRpnde", "rRtnde", "rRpnie", "rRtnie", "Rte")
      }
    }
    names(effect.pe) <- names(effect.se) <- names(effect.ci.low) <- names(effect.ci.high) <-
      names(effect.pval) <- effect_name
    out$effect.pe <- effect.pe
    out$effect.se <- effect.se
    out$effect.ci.low <- effect.ci.low
    out$effect.ci.high <- effect.ci.high
    out$effect.pval <- effect.pval
    
  } else if (model == "ne") {
    ###################################################################################################
    #########################################Natural Effect Model######################################
    ###################################################################################################
    if (!identical(class(yreg), c("glm", "lm"))) stop("model = 'ne' only supports yreg fitted via glm()")
    if (is.character(data[, exposure])) data[, exposure] <- as.factor(data[, exposure])
    
    environment(est.ne) <- environment()
    if (!multimp) {
      # point estimates of causal effects
      est <- est.ne(data = data, indices = NULL, outReg = TRUE, full = full)
      effect.pe <- est$est
      n_effect <- length(effect.pe)
      out$reg.output <- est$reg.output
      # bootstrap results
      boots <- boot(data = data, statistic = est.ne, R = nboot, outReg = FALSE, full = full)
      # bootstrap CIs
      environment(boot.ci) <- environment()
      effect.ci <- boot.ci(boots = boots)
      effect.ci.low <- effect.ci[, 1]
      effect.ci.high <- effect.ci[, 2]
      # bootstrap p-values
      effect.pval <- sapply(1:n_effect, function(x) boot.pval(boots = boots$t[, x], pe = effect.pe[x]))
    } else {
      # imputed data sets
      data_imp <- complete(do.call(mice, args_mice), action = "all")
      m <- length(data_imp)
      # estimate causal effects for each imputed data set
      est_imp <- lapply(1:m, function(x) est.ne(data = data_imp[[x]], indices = NULL, outReg = TRUE, full = full))
      est_imp_df <- do.call(rbind, lapply(1:m, function(x) est_imp[[x]]$est))
      effect.pe <- colMeans(est_imp_df)
      n_effect <- length(effect.pe)
      out$reg.output <- lapply(1:m, function(x) est_imp[[x]]$reg.output)
      boot.step <- function(data = NULL, indices = NULL) {
        data_boot <- data[indices, ]
        args_mice$data <- data_boot
        data_imp <- complete(do.call(mice, args_mice), action = "all")
        curVal <- get("counter", envir = env)
        assign("counter", curVal + 1, envir = env)
        setTxtProgressBar(get("progbar", envir = env), curVal + 1)
        return(colMeans(do.call(rbind, lapply(1:m, function(x)
          est.ne(data = data_imp[[x]], outReg = FALSE, full = full)))))
      }
      environment(boot.step) <- environment()
      # bootstrap results
      boots <- boot(data = data, statistic = boot.step, R = nboot)
      # bootstrap CIs
      environment(boot.ci) <- environment()
      effect.ci <- boot.ci(boots = boots)
      effect.ci.low <- effect.ci[, 1]
      effect.ci.high <- effect.ci[, 2]
      # bootstrap p-values
      effect.pval <- sapply(1:n_effect, function(x) boot.pval(boots = boots$t[, x], pe = effect.pe[x]))
    }
    
    if (family_yreg$family %in% c("gaussian", "inverse.gaussian", "Gamma", "quasi")) {
      # standard errors by bootstrapping
      effect.se <- sapply(1:n_effect, function(x) sd(boots$t[, x]))
      # effect names
      if (full) {
        if (EMint) effect_name <- c("cde", "pnde", "tnde", "pnie", "tnie", "te", 
                                    "intref", "intmed", "cde(prop)", "intref(prop)", "intmed(prop)", "pnie(prop)",
                                    "pm", "int", "pe")
        if (!EMint) effect_name <- c("cde", "pnde", "tnde", "pnie", "tnie", "te", "pm")
      }
      if (!full) effect_name <- c("cde", "pnde", "tnde", "pnie", "tnie", "te")
    } else {
      # transform standard errors of effects in log scale
      effect.se <- sapply(1:n_effect, function(x) ifelse(x <= 6, sd(exp(boots$t[, x])), sd(boots$t[, x])))
      # transform effects in log ratio scale into effects in ratio scale
      effect.pe[1:6] <- exp(effect.pe[1:6])
      effect.ci.low[1:6] <- exp(effect.ci.low[1:6])
      effect.ci.high[1:6] <- exp(effect.ci.high[1:6])
      # effect names
      if (full) {
        if (EMint) effect_name <- c("Rcde", "Rpnde", "Rtnde", "Rpnie", "Rtnie", "Rte", 
                                    "ERcde", "ERintref", "ERintmed", "ERpnie",
                                    "ERcde(prop)", "ERintref(prop)", "ERintmed(prop)", "ERpnie(prop)",
                                    "pm", "int", "pe")
        if (!EMint) effect_name <- c("Rcde", "Rpnde", "Rtnde", "Rpnie", "Rtnie", "Rte", "pm")
      }
      if (!full) effect_name <- c("Rcde", "Rpnde", "Rtnde", "Rpnie", "Rtnie", "Rte")
    }
    
    names(effect.pe) <- names(effect.se) <- names(effect.ci.low) <- names(effect.ci.high) <-
      names(effect.pval) <- effect_name
    out$effect.pe <- effect.pe
    out$effect.se <- effect.se
    out$effect.ci.low <- effect.ci.low
    out$effect.ci.high <- effect.ci.high
    out$effect.pval <- effect.pval
  }
  
  return(out)
}


boot.ci <- function(boots) {
  if (boot.ci.type == "per") {
    effect.ci.low <- sapply(1:n_effect, function(x) quantile(x = boots$t[, x], probs = 0.025, na.rm = TRUE))
    effect.ci.high <- sapply(1:n_effect, function(x) quantile(x = boots$t[, x], probs = 0.975, na.rm = TRUE))
    out <- cbind(effect.ci.low, effect.ci.high)
  } else if (boot.ci.type == "bca") {
    boot.ci.bca <- function(x) {
      xbar <- mean(x, na.rm = TRUE)
      z0 <- qnorm(length(which(x < xbar))/nboot)
      U <- x - xbar
      a <- sum(U^3, na.rm = TRUE)/(6*(sum(U^2, na.rm = TRUE))^(3/2))
      alpha1 <- pnorm(z0 + (z0 + qnorm(0.025))/(1 - a*(z0 + qnorm(0.025))))
      alpha2 <- pnorm(z0 + (z0 + qnorm(0.975))/(1 - a*(z0 + qnorm(0.975))))
      quantile(x, c(alpha1, alpha2)) 
    }
    out <- do.call(rbind, lapply(1:n_effect, function(x) boot.ci.bca(x = boots$t[, x])))
  }
  return(out)
}


boot.pval <- function(boots, pe){
  boots_noNA <- boots[which(!is.na(boots))]
  if (pe == 0) out <- 1
  if (pe != 0) out <- 2 * min(sum(boots_noNA > 0), sum(boots_noNA < 0)) / length(boots_noNA)
  return(out)
}




est.gformula <- function(data = NULL, indices = NULL, outReg = FALSE, full = TRUE) {
  if (is.null(indices)) indices <- 1:n
  # resample data
  data <- data[indices, ]
  
  # for case control study
  # method 1: weight subjects with y=1 by yprevalence/p(y=1) and weight subjects with y=0 by (1-yprevalence)/p(y=0)
  # method 2: fit yreg with all data and fit other regs on data among controls
  # use method 1 when yprevalence is provided
  # when yprevalence is not provided but the outcome is rare, use method 2
  if (casecontrol && !is.null(yprevalence)) {
    # method 1 for a case control design
    prob1 <- mean(data[, outcome] == y_case, na.rm = TRUE)
    w4casecon <- as.vector(ifelse(data[, outcome] == y_case, yprevalence / prob1, (1 - yprevalence) / (1 - prob1)))
    # weights for yreg
    if (!is.null(weights_yreg)) weights_yreg <- weights_yreg[indices] * w4casecon
    if (is.null(weights_yreg)) weights_yreg <- w4casecon
    # update yreg
    call_yreg$weights <- weights_yreg
    call_yreg$data <- data
    if (outReg && (inherits(yreg, "rcreg") | inherits(yreg, "simexreg"))) call_yreg$variance <- TRUE
    yreg <- eval.parent(call_yreg)
    # update mreg
    for (p in 1:length(mediator)) {
      if (!is.null(weights_mreg[[p]])) weights_mreg[[p]] <- weights_mreg[[p]][indices] * w4casecon
      if (is.null(weights_mreg[[p]])) weights_mreg[[p]] <- w4casecon
      call_mreg[[p]]$weights <- weights_mreg[[p]]
      call_mreg[[p]]$data <- data
      if (outReg && (inherits(mreg[[p]], "rcreg") | inherits(mreg[[p]], "simexreg"))) call_mreg[[p]]$variance <- TRUE
      mreg[[p]] <- eval.parent(call_mreg[[p]])
    }
    if (length(postc) != 0) {
      # update postcreg
      for (p in 1:length(postc)) {
        if (!is.null(weights_postcreg[[p]])) weights_postcreg[[p]] <- weights_postcreg[[p]][indices] * w4casecon
        if (is.null(weights_postcreg[[p]])) weights_postcreg[[p]] <- w4casecon
        call_postcreg[[p]]$weights <- weights_postcreg[[p]]
        call_postcreg[[p]]$data <- data
        if (outReg && (inherits(postcreg[[p]], "rcreg") | inherits(postcreg[[p]], "simexreg"))) call_postcreg[[p]]$variance <- TRUE
        postcreg[[p]] <- eval.parent(call_postcreg[[p]])
      }
    }
    rm(prob1, w4casecon)
  } else if (casecontrol && yrare) {
    # method 2 for a case control design
    # data from controls
    control_indices <- which(data[, outcome] == y_control)
    # update yreg
    call_yreg$weights <- weights_yreg[indices]
    call_yreg$data <- data
    if (outReg && (inherits(yreg, "rcreg") | inherits(yreg, "simexreg"))) call_yreg$variance <- TRUE
    yreg <- eval.parent(call_yreg)
    # update mreg
    for (p in 1:length(mediator)) {
      call_mreg[[p]]$weights <- weights_mreg[[p]][indices][control_indices]
      call_mreg[[p]]$data <- data[control_indices, ]
      if (outReg && (inherits(mreg[[p]], "rcreg") | inherits(mreg[[p]], "simexreg"))) call_mreg[[p]]$variance <- TRUE
      mreg[[p]] <- eval.parent(call_mreg[[p]])
    }
    # update postcreg
    if (length(postc) != 0) {
      for (p in 1:length(postc)) {
        # update postcreg[[p]]
        call_postcreg[[p]]$weights <- weights_postcreg[[p]][indices][control_indices]
        call_postcreg[[p]]$data <- data[control_indices, ]
        if (outReg && (inherits(postcreg[[p]], "rcreg") | inherits(postcreg[[p]], "simexreg"))) call_postcreg[[p]]$variance <- TRUE
        postcreg[[p]] <- eval.parent(call_postcreg[[p]])
      }
    }
    rm(control_indices)
  } else {
    # not a case control design
    # update yreg
    call_yreg$weights <- weights_yreg[indices]
    call_yreg$data <- data
    if (outReg && (inherits(yreg, "rcreg") | inherits(yreg, "simexreg"))) call_yreg$variance <- TRUE
    yreg <- eval.parent(call_yreg)
    # update mreg
    for (p in 1:length(mediator)) {
      call_mreg[[p]]$weights <- weights_mreg[[p]][indices]
      call_mreg[[p]]$data <- data
      if (outReg && (inherits(mreg[[p]], "rcreg") | inherits(mreg[[p]], "simexreg"))) call_mreg[[p]]$variance <- TRUE
      mreg[[p]] <- eval.parent(call_mreg[[p]])
    }
    # update postcreg
    if (length(postc) != 0) {
      for (p in 1:length(postc)) {
        call_postcreg[[p]]$weights <- weights_postcreg[[p]][indices]
        call_postcreg[[p]]$data <- data
        if (outReg && (inherits(postcreg[[p]], "rcreg") | inherits(postcreg[[p]], "simexreg"))) call_postcreg[[p]]$variance <- TRUE
        postcreg[[p]] <- eval.parent(call_postcreg[[p]])
      }
    }
  }
  
  # output list
  out <- list()
  if (outReg) {
    out$reg.output$yreg <- yreg
    out$reg.output$mreg <- mreg
    if (length(postc) != 0) out$reg.output$postcreg <- postcreg
  }
  
  # the index of the reference level for a categorical outcome
  if ((is_glm_yreg && (family_yreg$family %in% c("binomial", "quasibinomial", "multinom") |
                       startsWith(family_yreg$family, "Ordered Categorical"))) |
      is_multinom_yreg | is_polr_yreg) {
    y_lev <- levels(droplevels(as.factor(data[, outcome])))
    yval_index <- switch((yval %in% y_lev) + 1, "1" = NULL, "2" = which(y_lev == yval))
    rm(y_lev)
  }
  
  # simulate A
  if (is.factor(data[, exposure])) {
    a_sim <- factor(c(rep(a, n)), levels = a_lev)
    astar_sim <- factor(c(rep(astar, n)), levels = a_lev)
  } else {
    a_sim <- c(rep(a, n))
    astar_sim <- c(rep(astar, n))
  }
  
  # simulate C
  basec_sim <- data[, basec]
  
  # design matrices for simulating postc[p]
  postcdesign_a <- data.frame(a_sim, basec_sim)
  postcdesign_astar <- data.frame(astar_sim, basec_sim)
  colnames(postcdesign_a) <- colnames(postcdesign_astar) <- c(exposure, basec)
  postc_a <- postc_astar <- data.frame(matrix(nrow = n, ncol = length(postc)))
  colnames(postc_a) <- colnames(postc_astar) <- postc
  
  if (length(postc) != 0) {
    # simulating postc[p]
    for (p in 1:length(postc)) {
      # predict postc[p]
      type <- ifelse(is_multinom_postcreg[p] | is_polr_postcreg[p], "probs", "response")
      postcpred_a <- predict(postcreg[[p]], newdata = postcdesign_a, type = type)
      postcpred_astar <- predict(postcreg[[p]], newdata = postcdesign_astar, type = type)
      full_index <- which(rowSums(is.na(postcdesign_a))==0)
      n_full <- length(full_index)
      # categorical L
      if ((is_glm_postcreg[p] && ((family_postcreg[[p]]$family %in% c("binomial", "multinom")) |
                                  startsWith(family_postcreg[[p]]$family, "Ordered Categorical")))|
          is_multinom_postcreg[p] | is_polr_postcreg[p]) {
        l_lev <- levels(droplevels(as.factor(data[, postc[p]])))
        prob_a <- as.matrix(postcpred_a)
        prob_astar <- as.matrix(postcpred_astar)
        if (dim(prob_a)[2] == 1) {
          mid_a <- l_lev[rbinom(n_full, size = 1, prob = prob_a[full_index, 1]) + 1]
          mid_astar <- l_lev[rbinom(n_full, size = 1, prob = prob_astar[full_index, 1]) + 1]
        } else {
          mid_a <- l_lev[apply(prob_a[full_index,], 1, FUN = function(x) apply(t(rmultinom(1, 1, prob = x)), 1, which.max))]
          mid_astar <- l_lev[apply(prob_astar[full_index,], 1, FUN = function(x) apply(t(rmultinom(1, 1, prob = x)), 1, which.max))]
        }
        
        if (is.numeric(data[, postc[p]])) {
          mid_a <- as.numeric(l_lev[mid_a])
          mid_astar <- as.numeric(l_lev[mid_astar])
        } 
        
        rm(prob_a, prob_astar, l_lev)
        # linear L
      } else if ((is_lm_postcreg[p] | is_glm_postcreg[p]) && family_postcreg[[p]]$family == "gaussian") {
        error <- rnorm(n_full, mean = 0, sd = sigma(postcreg[[p]]))
        mid_a <- postcpred_a[full_index] + error
        mid_astar <- postcpred_astar[full_index] + error
        rm(error)
        # gamma L
      } else if (is_glm_postcreg[p] && family_postcreg[[p]]$family == "Gamma") {
        shape_postcreg <- MASS::gamma.shape(postcreg[[p]])$alpha
        mid_a <- rgamma(n_full, shape = shape_postcreg, scale = postcpred_a[full_index]/shape_postcreg)
        mid_astar <- rgamma(n_full, shape = shape_postcreg, scale = postcpred_astar[full_index]/shape_postcreg)
        rm(shape_postcreg)
        # inverse gaussian L
      } else if (is_glm_postcreg[p] && family_postcreg[[p]]$family == "inverse.gaussian") {
        lambda <- 1/summary(postcreg[[p]])$dispersion
        mid_a <- SuppDists::rinvGauss(n_full, nu = postcpred_a[full_index], lambda = lambda)
        mid_astar <- SuppDists::rinvGauss(n_full, nu = postcpred_astar[full_index], lambda = lambda)
        rm(lambda)
        # poisson L
      } else if (is_glm_postcreg[p] && family_postcreg[[p]]$family == "poisson") {
        mid_a <- rpois(n_full, lambda = postcpred_a[full_index])
        mid_astar <- rpois(n_full, lambda = postcpred_astar[full_index])
        # quasipoisson L
      } else if (is_glm_postcreg[p] && family_postcreg[[p]]$family == "quasipoisson") {
        dispersion <- summary(postcreg[[p]])$dispersion
        if (dispersion > 1) {
          mid_a <- sapply(1:n_full, function(i) predint::rqpois(1, lambda = postcpred_a[full_index][i], phi = dispersion)[,1])
          mid_astar <- sapply(1:n_full, function(i) predint::rqpois(1, lambda = postcpred_astar[full_index][i], phi = dispersion)[,1])
        } else {
          mid_a <- postcpred_a[full_index]
          mid_astar <- postcpred_astar[full_index]
        }
        rm(dispersion)
        # negative binomial L
      } else if (is_glm_postcreg[p] && startsWith(family_postcreg[[p]]$family, "Negative Binomial")) {
        theta <- summary(postcreg[[p]])$theta
        mid_a <- MASS::rnegbin(n_full, mu = postcpred_a[full_index], theta = theta)
        mid_astar <- MASS::rnegbin(n_full, mu = postcpred_astar[full_index], theta = theta)
        rm(theta)
      } else stop(paste0("Unsupported postcreg[[", p, "]]"))
      
      postc_a[full_index, p] <- mid_a
      postc_astar[full_index, p] <- mid_astar
      
      if (is.factor(data[, postc[p]])) {
        l_lev <- levels(droplevels(as.factor(data[, postc[p]])))
        postc_a[, p] <- factor(postc_a[, p], levels = l_lev)
        postc_astar[, p] <- factor(postc_astar[, p], levels = l_lev)
      }
    }
    
    rm(postcdesign_a, postcdesign_astar, type, postcpred_a, postcpred_astar, mid_a, mid_astar, full_index, n_full)
  }
  
  # design matrices for simulating mediator[p]
  mdesign_a <- data.frame(a_sim, basec_sim, postc_a)
  mdesign_astar <- data.frame(astar_sim, basec_sim, postc_astar)
  colnames(mdesign_a) <- colnames(mdesign_astar) <- c(exposure, basec, postc)
  m_a <- m_astar <- data.frame(matrix(nrow = n, ncol = length(mediator)))
  colnames(m_a) <- colnames(m_astar) <- mediator
  
  # simulating mediator[p]
  for (p in 1:length(mediator)) {
    # predict mediator[p]
    type <- ifelse(is_multinom_mreg[p] | is_polr_mreg[p], "probs", "response")
    mpred_a <- predict(mreg[[p]], newdata = mdesign_a, type = type)
    mpred_astar <- predict(mreg[[p]], newdata = mdesign_astar, type = type)
    full_index <- which(rowSums(is.na(mdesign_a))==0)
    n_full <- length(full_index)
    # categorical M
    if ((is_glm_mreg[p] && ((family_mreg[[p]]$family %in% c("binomial", "multinom")) |
                            startsWith(family_mreg[[p]]$family, "Ordered Categorical")))|
        is_multinom_mreg[p] | is_polr_mreg[p]) {
      m_lev <- levels(droplevels(as.factor(data[, mediator[p]])))
      prob_a <- as.matrix(mpred_a)
      prob_astar <- as.matrix(mpred_astar)
      if (dim(prob_a)[2] == 1) {
        # simulate mediator[p] for exposure=a
        mid_a <- m_lev[rbinom(n_full, size = 1, prob = prob_a[full_index, 1]) + 1]
        # simulate mediator[p] for exposure=astar
        mid_astar <- m_lev[rbinom(n_full, size = 1, prob = prob_astar[full_index, 1]) + 1]
      } else {
        mid_a <- m_lev[apply(prob_a[full_index,], 1, FUN = function(x) apply(t(rmultinom(1, 1, prob = x)), 1, which.max))]
        mid_astar <- m_lev[apply(prob_astar[full_index,], 1, FUN = function(x) apply(t(rmultinom(1, 1, prob = x)), 1, which.max))]
      }
      
      if (is.numeric(data[, mediator[p]])) {
        mid_a <- as.numeric(mid_a)
        mid_astar <- as.numeric(mid_astar)
      }
      
      rm(prob_a, prob_astar, m_lev)
      # linear M
    } else if ((is_lm_mreg[p] | is_glm_mreg[p]) && family_mreg[[p]]$family == "gaussian") {
      error <- rnorm(n_full, mean = 0, sd = sigma(mreg[[p]]))
      mid_a <- mpred_a[full_index] + error
      mid_astar <- mpred_astar[full_index] + error
      rm(error)
      # gamma M
    } else if (is_glm_mreg[p] && family_mreg[[p]]$family == "Gamma") {
      shape_mreg <- MASS::gamma.shape(mreg[[p]])$alpha
      mid_a <- rgamma(n_full, shape = shape_mreg, scale = mpred_a[full_index]/shape_mreg)
      mid_astar <- rgamma(n_full, shape = shape_mreg, scale = mpred_astar[full_index]/shape_mreg)
      rm(shape_mreg)
      # inverse gaussian M
    } else if (is_glm_mreg[p] && family_mreg[[p]]$family == "inverse.gaussian") {
      lambda <- 1/summary(mreg[[p]])$dispersion
      mid_a <- SuppDists::rinvGauss(n_full, nu = mpred_a[full_index], lambda = lambda)
      mid_astar <- SuppDists::rinvGauss(n_full, nu = mpred_astar[full_index], lambda = lambda)
      rm(lambda)
      # poisson M
    } else if (is_glm_mreg[p] && family_mreg[[p]]$family == "poisson") {
      mid_a <- rpois(n_full, lambda = mpred_a[full_index])
      mid_astar <- rpois(n_full, lambda = mpred_astar[full_index])
      # quasipoisson M
    } else if (is_glm_mreg[p] && family_mreg[[p]]$family == "quasipoisson") {
      dispersion <- summary(mreg[[p]])$dispersion
      if (dispersion > 1) {
        mid_a <- sapply(1:n_full, function(i) predint::rqpois(1, lambda = mpred_a[full_index][i], phi = dispersion)[,1])
        mid_astar <- sapply(1:n_full, function(i) predint::rqpois(1, lambda = mpred_astar[full_index][i], phi = dispersion)[,1])
      } else {
        mid_a <- mpred_a[full_index]
        mid_astar <- mpred_astar[full_index]
      }
      rm(dispersion)
      # negative binomial M
    } else if ( is_glm_mreg[p] && startsWith(family_mreg[[p]]$family, "Negative Binomial")) {
      theta <- summary(mreg[[p]])$theta
      mid_a <- MASS::rnegbin(n_full, mu = mpred_a[full_index], theta = theta)
      mid_astar <- MASS::rnegbin(n_full, mu = mpred_astar[full_index], theta = theta)
      rm(theta)
    } else stop(paste0("Unsupported mreg[[", p, "]]"))
    
    # randomly shuffle values of simulated mediator[p] if postc is not empty
    if (length(postc) != 0) {
      m_a[full_index, p] <- sample(mid_a, replace = FALSE)
      m_astar[full_index, p] <- sample(mid_astar, replace = FALSE)
    } else {
      m_a[full_index, p] <- mid_a
      m_astar[full_index, p] <- mid_astar
    }
    
    if (is.factor(data[, mediator[p]])) {
      m_lev <- levels(droplevels(as.factor(data[, mediator[p]])))
      m_a[, p] <- factor(m_a[, p], levels = m_lev)
      m_astar[, p] <- factor(m_astar[, p], levels = m_lev)
    }
    
  }
  rm(mdesign_a, mdesign_astar, type, mpred_a, mpred_astar, mid_a, mid_astar, full_index, n_full)
  
  # simulate mstar for cde
  mstar_sim <- do.call(cbind, lapply(1:length(mediator), function(x)
    if (is.factor(data[, mediator[x]])) {
      data.frame(factor(rep(mval[[x]], n), levels = levels(data[, mediator[x]])))
    } else data.frame(rep(mval[[x]], n))))
  
  # design matrices for outcome simulation
  ydesign0m <- data.frame(astar_sim, mstar_sim, basec_sim, postc_astar)
  ydesign1m <- data.frame(a_sim, mstar_sim, basec_sim, postc_a)
  ydesign00 <- data.frame(astar_sim, m_astar, basec_sim, postc_astar)
  ydesign01 <- data.frame(astar_sim, m_a, basec_sim, postc_astar)
  ydesign10 <- data.frame(a_sim, m_astar, basec_sim, postc_a)
  ydesign11 <- data.frame(a_sim, m_a, basec_sim, postc_a)
  rm(a_sim, astar_sim, m_a, m_astar, mstar_sim, basec_sim, postc_a, postc_astar)
  colnames(ydesign0m) <- colnames(ydesign1m) <- colnames(ydesign00) <- colnames(ydesign01) <-
    colnames(ydesign10) <- colnames(ydesign11) <- c(exposure, mediator, basec, postc)
  
  # predict Y
  type <- ifelse(is_coxph_yreg, "risk", ifelse(is_multinom_yreg | is_polr_yreg, "probs", "response"))
  EY0m_pred <- as.matrix(predict(yreg, newdata =  ydesign0m, type = type))
  EY1m_pred <- as.matrix(predict(yreg, newdata =  ydesign1m, type = type))
  EY00_pred <- as.matrix(predict(yreg, newdata =  ydesign00, type = type))
  EY01_pred <- as.matrix(predict(yreg, newdata =  ydesign01, type = type))
  EY10_pred <- as.matrix(predict(yreg, newdata =  ydesign10, type = type))
  EY11_pred <- as.matrix(predict(yreg, newdata =  ydesign11, type = type))
  rm(type, ydesign0m, ydesign1m, ydesign00, ydesign01, ydesign10, ydesign11)
  
  # weights of yreg
  weightsEY <- as.vector(model.frame(yreg)$'(weights)')
  if (is.null(weightsEY)) weightsEY <- rep(1, n)
  
  # categorical Y
  if ((is_glm_yreg && ((family_yreg$family %in% c("binomial", "quasibinomial", "multinom")) |
                       startsWith(family_yreg$family, "Ordered Categorical")))|
      is_multinom_yreg | is_polr_yreg) {
    if (!is.null(yval_index)) {
      if (dim(EY0m_pred)[2] == 1) {
        EY0m <- weighted_mean(cbind(1 - EY0m_pred, EY0m_pred)[, yval_index], w = weightsEY)
        EY1m <- weighted_mean(cbind(1 - EY1m_pred, EY1m_pred)[, yval_index], w = weightsEY)
        EY00 <- weighted_mean(cbind(1 - EY00_pred, EY00_pred)[, yval_index], w = weightsEY)
        EY01 <- weighted_mean(cbind(1 - EY01_pred, EY01_pred)[, yval_index], w = weightsEY)
        EY10 <- weighted_mean(cbind(1 - EY10_pred, EY10_pred)[, yval_index], w = weightsEY)
        EY11 <- weighted_mean(cbind(1 - EY11_pred, EY11_pred)[, yval_index], w = weightsEY)
      } else {
        EY0m <- weighted_mean(EY0m_pred[, yval_index], w = weightsEY)
        EY1m <- weighted_mean(EY1m_pred[, yval_index], w = weightsEY)
        EY00 <- weighted_mean(EY00_pred[, yval_index], w = weightsEY)
        EY01 <- weighted_mean(EY01_pred[, yval_index], w = weightsEY)
        EY10 <- weighted_mean(EY10_pred[, yval_index], w = weightsEY)
        EY11 <- weighted_mean(EY11_pred[, yval_index], w = weightsEY)
      }
    } else EY0m <- EY1m <- EY00 <- EY01 <- EY10 <- EY11 <- 0
  } else {
    # non-categorical Y
    EY0m <- weighted_mean(EY0m_pred, w = weightsEY)
    EY1m <- weighted_mean(EY1m_pred, w = weightsEY)
    EY00 <- weighted_mean(EY00_pred, w = weightsEY)
    EY01 <- weighted_mean(EY01_pred, w = weightsEY)
    EY10 <- weighted_mean(EY10_pred, w = weightsEY)
    EY11 <- weighted_mean(EY11_pred, w = weightsEY)
  }
  rm(weightsEY, EY0m_pred, EY1m_pred, EY00_pred, EY01_pred, EY10_pred, EY11_pred)
  
  # output causal effects in additive scale for continuous Y
  if ((is_lm_yreg | is_glm_yreg) &&
      (family_yreg$family %in% c("gaussian", "inverse.gaussian", "Gamma", "quasi"))) {
    cde <- EY1m - EY0m
    pnde <- EY10 - EY00
    tnde <- EY11 - EY01
    pnie <- EY01 - EY00
    tnie <- EY11 - EY10
    te <- tnie + pnde
    if (full) {
      pm <- tnie / te
      if (EMint) {
        intref <- pnde - cde
        intmed <- tnie - pnie
        cde_prop <- cde/te
        intref_prop <- intref/te
        intmed_prop <- intmed/te
        pnie_prop <- pnie/te
        int <- (intref + intmed)/te
        pe <- (intref + intmed + pnie)/te
        est <- c(cde, pnde, tnde, pnie, tnie, te, intref, intmed, 
                 cde_prop, intref_prop, intmed_prop, pnie_prop, pm, int, pe)
      } else est <- c(cde, pnde, tnde, pnie, tnie, te, pm)
    } else est <- c(cde, pnde, tnde, pnie, tnie, te)
  } else {
    # output causal effects in ratio scale for non-continuous Y
    
    ## output effects on the odds ratio scale for logistic regressions
    if (is_glm_yreg && family_yreg$family %in% c("binomial", "quasibinomial") &&
        family_yreg$link == "logit") {
      logRRcde <- log(EY1m/(1-EY1m)) - log(EY0m/(1-EY0m))
      logRRpnde <- log(EY10/(1-EY10)) - log(EY00/(1-EY00))
      logRRtnde <- log(EY11/(1-EY11)) - log(EY01/(1-EY01))
      logRRpnie <- log(EY01/(1-EY01)) - log(EY00/(1-EY00))
      logRRtnie <- log(EY11/(1-EY11)) - log(EY10/(1-EY10))
      ## otherwise on the risk ratio scale
    } else {
      logRRcde <- log(EY1m) - log(EY0m)
      logRRpnde <- log(EY10) - log(EY00)
      logRRtnde <- log(EY11) - log(EY01)
      logRRpnie <- log(EY01) - log(EY00)
      logRRtnie <- log(EY11) - log(EY10)
    }
    
    logRRte <- logRRtnie + logRRpnde
    if (full) {
      pm <- (exp(logRRpnde) * (exp(logRRtnie) - 1)) / (exp(logRRte) - 1)
      if (EMint) {
        ERRcde <- (EY1m-EY0m)/EY00
        ERRintref <- exp(logRRpnde) - 1 - ERRcde
        ERRintmed <- exp(logRRtnie) * exp(logRRpnde) - exp(logRRpnde) - exp(logRRpnie) + 1
        ERRpnie <- exp(logRRpnie) - 1
        ERRte <- exp(logRRte) - 1
        ERRcde_prop <- ERRcde/ERRte
        ERRintmed_prop <- ERRintmed/ERRte
        ERRintref_prop <- ERRintref/ERRte
        ERRpnie_prop <- ERRpnie/ERRte
        int <- (ERRintref + ERRintmed)/ERRte
        pe <- (ERRintref + ERRintmed + ERRpnie)/ERRte
        est <- c(logRRcde, logRRpnde, logRRtnde, logRRpnie, logRRtnie, logRRte, 
                 ERRcde, ERRintref, ERRintmed, ERRpnie,
                 ERRcde_prop, ERRintref_prop, ERRintmed_prop, ERRpnie_prop,
                 pm, int, pe)
      } else est <- c(logRRcde, logRRpnde, logRRtnde, logRRpnie, logRRtnie, logRRte, pm)
    } else est <- c(logRRcde, logRRpnde, logRRtnde, logRRpnie, logRRtnie, logRRte)
    
  } 
  
  # progress bar
  if (!multimp) {
    curVal <- get("counter", envir = env)
    assign("counter", curVal + 1, envir = env)
    setTxtProgressBar(get("progbar", envir = env), curVal + 1)
  }
  if (outReg) out$est <- est
  if (!outReg) out <- est
  return(out)
}




#power function
aim_3_power <- function(xy, xm, my, int, cov, sample_size, alpha, nsim, sd, srmu, lambda){
  
  cde_sig_results <- c()
  pie_sig_results <- c()
  pai_sig_results <- c()
  
  for (i in 1:nsim) {
    
    #simulate data
    aim_3_data <- zcta_msp %>%
      mutate(c = rnorm(sample_size, mean = 1, sd = 1),
                       sr = rnorm(sample_size, srmu, sqrt(srmu*(1-srmu)/sample_size)), 
             ps = ifelse(pop_2020==0,
                                   NA_integer_,
                                   (xm*sr+rpois(sample_size, lambda = lambda))/pop_2020),
      health = xy*sr + my*ps + cov*c + int*ps*sr + rnorm(sample_size, mean = 0, sd = sd)) %>%
      drop_na() %>%
      st_drop_geometry()
    
    #Counterfactual Mediation Analysis - gformula
    gformula <- cmest(data = aim_3_data, 
                          model = "gformula", 
                          outcome = "health", 
                          exposure = "sr",
                          mediator = "ps", 
                          basec = c("c"), 
                          EMint = TRUE,
                          mreg = list("linear"), 
                          yreg = "linear",
                          astar = 0, a = 1, 
                          mval = list(1),
                          estimation = "imputation", 
                          inference = "bootstrap", 
                          nboot = 2)
    
    #extract p-values
    cde_sig_results[i] <- as.numeric(gformula$effect.pval["cde"]) <= alpha
    pie_sig_results[i] <- as.numeric(gformula$effect.pval["pnie"]) <= alpha
    pai_sig_results[i] <- as.numeric(gformula$effect.pval["int"]) <= alpha
  }
  
  cde <- cde_sig_results %>% mean() 
  pie <- pie_sig_results %>% mean()
  pai <- pai_sig_results %>% mean()
  
  c(cde, pie, pai) %>%
    return()
}


aim_3_power(xy = .1, 
            xm = .0001, 
            my = .1,
            cov = .3,
            int = .1,
            srmu = .5,
            lambda = 1,
            sample_size = 50, 
            alpha = .05, 
            nsim = 100, 
            sd = .3)






