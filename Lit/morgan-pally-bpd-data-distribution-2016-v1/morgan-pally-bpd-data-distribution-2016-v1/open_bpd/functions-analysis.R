# Primary author: Stephen L. Morgan

#### define function to create a data frame for a particular window of time ####

window.df.and.zoo <- function(df, outcome, end.date, start.date, end.pre.date) {

  ### Function to generate windowed zoo and data frame objects
  ###   for a full time period (start through end)
  ###   and for a pre-treatment period (start through pre)
  ### Inputs
  ###  df: data frame to be windowed
  ###  outcome: column number of a crime or arrest category (integer)
  ###  end.date: last date for full time period (date format)
  ###  start.date: first date for full time period (date format)
  ###  end.pre.date: last date for pre-treatment period (date format)
  ### Returns:
  ###  list containing 2 data frame objects and corresponding 2 zoo objects
  ### Example (where 172 is column for "homicide): 
  ###  window.df.and.zoo(crimes.type.week, 172, "2015-12-28", "2010-03-01", "2014-08-04")

  ## select variables from data frame and filter
  
  df.windowed.pre <- select(df, week.first, y = outcome, dark.before.12, tmax.f, 
                            precip.in, snow.in, season, school) %>%
                     dplyr::filter(week.first >= start.date,
                            week.first <= end.pre.date)
  df.windowed.all <- select(df, week.first, y = outcome, dark.before.12, tmax.f, 
                            precip.in, snow.in, season, school) %>%
                     dplyr::filter(week.first >= start.date,
                            week.first <= end.date)
  
  ## convert to windowed zoo objects
  
  zoo.windowed.pre <- zoo(df.windowed.pre, df.windowed.pre$week.first)
  zoo.windowed.all <- zoo(df.windowed.all, df.windowed.all$week.first)
  zoo.windowed.pre$week.first <- NULL
  zoo.windowed.all$week.first <- NULL
  
  ## return list
  
  list.df.and.zoo <- list(df.windowed.pre, df.windowed.all,
                          zoo.windowed.pre, zoo.windowed.all)
  return(list.df.and.zoo)
}

#### define function to estimate period effects as "interventions" ####

intervention.models <- function(dataset, i, end.date, start.date, end.pre.date) {

  ### Estimate five linear models from windowed dataset, return as a list
  ###   calls on window.df.zoo() function
  ### Inputs
  ###  dataset: data frame with data to be modeled
  ###  outcome: column number of a crime or arrest category (integer)
  ###  end.date: last date for full time period (date format)
  ###  start.date: first date for full time period (date format)
  ###  end.pre.date: last date for pre-treatment period (date format)
  ### Returns:
  ###  list containing five linear models and mean predicted value for 
  ###   the last 52 weeks of the pre-Ferguson period

  ## set up data, specifying outcome as i and setting ending date for window
  ##  use window.df.and.zoo function defined elsewhere
  
  w <- window.df.and.zoo(dataset, i, end.date, start.date, end.pre.date)
  df.windowed.pre <- w[[1]]
  df.windowed.all <- w[[2]]
  zoo.windowed.pre <- w[[3]]
  zoo.windowed.all <- w[[4]]
  
  df.windowed.pre.original <- df.windowed.pre
  
  ## create linear time variables
  
  df.windowed.pre$t <- 1:length(df.windowed.pre$y)
  df.windowed.all$t <- 1:length(df.windowed.all$y)
  
  ## create linear spline for temperature
  
  parameterize.spline <- function(x, c) ifelse (x > c, x - c, 0)
  tmax.f.knots <- c(0, 50, 60, 70, 80) 
  df.windowed.pre$tmax.f.spline <- outer(df.windowed.pre$tmax.f, 
                                         tmax.f.knots, parameterize.spline)
  df.windowed.all$tmax.f.spline <- outer(df.windowed.all$tmax.f, 
                                         tmax.f.knots, parameterize.spline)

  ## create spike and period type variables
  
  df.windowed.pre$period.ferguson <- as.numeric(df.windowed.pre$week.first >= as.Date("2014-08-11"))
  df.windowed.all$period.ferguson <- as.numeric(df.windowed.all$week.first >= as.Date("2014-08-11"))
  
  df.windowed.pre$period.gray <- as.numeric(df.windowed.pre$week.first >= as.Date("2015-04-20"))
  df.windowed.all$period.gray <- as.numeric(df.windowed.all$week.first >= as.Date("2015-04-20"))
  
  df.windowed.pre$period.davis <-
    as.numeric(df.windowed.pre$week.first >= as.Date("2015-07-13"))
  df.windowed.all$period.davis <- 
    as.numeric(df.windowed.all$week.first >= as.Date("2015-07-13"))
  
  df.windowed.pre$spike.gray <- as.numeric(df.windowed.pre$week.first == as.Date("2015-04-27"))
  df.windowed.all$spike.gray <- as.numeric(df.windowed.all$week.first == as.Date("2015-04-27"))
  
  ## specify model, estimate, and store selected output
  
  naive.int.model.formula <- as.formula(paste("y ~ t + 
                                              period.ferguson + 
                                              period.gray + spike.gray + period.davis"))
  covariate.model.formula <- as.formula(paste("y ~ t + tmax.f.spline + snow.in +
                                              precip.in + dark.before.12 + school"))
  full.int.model.formula <- as.formula(paste("y ~ t + 
                                             period.ferguson + 
                                             period.gray + spike.gray + period.davis +
                                             tmax.f.spline + snow.in +
                                             precip.in + dark.before.12 + school"))
  constrained.int.model.formula <- as.formula(paste("y.diff ~  
                                                    period.ferguson + 
                                                    period.gray + spike.gray + period.davis"))
  
  ls.naive.int <- lm(naive.int.model.formula, df.windowed.all)
  ls.pre <- lm(covariate.model.formula, df.windowed.pre)
  ls.all <- lm(covariate.model.formula, df.windowed.all)
  ls.full.int <- lm(full.int.model.formula, df.windowed.all)
  df.windowed.all$y.predicted <- predict(ls.pre, df.windowed.all,
                                         type = "response")
  df.windowed.all$y.diff <- df.windowed.all$y - df.windowed.all$y.predicted
  ls.constrained.int <- lm(constrained.int.model.formula, df.windowed.all)

  ## calculate mean predicted value for prior 52 weeks
  
  df.windowed.pre.original$y.predicted <- predict(ls.pre, df.windowed.pre, 
                                               type = "response")
  df.windowed.pre52 <- dplyr::filter(df.windowed.pre.original,
                              week.first >= as.Date("2013-07-29"),
                              week.first <= end.pre.date.for.window)
  mean.pre52.predicted <- mean(df.windowed.pre52$y.predicted)

  ## return list
  
  return(list(ls.pre, ls.all, ls.naive.int, ls.constrained.int, ls.full.int,
              mean.pre52.predicted))
}

#### define functions for least squares model estimation, followed by plot ####

fit.least.squares.pre.ferg <- function(dataset, i, end.date, start.date, end.pre.date) {

  ## set up data, specifying outcome as i and setting dates

  names <- c("Ferguson Protests Begin", "Gray Protests Begin", "Davis Appointed")
  doi <- c(as.Date("2014-08-11"), as.Date("2015-04-20"), as.Date("2015-07-13"))
  w <- window.df.and.zoo(dataset, i, end.date, start.date, end.pre.date)
  df.windowed.pre <- w[[1]]
  df.windowed.all <- w[[2]]
  zoo.windowed.pre <- w[[3]]
  zoo.windowed.all <- w[[4]]
  
  ## create linear time variables
  
  df.windowed.pre$t <- 1:length(df.windowed.pre$y)
  df.windowed.all$t <- 1:length(df.windowed.all$y)
  
  ## create linear spline for temperature
  
  parameterize.spline <- function(x, c) ifelse (x > c, x - c, 0)
  tmax.f.knots <- c(0, 50, 60, 70, 80) 
  df.windowed.pre$tmax.f.spline <- outer(df.windowed.pre$tmax.f, 
                                         tmax.f.knots, parameterize.spline)
  df.windowed.all$tmax.f.spline <- outer(df.windowed.all$tmax.f, 
                                         tmax.f.knots, parameterize.spline)
  
  ## specify, estimate, and store model
  
  model.formula <- as.formula(paste("y ~ t + tmax.f.spline + snow.in +
                                    precip.in + dark.before.12 + school"))
  fitted.model <- lm(model.formula, df.windowed.pre)

  ## generate dummy variable for post period, and assemble
  ##   predicted and observed values for the full windowed series
  
  df.windowed.all$post <- df.windowed.all$week.first > as.Date("2014-08-10")
  df.windowed.all$y.predicted <- predict(fitted.model, df.windowed.all)
  predicted.and.observed <- as.data.frame(cbind(df.windowed.all$post, 
                                                df.windowed.all$y.predicted, 
                                                df.windowed.all$y))
  colnames(predicted.and.observed) <- c("post", "predicted", "observed")
  predicted.and.observed$week.first <- as.Date(df.windowed.all$week.first)
  
  ##  plot using functions defined in this file
  
  plot.post <- plot.pre.ferg.model(predicted.and.observed, 
                                   variable.names[i], names, doi)
  plot(plot.post)
  plot.with.smooth <- plot.pre.ferg.model.with.smooth(predicted.and.observed, 
                                                      variable.names[i], names, doi)
  plot(plot.with.smooth)                               
  
  ## plot autocorrelation functions
  
  acf.y <- acf(df.windowed.pre$y, lag.max = 52, plot = FALSE)
  acf.resid <- acf(fitted.model$residuals, lag.max = 52, plot = FALSE)
  
  par(mfrow=c(2, 1))
  plot(acf.y, main = " ")
  title(main = "Autocorrelation Function for the Observed Outcome for Model (2)",
        font.main = 1)
  plot(acf.resid, main = " ")
  title(main = "Autocorrelation Function for the Residuals from Model (2)",
        font.main = 1)
  par(mfrow=c(1, 1))
  
}

#### define function for poisson model estimation, followed by plot ####

fit.poisson.pre.ferg <- function(dataset, i, end.date, start.date, end.pre.date) {
  
  ## set up data, specifying outcome as i and setting dates
  
  names <- c("Ferguson Protests Begin", "Gray Protests Begin", "Davis Appointed")
  doi <- c(as.Date("2014-08-11"), as.Date("2015-04-20"), as.Date("2015-07-13"))
  w <- window.df.and.zoo(dataset, i, end.date, start.date, end.pre.date)
  df.windowed.pre <- w[[1]]
  df.windowed.all <- w[[2]]
  zoo.windowed.pre <- w[[3]]
  zoo.windowed.all <- w[[4]]
  
  ## create linear time variables
  
  df.windowed.pre$t <- 1:length(df.windowed.pre$y)
  df.windowed.all$t <- 1:length(df.windowed.all$y)
  
  ## create linear spline for temperature
  
  parameterize.spline <- function(x, c) ifelse (x > c, x - c, 0)
  tmax.f.knots <- c(0, 50, 60, 70, 80) 
  df.windowed.pre$tmax.f.spline <- outer(df.windowed.pre$tmax.f, 
                                         tmax.f.knots, parameterize.spline)
  df.windowed.all$tmax.f.spline <- outer(df.windowed.all$tmax.f, 
                                         tmax.f.knots, parameterize.spline)
  
  ## specify, estimate, and store model
  
  model.formula <- as.formula(paste("y ~ t + tmax.f.spline + snow.in +
                                    precip.in + dark.before.12 + school"))
  fitted.model <- glm(model.formula, df.windowed.pre, family = poisson)
  
  ## display estimated model
  
  print(summary(fitted.model))
  
  ## generate dummy variable for post period, and assemble
  ##   predicted and observed values for the full windowed series
  
  df.windowed.all$post <- df.windowed.all$week.first > as.Date("2014-08-10")
  df.windowed.all$y.predicted <- predict(fitted.model, df.windowed.all,
                                         type = "response")
  predicted.and.observed <- as.data.frame(cbind(df.windowed.all$post, 
                                                df.windowed.all$y.predicted, 
                                                df.windowed.all$y))
  colnames(predicted.and.observed) <- c("post", "predicted", "observed")
  predicted.and.observed$week.first <- as.Date(df.windowed.all$week.first)
  
  ##  plot using functions defined in this file
  
  plot.post <- plot.pre.ferg.model(predicted.and.observed, 
                                   variable.names[i], names, doi)
  plot(plot.post)  
  plot.with.smooth <- plot.pre.ferg.model.with.smooth(predicted.and.observed, 
                                                      variable.names[i], names, doi)
  plot(plot.with.smooth)                               
}

#### define functions to plot predicted and observed data ####

plot.pre.ferg.model <- function(w, outcome.name, names, dates) {
  
  ## assign colors to pre and post periods of data, based on post dummy
  
  w$cplot <- as.character(gsub(0, "black", (w$post)) %>% gsub("1", "red",.))
  
  ## plot figure
  
  p <- ggplot(data= w, aes(x = week.first)) + 
    geom_point(aes(y = observed), color = "grey") +
    geom_line(aes(y = predicted), color = w$cplot) +
    theme_bw() +
    theme(plot.title = element_text(size=14, face="bold", vjust=2,
                                    family="Palatino"),
          axis.title.y = element_text(family="Palatino"),
          axis.title.x = element_text(family="Palatino"),
          axis.text.y = element_text(family="Palatino"),
          axis.text.x = element_text(family="Palatino")) +
    ylab("Count per week") + xlab("") + labs(title = outcome.name) +
    geom_vline(xintercept = as.numeric(dates), colour="black", 
               show_guide = TRUE) + 
    annotate("text", label = names, x = dates, y = 0, size = 4, 
             colour = "black", angle = 90, vjust = -0.4, hjust = 0,
             family="Palatino")
  return(p)
}

plot.pre.ferg.model.with.smooth <- function (w, outcome.name, names, dates) {
  
  ## calculate smoothed values for observed data
  
  ma.three <- rep(1/3, 3)
  w$observed.smoothed <- stats::filter(w$observed, ma.three, sides=2)
  
  ## assign colors to pre and post periods of data, based on post dummy
  
  w$cplot <- as.character(gsub(0, "black", (w$post)) %>% gsub("1", "red",.))
  
  ## plot figure
  
  p <- ggplot(data= w, aes(x = week.first)) + 
    geom_line(aes(y = observed.smoothed), color = "lightblue") +
    geom_line(aes(y = predicted), color = w$cplot) +
    theme_bw() +
    theme(plot.title = element_text(size=14, face="bold", vjust=2,
                                    family="Palatino"),
          axis.title.y = element_text(family="Palatino"),
          axis.title.x = element_text(family="Palatino"),
          axis.text.y = element_text(family="Palatino"),
          axis.text.x = element_text(family="Palatino")) +
    ylab("Count per week") + xlab("") + labs(title = outcome.name) +
    geom_vline(xintercept = as.numeric(dates), colour="black", 
               show_guide = TRUE) + 
    annotate("text", label = names, x = dates, y = 0, size = 4, 
             colour = "black", angle = 90, vjust = -0.4, hjust = 0,
             family="Palatino")
  return(p)
}



