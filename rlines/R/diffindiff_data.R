#' Take region names and a comparison date and does diff-in-diff
#' 
#' @param target.region: The region of interest, matching a backpage domain
#' @param comparison.region.set: A set of regions (e.g. c('nova','abilene')) to compare to
#' @param event.date: a YYYY-MM-DD date string for the actual event date
#' @param logged: A boolean as to whether to perform a log transform on the data first
#' @param input_data: The dataframe to use for diff-in-diff - needs to have 
#' @param date.var: String name for the date variable in the input_data dataframe.
#' @param group.var:  String name for the group (i.e. city name) variable in the input_data dataframe
#' @param var.of.interest: String name for the variable of interest (e.g. ad counts) in the input_data dataframe
#' @keywords max
#' @export
#' @examples
#' diffindiff_data()


diffindiff_data<-function(target.region, comparison.region.set, event.date, logged=FALSE, normalize=FALSE, input_data, date.var="monthdate", group.var="region", var.of.interest="counts", standard.errors="permutation"){
  #We need the reshape2 library to use the melt function
  library(reshape2)
  library(jsonlite)
  #Convert data to data frame
  input_data <- as.data.frame(input_data)
  #We need counts to be numberic
  input_data$counts <- as.numeric(input_data$counts)

  #I believe this is called a character vector. Whatever it's call, this is the form comparison.region.set needs to be in
  comparison.region.set <- c(comparison.region.set$region)
  
  post.var<-"post" # Set the variable name for "post"


  #Gives us counts for target and comparisons by date.  
  data<-twolines_data(target.region=target.region, comparison.region.set=comparison.region.set, data=input_data, date.var=date.var, group.var=group.var, var.of.interest=var.of.interest)
  print(data)
  print(input_data)
  #We need to make monthdate a date, and not a string.
  data[date.var]<-as.Date(data[[date.var]], "%Y-%m-%d")
  
  #Ensure event.date is type date
  event.date <- as.Date(event.date, "%Y-%m-%d")

  #Mark everything after the event date
  data[post.var] = data[[date.var]] > event.date
  #Transform the data by region
  data<-melt(data, id=c(date.var,post.var), variable.name=group.var, value.name=var.of.interest)
  #To look at results in term of percentage if logged is true
  if (logged){
    data[var.of.interest]<-log(1+data[var.of.interest])
  }
  #Both of these comes out as NaN. Why take data$group when there is no such column?
  pre.target.avg<-mean(data[data$post==FALSE & data$region == "Target","counts"])
  pre.comparison.avg<-mean(data[data$post==FALSE & data$region == "Comparison","counts"])
  #What is normalize for?
  if (normalize){
    data$counts[data$group == "Comparison"] <- data$counts[data$group == "Comparison"] * pre.target.avg/pre.comparison.avg
  }
  #print(data)
  form<-paste(var.of.interest," ~ ", post.var, "*", group.var)
  print(form)
  if (standard.errors == "naive"){
    # This is the naive approach, where we use actual standard errors from the regression
    data <- within(data, region <- relevel(data$region, "Comparison")) # Set comparison as base group
    #Create a string for an equation
    
    model<-lm(formula=form, data=data)
    # The idea for this model is for the results to be in the order:
    #  1: (Intercept)           
    #  2: postTRUE            
    #  3: groupTarget         
    #  4: postTRUE:groupTarget
    model.results<-coef(summary(model))
    vcov.matrix<-vcov(model)
    browser()
    dd<-list(
      b=model.results[4,"Estimate"],
      se=model.results[4, "Std. Error"],
      t=model.results[4,"t value"]
    )
    dd$p<-2*pt(-abs(dd$t),df=model$df.residual-1)
    ###### Code to compute target and comparison differences for simple model
    #   target.change.vec<-c(0,1,0,1)
    #   # The target change is the sum of the 2nd and 4th variables (set target=True and change post from 1 to 0)
    #   b.target<-target.change.vec %*% model.results[,"Estimate"]
    #   se.target<-sqrt(target.change.vec %*% vcov.matrix %*% target.change.vec)
    #   target.change<-list(
    #     b=b.target[1,1],
    #     se=se.target[1,1],
    #     t=b.target[1,1]/se.target[1,1]
    #     )
    #   comparison.vec<-c(0,1,0,0)
    #   # The comparison group is the sum of the 1st and 4th variables
    #   b.comparison<-comparison.vec %*% model.results[,"Estimate"]
    #   se.comparison<-sqrt(comparison.vec %*% vcov.matrix %*% comparison.vec)
    #   comparison.change<-list(
    #     b=b.comparison[1,1],
    #     se=se.comparison[1,1],
    #     t=b.comparison[1,1]/se.comparison[1,1]
    #   )
    #   comparison.change$p<-2*pt(-abs(comparison.change$t),df=df-1)
    #   # Note: we have to compute p manually since non-unit vectors will have covariance terms in SE
    #   # and hence won"t be in the main results
  } else if(standard.errors == "permutation"){
    # Here we use permutation based standard errors
    
    # To do this, we're starting with the input data with observations
    # at the MSA-month level, and defining POST and Treatment variables
    model.data<-input_data
    model.data$monthdate <-as.Date(model.data$monthdate, "%Y-%m-%d")
    model.data$new.region<-model.data$region == target.region
    model.data$region<-model.data$new.region
    model.data$new.region<-NULL
    #data$monthdate<-as.Date(data[[date.var]], "%Y-%m-%d")
    model.data$post<-model.data$monthdate > event.date
    print(model.data)
    
    
    # Determine permutation list of dates
    date.list<-unique(model.data$monthdate) # find unique dates
    date.list<-date.list[2:(length(date.list)-1)]
    # Restrict unique dates to those that are not the very ends so that
    # difference in differences is well defined.
    max.dates<-40
    if (length(date.list) > max.dates){
      event.date.sample<-sample(date.list,max.dates, replace=TRUE)
    } else{
      event.date.sample<-date.list
    }
    
    # Take a sample of the event dates
    results<-data.frame(monthdate=event.date.sample)
    
    estimate_dd_model<-function(data, form){
      return(lm(formula=form, data=data))
    }
    results=c()
    for (ed in event.date.sample) {
      # We're going to only compute our regression for each event date once
      model.data$post <- model.data$monthdate > ed
      model.results<-estimate_dd_model(data=model.data, form=form)
      results.matrix<-coef(summary(model.results))
      results<-c(results,results.matrix[1,"Estimate"])
      
      # Needs to be 4 for the dd estimate
    } 
    result.df<-data.frame(date=event.date.sample, results=results)
    point.estimate<-result.df[result.df$date==event.date,'results']  
    number.of.results.larger <-sum(result.df$results > point.estimate)
    number.of.results.smaller <-sum(result.df$results < point.estimate)
    p.smaller <- (1+number.of.results.smaller)/dim(result.df)
    p.larger <- (1+number.of.results.larger)/dim(result.df)
    dd<-list(
      b=point.estimate,
      p=min(p.smaller,p.larger)
    )
    #model<-estimate_dd_model(data=model.data, form=form)
    #print(model)
  }


  
  # Clean up data to send back
  comparison<-data[data[group.var] == "Comparison",c(date.var,var.of.interest)]
  comparison[date.var] <- strftime(comparison[[date.var]],"%Y-%m-%d")
  target<-data[data[group.var] == "Target",c(date.var,var.of.interest)]
  target[date.var] <- strftime(target[[date.var]],"%Y-%m-%d")
  print(data)
  data<-reshape2::dcast(data=data, formula=paste(date.var,"~",group.var), value=eval(parse(text=var.of.interest)))
  return(list(data=data,
              comparison=comparison,
              target=target,
              #model=model, 
              diff_in_diff=dd
              #target_diff=target.change,
              #comparison_diff=comparison.change
              ))
}
