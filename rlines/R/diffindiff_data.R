#' Take region names and a comparison date and does diff-in-diff
#' 
#' @param target.region: The region of interest, matching a backpage domain
#' @param comparison.region.set: A set of regions (e.g. c('nova','abilene')) to compare to
#' @param event.date: a YYYY-MM-DD date string for the actual event date
#' @param logged: A boolean as to whether to perform a log transform on the plotting.data first
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
  #Convert plotting.data to plotting.data frame
  input_data <- data.frame(input_data)
  
  #We need counts to be numeric, and for the date to be a Date datatype
  input_data$counts <- as.numeric(input_data$counts) 
  input_data[date.var] <-as.Date(input_data[[date.var]], "%Y-%m-%d")

  #I believe this is called a character vector. Whatever it's call, this is the form comparison.region.set needs to be in
  comparison.region.set <- c(comparison.region.set$region)
  
  post.var<-"post" # Set the variable name for "post"


  #Gives us counts for target and comparisons by date.
  # This data will be returned to potentially plot lines on a graph
  plotting.data<-twolines_data(target.region=target.region, comparison.region.set=comparison.region.set, data=input_data, date.var=date.var, group.var=group.var, var.of.interest=var.of.interest)
  #We need to make monthdate a date, and not a string.
  plotting.data[date.var]<-as.Date(plotting.data[[date.var]], "%Y-%m-%d")
  
  #Ensure event.date, the date of our event of interest is type Date
  event.date <- as.Date(event.date, "%Y-%m-%d")

  #Mark everything after the event date in the plotting data
  plotting.data[post.var] = plotting.data[[date.var]] > event.date
  #Transform the plotting.data by region
  plotting.data<-melt(plotting.data, id=c(date.var,post.var), variable.name=group.var, value.name=var.of.interest)
  
  #To look at results in term of percentage if logged is true
  if (logged){
    plotting.data[var.of.interest]<-log(1+plotting.data[var.of.interest])
    input_data$counts <- log(1+input_data$counts)
  }

  #Normalize will move the comparison and target lines so that they display
  # well on the same axis on a graph. This only applies to the plotting.data data frame and doesn't
  # affect analysis.
  if (normalize){
    #Both of these comes out as NaN. Why take plotting.data$group when there is no such column?
    pre.target.avg<-mean(plotting.data[plotting.data$post==FALSE & plotting.data$region == "Target","counts"])
    pre.comparison.avg<-mean(plotting.data[plotting.data$post==FALSE & plotting.data$region == "Comparison","counts"])
    plotting.data$counts[plotting.data$group == "Comparison"] <- plotting.data$counts[plotting.data$group == "Comparison"] * pre.target.avg/pre.comparison.avg
  }
  
  # Create an R formula object representing the model we will fit for the difference-in-differences estimate
  form<-paste(var.of.interest," ~ ", post.var, "*", group.var)

  if (standard.errors == "naive"){
    # This is the "naive" approach, where we use actual standard errors from the regression. 
    
    # The next lines replace the "region" variable with "Comparison" and "Target" values.
    # For sure there is a smoother way to do this.
    input_data$new <- "Comparison"
    input_data[input_data$region == target.region, c('new')]<-"Target"
    input_data$new<- as.factor(input_data$new)
    input_data[group.var]<-input_data$new
    input_data$new <- NULL
    input_data <- within(input_data, region <- relevel(input_data$region, "Comparison")) # Set comparison as base group
    
    # create the 'post' variable
    input_data$post <- input_data[[date.var]] > event.date

    #Fit our model to the data
    model<-lm(formula=form, data=input_data)
    
    # The idea for this model is for the results to be in the order:
    #  1: (Intercept)           
    #  2: postTRUE            
    #  3: groupTarget         
    #  4: postTRUE:groupTarget
    # Then, the 4th element is the difference in differences estimate.
    
    model.results<-coef(summary(model)) # Grab results from the model
    vcov.matrix<-vcov(model)
    dd<-list(
      b=model.results[4,"Estimate"],
      se=model.results[4, "Std. Error"],
      t=model.results[4,"t value"]
    )
    dd$p<-2*pt(-abs(dd$t),df=model$df.residual-1)
    # Put the results of the model into the 'dd' object for export.
    # Note:
    # b - estimate of the actual difference in differences effect - this is interpreted as the effect of an event on counts in the target region, assuming that otherwise the target region would have behaved like the comparison region.
    # se - the statistical uncertainty around the estiamte
    # t - the t-statistic... if you have to ask, go for a nice explanation on wikipedia!
    # p - the probability that an effect this large would occur if the REAL effect of the intervention were 0, just due to sampling chance
    ###### Code to compute target and comparison differences for simple model. This is no longer used.
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
    # Here we use permutation based standard errors, permuted over possible valid event dates.
    # The basic idea is to check all the different POTENTIAL events that could have occured and 
    # answer the question: how big an effect is this the actual event, compared to all other potential events.
    # The p-value that comes back represents how many potential events have bigger (in absolute value)
    # effects than the entered event.
    
    # To do this, we're starting with the input plotting.data with observations
    # at the MSA-month level, and defining POST and Treatment variables
    
    # Set up the model dataframe here
    model.data<-input_data
    model.data$monthdate <-as.Date(model.data$monthdate, "%Y-%m-%d")
    model.data$new.region<-model.data$region == target.region
    model.data$region<-model.data$new.region
    model.data$new.region<-NULL
    model.data$post<-model.data$monthdate > event.date
    
    # Determine permutation list of dates
    date.list<-unique(model.data$monthdate) # find unique dates
    date.list<-date.list[2:(length(date.list)-1)] 
    # We can't have events on the first or last day...
    # otherwise difference in differences is not well defined.
    
    max.dates<-40
    # If we have more than 40 potential dates, let's just take a sample. Otherwise, we will
    # use all the dates we have. The reason to take a sample is that fitting the model >100 times might be
    # computationally difficult, but will have diminishing usefulness.
    if (length(date.list) > max.dates){
      event.date.sample<-sample(date.list,max.dates, replace=TRUE)
    } else{
      event.date.sample<-date.list
    }
    
    estimate_dd_model<-function(data, form){
      # This is a local function which will actually estimate the model.
      # We will loop over this for each potential event date
      return(lm(formula=form, data=plotting.data))
    }
    results=c()
    for (ed in event.date.sample) {
      # We're going to compute our regression for each date in the sample
      
      model.data$post <- model.data$monthdate > ed
      model.results<-estimate_dd_model(data=model.data, form=form)
      results.matrix<-coef(summary(model.results))
      results<-c(results,results.matrix[4,"Estimate"])
      
      # Needs to be 4 for the dd estimate
    } 
    result.df<-data.frame(date=event.date.sample, results=results) # Put our results into a data.frame
    point.estimate<-result.df[result.df$date==event.date,'results'] # Our actual estimate of dd
    
    # Compute how many effects for different dates are bigger and smaller than our estimate.
    number.of.results.larger <-sum(result.df$results > point.estimate)
    number.of.results.smaller <-sum(result.df$results < point.estimate)
    
    # The p-values here are just the fraction of results which are smaller (bigger) than our result
    p.smaller <- (1+number.of.results.smaller)/dim(result.df)
    p.larger <- (1+number.of.results.larger)/dim(result.df)
    # Note: we'll send back the SMALLER p value of the two, since the estimate is either closer to being 
    # above or below normal, by chance.
    
    dd<-list(
      b=point.estimate,
      p=min(p.smaller,p.larger)
    )
    # Create the dd return object
  }


  
  # Clean up plotting.data to send back
  comparison<-plotting.data[plotting.data[group.var] == "Comparison",c(date.var,var.of.interest)]
  comparison[date.var] <- strftime(comparison[[date.var]],"%Y-%m-%d")
  target<-plotting.data[plotting.data[group.var] == "Target",c(date.var,var.of.interest)]
  target[date.var] <- strftime(target[[date.var]],"%Y-%m-%d")
  print(plotting.data)
  plotting.data<-reshape2::dcast(data=plotting.data, formula=paste(date.var,"~",group.var), value=eval(parse(text=var.of.interest)))
  # Reshape the plotting.data dataframe to return like:
  # monthdate, Target, Comparison
  # 2011-08-01, 1, 8
  # Instead of:
  # monthdate, region, counts
  # 2011-08-01, Target, 1
  # 2011-08-01, Comparison, 8
  return(list(data=plotting.data,
              comparison=comparison,
              target=target,
              diff_in_diff=dd
              ))
}
