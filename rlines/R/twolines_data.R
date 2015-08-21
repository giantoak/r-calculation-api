#' Take region names and lookup series from a panel
#' 
#' @param target.region: The region of interest, matching a backpage domain
#' @param comparison.region.set: A set of regions (e.g. c('nova','abilene')) to compare to
#' @param data: The data frame to split up - needs columns 'region', 'MonthDate', and 'counts'
#' @keywords max
#' @export
#' @examples
#' twolines_data()

twolines_data<-function(target.region, comparison.region.set, data,date.var="monthdate", var.of.interest="counts",  group.var="region"){
  library(plyr)
  target.varname<-"Target"
  comparison.varname<-"Comparison"
  #Get only the target data
  output<-data[data[group.var] == target.region,c(date.var,var.of.interest)]
  #Rename the columsn
  names(output)<-c(date.var,target.varname)
  #Get the comparison data
  comparison<-data[data[[group.var]] %in% comparison.region.set,c(date.var,var.of.interest)]
  #Sum the comparison data
  comparison<-ddply(comparison, date.var, summarize, counts=sum(counts))
  #Rename the columns
  names(comparison)<-c(date.var,comparison.varname)
  #Merge the data
  output<-merge(output,comparison)
  return(output)
}
