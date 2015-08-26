#' A wrapper for the diffindiff function. This is needed because when passing arguements from python they are converted to a list. This 
#' This will convert them back
#' 
#' @param target.region: The region of interest, matching a backpage domain
#' @param comparison.region.set: A set of regions (e.g. c('nova','abilene')) to compare to
#' @param event.date: a YYYY-MM-DD date string for the actual event date
#' @param input_data: The dataframe to use for diff-in-diff - needs to have 
#' @keywords diffindiff_wrapper
#' @export
#' @examples
#' diffindiff_wrapper()

diffindiff_wrapper <- function(target.region, comparison.region.set, event.date, input_data, standard.errors="naive"){
	library(jsonlite)
	input_data <- fromJSON(input_data)
	comparison.region.set <- fromJSON(comparison.region.set)
	result <- diffindiff_data(target.region = target.region, comparison.region.set = comparison.region.set, event.date = event.date, input_data = input_data, standard.errors = standard.errors)
	return(result)
}
