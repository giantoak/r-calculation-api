#' A test function so we can examine how content looks when it's passed to opencpu
#' 
#' @param value: the value passed
#' @keywords show_contents
#' @export
#' @examples
#' show_contents()

show_contents <- function(target.region, comparison.region.set, event.date, input_data){
	library(jsonlite)
	input_data <- toJSON(input_data)
	input_data <- fromJSON(input_data)
	return(input_data)
}

