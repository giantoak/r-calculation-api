#' Find similar regions based on the coveriate in the given data set
#' 
#' @param data: The data set for the function to analyze.
#' @param treatment_col: Name of the column we want to analze (usually 'msaname')
#' @param treatment_selection: Which value in the treatment_col we are trying to find comparisons for.
#' @param covariate_cols: What we want to use in the comparision (rape, violent, population etc).
#' @keywords matchit_data
#' @export
#' @examples
#' matchit_data()

matchit_data<-function(data, treatment_col, treatment_selection, covariate_cols){
	library(MatchIt)
	target <- "target_var"
      	score.varname <- "score_var"
	num.matches <- 3

	data[target] <- data[[treatment_col]] == treatment_selection
       
      	covs_plus = do.call(paste, c(as.list(covariate_cols), sep=" + "))
      	matching.formula <- paste(target, "~", covs_plus)
      	match <- matchit(as.formula(matching.formula), data)
      
      	data[score.varname] <- match$distance
      	target.row <- data[data[[target]],]
      	data<-data[!data[target],]
      	data<-data[order(data[score.varname], decreasing=TRUE),]
      	results <- unique(data[[treatment_col]])[1:num.matches]
      	return(results)
}

