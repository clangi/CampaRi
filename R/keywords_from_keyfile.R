
#' @title Transform a keyfile in a list of keywords
#' @description
#'      Transform an input keyfile in a list (or table) of keywords which can be modified and used again in \code{\link{run_campari}}. 
#'      If \code{return_table} is \code{TRUE} a data.frame of the keywords is returned (it should be used mainly in the case of multivalues keywords, 
#'      that would be truncated to the first value otherwise).
#' @param key_file_input An already formatted keyfile. It can contain comments (#) and blank lines. Please note that in the case of multivariable keywords 
#' the algorithm will keep only the first value. Consider using \code{return_table} option in such a case.
#' @param return_table Defaults to \code{FALSE}. It will return a formatted table instead of a named list.
#' @param keyword_list If provided, the function will append the readed keywords to this list and return it.
#' @param return_string_of_arguments This argument will provide the possibility to have a string with already formatted function arguments (can be directly copied in \code{\link{run_campari}}).
#' @return A named list or a table of keywords readed from the keyfile in input
#' @seealso
#' \code{\link{run_campari}}, \code{\link{sapphire_plot}}.
#' @examples
#' \dontrun{
#'  keywords_from_keyfile("keyfile.key")
#' }
#'
#' @importFrom data.table fread
#' 
#' @export keywords_from_keyfile
#' 

keywords_from_keyfile <- function(key_file_input, return_table=FALSE, keyword_list=NULL, return_string_of_arguments=FALSE){
  
  if(!is.character(key_file_input) || !file.exists(key_file_input))
    stop('Inserted key_file_input does not exist.')
  
  if(!is.logical(return_table))
    stop('return_table must be a logical.')
  if(!is.logical(return_string_of_arguments))
    stop('return_string_of_arguments must be a logical.')
  if(return_table && return_string_of_arguments)
    stop('It is not possible to return a string format AND a table for return. Please choose one.')
  
  if(!is.null(keyword_list)){
    cat('Keyfile manually inserted for reformatting. Attention: all the standard checks will be overriden.\n')
    key_names <- names(keyword_list)
  }else{
    key_names<-NULL
  }
  
  keywords_table <- fread(input = paste0('sed "s/#.*//" ', key_file_input), blank.lines.skip = TRUE, header = FALSE, fill = TRUE) # fill is for multivalues keywords
  colnames(keywords_table) <- c('keyword', paste0('value', 1:(ncol(keywords_table)-1)))
  keyword_list <- c(keyword_list, c(keywords_table[,2])[[1]])
  names(keyword_list) <- c(key_names, c(keywords_table[,1])[[1]])
  cat('Keyfile successfully loaded.', nrow(keywords_table), 'new keywords found.\n')
  
  if(return_table){
    if(!is.null(keyword_list))
      warning('The return_table option will not consider the keyword_list to be appended. It will return only the keyfile table.')
    return(keywords_table)
  }else if(return_string_of_arguments){
    keywords_string <- paste0(names(keyword_list), "=", as.vector(unlist(keyword_list[names(keyword_list)])), collapse = ", ")
    return(keywords_string)
  }else{
    return(keyword_list)
  }
}
