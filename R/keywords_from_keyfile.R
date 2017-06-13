#' @title Transform a keyfile in a list of keywords
#' @description
#'      Transform an input keyfile in a list (or table) of keywords which can be modified and used again in \code{\link{run_campari}}. 
#'      If \code{return_table} is \code{TRUE} a data.frame of the keywords is returned (it should be used mainly in the case of multivalues keywords, 
#'      that would be truncated to the first value otherwise).
#' @param key_file_input An already formatted keyfile. It can contain comments (#) and blank lines. Please note that in the case of multivariable keywords 
#' the algorithm will keep only the first value. Consider using \code{return_table} option in such a case. This variable can be also set to a character string 
#' with correct formatting if the \code{key_file_is_keywords} is turned on (\code{TRUE}).
#' @param return_table Defaults to \code{FALSE}. It will return a formatted table instead of a named list.
#' @param keyword_list If provided, the function will append the readed keywords to this list and return it.
#' @param return_string_of_arguments This argument will provide the possibility to have a string with already formatted function arguments (can be directly copied in \code{\link{run_campari}}).
#' @param keyword_list_first It defines if the inserted keyword_list must be inserted afterwords or before (\code{TRUE}) the loaded keyfile keywords.
#' @param key_file_is_keywords If \code{TRUE} the function will treat the key_file_input as a text element with keyfile alike features (attention: for this no fault catch is active).
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

keywords_from_keyfile <- function(key_file_input, return_table=FALSE, keyword_list=NULL, return_string_of_arguments=FALSE,
                                  keyword_list_first=TRUE, key_file_is_keywords=FALSE){
  
  if(!is.character(key_file_input))
    stop('Inserted key_file_input is not a char string.')
  if(!file.exists(key_file_input) && !key_file_is_keywords)
    stop('Inserted key_file_input does not exists.')
  if(!is.logical(key_file_is_keywords))
    stop('key_file_is_keywords must be a logical.')
  if(key_file_is_keywords)
    warning("Using the option of directly insert the variables as a string (key_file_is_keywords) delivers to the user the duty to format them as a keyfile. No further fault detection on the format is performed.")
  
  if(!is.logical(return_table))
    stop('return_table must be a logical.')
  if(!is.logical(keyword_list_first))
    stop('keyword_list_first must be a logical.')
  if(!is.logical(return_string_of_arguments))
    stop('return_string_of_arguments must be a logical.')
  if(return_table && return_string_of_arguments)
    stop('It is not possible to return a string format AND a table for return. Please choose one.')
  
  if(!is.null(keyword_list)){
    if(keyword_list_first)
      cat('A predefined list of keywords have been inserted. The ones from the keyfile will be appended after those.\n')
    else
      cat('A predefined list of keywords have been inserted. The ones from the keyfile will be appended before those.\n')
    key_names <- names(keyword_list)
  }else{
    key_names<-NULL
  }
  
  if(!key_file_is_keywords) what_to_read <- paste0('sed "s/#.*//" ', key_file_input)
  else what_to_read <- key_file_input
  keywords_table <- fread(input = what_to_read, blank.lines.skip = TRUE, header = FALSE, fill = TRUE) # fill is for multivalues keywords
  colnames(keywords_table) <- c('keyword', paste0('value', 1:(ncol(keywords_table)-1)))
  if(keyword_list_first){
    keyword_list <- c(keyword_list, c(keywords_table[,2])[[1]])
    names(keyword_list) <- c(key_names, c(keywords_table[,1])[[1]])
  }else{
    keyword_list <- c(c(keywords_table[,2])[[1]], keyword_list)
    names(keyword_list) <- c(c(keywords_table[,1])[[1]], key_names)
  }
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
