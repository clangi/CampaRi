context('keyfile handling')

test_that('Keyfile handling', {
  a <- keywords_from_keyfile(key_file_input='KEYWORD 123', return_table=TRUE, return_string_of_arguments=FALSE,
                             keyword_list_first=TRUE, key_file_is_keywords=TRUE)
  expect_true(!is.null(a))
  a <- keywords_from_keyfile(key_file_input='KEYWORD 123', return_table=FALSE, return_string_of_arguments=TRUE,
                             keyword_list_first=FALSE, key_file_is_keywords=TRUE)
  expect_true(!is.null(a))
})