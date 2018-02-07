context('install_campari')

test_that('Test campari original installation from inst directory', {
  silent <- T
  # expect_error(install_campari(install_ncminer = T), NA)
  # expect_error(install_campari(install_threads = T, silent_built = F), NA)
  # expect_error(install_campari(install_threads = T, install_ncminer = T, silent_built = F), NA)
  # expect_error(install_campari(install_threads = T, install_mpi = T, silent_built = F), NA)
  # expect_error(install_campari(install_threads = T, install_mpi = T, silent_built = F), NA)
  expect_error(install_campari(install_ncminer = T, install_threads = T, install_mpi = T, silent_built = F))
  expect_error(install_campari(installation_location = 'sadsdasdsadsadsadsadsadsafdfsdf', install_ncminer = T, install_threads = T, install_mpi = T, silent_built = silent))
  
  if(!dir.exists('to_delete_just_in_a_moment')) dir.create('to_delete_just_in_a_moment')
  debugonce(install_campari)
  expect_error(install_campari(installation_location = 'to_delete_just_in_a_moment', install_ncminer = T, install_threads = F, silent_built = silent), NA)
  if(dir.exists('to_delete_just_in_a_moment')) unlink('to_delete_just_in_a_moment', recursive = T)
})