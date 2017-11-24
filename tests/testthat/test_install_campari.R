context('install_campari')

test_that('Test campari original installation from inst directory', {
  expect_error(install_campari(), NA)
  expect_error(install_campari(install_ncminer = T), NA)
  expect_error(install_campari(install_threads = T), NA)
  expect_error(install_campari(install_threads = T, install_ncminer = T), NA)
  expect_error(install_campari(install_threads = T, install_mpi = T), NA)
  expect_error(install_campari(install_threads = T, install_mpi = T), NA)
  expect_error(install_campari(install_ncminer = T, install_threads = T, install_mpi = T))
  expect_error(install_campari(installation_location = 'sadsdasdsadsadsadsadsadsafdfsdf', install_ncminer = T, install_threads = T, install_mpi = T))
})