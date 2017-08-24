context('writing a pdb')


test_that('writing a pdb', {
  data_f <- matrix(rnorm(n = 1000), nrow = 100, ncol = 10)
  expect_that(write_pdb(trj = data_f, file_name = "acciaio_inox.pdb"), not(throws_error()))
  expect_that(write_pdb(trj = data_f, file_name = "acciaio_inox.pdb", seq_name = "seq.in"), not(throws_error()))
  expect_that(write_pdb(trj = data_f[,1:9], file_name = "acciaio_inox.pdb", method = "bio3d"), not(throws_error()))
  expect_that(write_pdb(trj = data_f[,1:9], file_name = "acciaio_inox.pdb", method = "automatic_unsafe"), not(throws_error()))
  expect_that(write_pdb(trj = data_f, file_name = "acciaio_inox.pdb", method = "automatic_safe", dim_check = T), not(throws_error()))
  data_f <- 'ATOM  1    CH3    ACE  A    1 30.11 46.91 45.81     1    0
ATOM  2      C    ACE  A    1 31.33 46.73 46.64     1    0
ATOM  3      O    ACE  A    1 32.07 45.75 46.56     1    0
ATOM  4     1H    ACE  A    1 29.28 46.28 46.20     1    0
ATOM  5     2H    ACE  A    1 29.84 47.98 45.98     1    0
ATOM  6     3H    ACE  A    1 30.35 46.85 44.73     1    0'
  library(data.table)
  data_f <- fread(data_f, data.table = F)
  expect_that(write_pdb(trj = data_f, file_name = "acciaio_inox.pdb", method = "automatic_safe", dim_check = F), not(throws_error()))
  if(file.exists('acciaio_inox.pdb')) file.remove("acciaio_inox.pdb")
  if(file.exists('seq.in')) file.remove("seq.in")
})
