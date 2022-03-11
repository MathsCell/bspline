Sys.unsetenv("R_TESTS")
if (requireNamespace("RUnit", quietly=TRUE) && requireNamespace("bspline", quietly=TRUE)) {
  testSuite <- RUnit::defineTestSuite(
    name = "bspline unit tests",
    dirs = system.file("unitTests", package = "bspline"),
    testFuncRegexp = "^[Tt]est.+"
  )
  tests <- RUnit::runTestSuite(testSuite)

  RUnit::printTextProtocol(tests)

  if (RUnit::getErrors(tests)$nFail > 0) stop("RUnit test failure")
  if (RUnit::getErrors(tests)$nErr > 0) stop("Errors in RUnit tests")
}
