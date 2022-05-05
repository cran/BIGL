library(testthat)
library(BIGL)

#if (Sys.getenv("JENKINS_URL") != "")
#  options(testthat.default_check_reporter = "junit", testthat.output_file = "results.xml")

test_check("BIGL")
