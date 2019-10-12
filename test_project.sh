
# Use the 'testthat' package to run all tests in "tests/".
# Output is printed to screen and a separate log file in "logs/tests/"

testLogFile="logs/tests/project-test-log-$(date +"%Y-%m-%d_%H:%M:%S")"
Rscript -e "n; ProjectTemplate::test.project()" | tee $testLogFile
