# Running Test Case Simulations
Each of the test case directories contains the following files to support
running test case simulations:

* `params.in` - Simulation parameters file
* `<test_case>.run` - Job submission script to launch the `<test_case>` test
  case

The files above need to be modified before they can be used with and to launch
the simulation.  Each are templates that have the following placeholders that
need to be replaced:

* `PATH_TO_SCRATCH` - Location of the simulation's outputs.  For simulations
  running on Notre Dame systems, this must be beneath `/scratch365`
* `PATH_TO_CODE` - Location of the NTLP repository's working copy

The following commands will modify both `params.in` and `<test_case>.run`.  They
assume they are run from the `test_cases/` directory.  Change `TEST_CASE` as
appropriate.

```shell
TEST_CASE=channel
USER=agrace4
HOME=~
sed -i -e "s#PATH_TO_SCRATCH#/scratch365/${USER}#" ${TEST_CASE}/params.in
sed -i -e "s#PATH_TO_SCRATCH#/scratch365/${USER}#" \
       -e "s#PATH_TO_CODE#${HOME}/NTLP#" ${TEST_CASE}/${TEST_CASE}.run
```
