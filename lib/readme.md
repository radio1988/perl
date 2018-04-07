# the perl libs used by perl scripts
1. follow instructions in ../set_perl_lib.bash_profile to set perl lib path
2. alternatively, create a softlink of the specific pm file (e.g. READ_SAM.pm) in the folder where you run your perl scripts to make them work
  - to check which pm file to use, see `use XXX` in the perl code
