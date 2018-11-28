#!/bin/bash

ulimit -s unlimited
ulimit -n 4096
ulimit -i 63040
ulimit -q 819200

./pinc FBinstability_SI_enhance.ini
