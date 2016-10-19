#!/bin/bash

valgrind --leak-check=full --suppressions=pinc.supp --show-reachable=yes "$@"
