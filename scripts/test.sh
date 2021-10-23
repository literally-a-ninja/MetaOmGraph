#!/bin/bash

export BUILD_MODE=2

# Bash framework
source "$(realpath "$(dirname "$0")")/bootstrap.sh"

debug "---------------------------------------------------"
debug "                   BEGIN TESTING                   "
debug "---------------------------------------------------"

mvn test

debug "---------------------------------------------------"
debug "                    END TESTING                    "
debug "---------------------------------------------------"
