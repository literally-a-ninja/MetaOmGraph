#!/usr/bin/env bash

export BUILD_MODE=2

# Bash framework
# -------------------------------
source "$(realpath "$(dirname "$0")")/bootstrap.sh"

# Testing!
# -------------------------------
debug "---------------------------------------------------"
debug "                   BEGIN TESTING                   "
debug "---------------------------------------------------"

mvn test

debug "---------------------------------------------------"
debug "                    END TESTING                    "
debug "---------------------------------------------------"
