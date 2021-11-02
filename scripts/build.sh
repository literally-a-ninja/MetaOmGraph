#!/bin/bash

export BUILD_MODE=2

# Bash framework
# -------------------------------
source "$(realpath "$(dirname "$0")")/bootstrap.sh"

# App libraries
# -------------------------------
source "$DIR_SCRIPT/bootstrap_libraries.sh"

# Let's build the app!
# -------------------------------
important "Building app..."
mvn verify clean package >/dev/null
