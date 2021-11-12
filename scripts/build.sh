#!/usr/bin/env bash

export BUILD_MODE=2

# Bash framework
# -------------------------------
source "$(realpath "$(dirname "$0")")/bootstrap.sh"
source "${DIR_BOOT}/libraries.sh"

## Libraries
## ===================== ##
cd "$DIR_ROOT"
source "$DIR_SCRIPT/bootstrap_libraries.sh"

ok "Build system is READY, building in three seconds."
sleep 3

## Build
## ===================== ##
important "Building app..."
mvn verify package
