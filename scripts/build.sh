#!/bin/bash

export BUILD_MODE=2

# Bash framework
# -------------------------------
source "$(realpath "$(dirname "$0")")/bootstrap.sh"

# Let's build the app!
# -------------------------------
important "Building app..."
mvn verify clean package >/dev/null
