#!/bin/bash

# Bash framework
source "$(realpath "$(dirname "$0")")/bootstrap.sh"


info "MVN: Generating build system objects..."
mvn validate initialize generate-sources process-sources generate-resources process-resources

info "MVN: Installing application..."
mvn clean package
