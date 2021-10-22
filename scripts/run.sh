#!/bin/bash

# Bash framework
source "$(realpath "$(dirname "$0")")/bootstrap.sh"

docker build .
docker-gui .
