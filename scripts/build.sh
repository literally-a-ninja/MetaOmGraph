#!/usr/bin/env bash

export BUILD_MODE=2

## Bash framework
## ===================== ##
source "$(realpath "$(dirname "$0")")/bootstrap.sh"
source "${DIR_BOOT}/libraries.sh"

ok "Build system is READY, building in three seconds."
sleep 3

## Build
## ===================== ##
important "Building app..."
mvn package

[[ -z "$CI_COMMIT_REF_SLUG" ]] \
  && version='latest' \
  || version="$CI_COMMIT_REF_SLUG"

# Make sure dist directory exists.
if [[ ! -d "$DIR_ROOT/dist" ]]; then
  mkdir "$DIR_ROOT/dist"
fi


# 1. Simple export jar + zip
name="metaomgraph4-jvm-$version"
cp "$DIR_ROOT/target/metaomgraph4-jar-with-dependencies.jar" "$DIR_ROOT/dist/metaomgraph4.jar"
cp "$DIR_ROOT/target/launcher.jar" "$DIR_ROOT/dist/launcher.jar"

cd "$DIR_ROOT/dist/"
zip "$DIR_ROOT/dist/$name.zip" "metaomgraph4.jar" "launcher.jar"

# 1. Export as zip archive
source "$DIR_ROOT/build/dpkg/target.sh"

# 2. Export as deb (linux)
source "$DIR_ROOT/build/zip/target.sh"

# 3. Export as package (macos)
# TODO
