#!/usr/bin/env bash

# Assume bash framework is already loaded.
# source "$(realpath "$(dirname "$0")")/bootstrap.sh"
# source "${DIR_BOOT}/libraries.sh"

name="metaomgraph4-jvm-$version"
cp "$DIR_ROOT/target/metaomgraph4-jar-with-dependencies.jar" "$DIR_ROOT/dist/metaomgraph4.jar"
cp "$DIR_ROOT/target/launcher.jar" "$DIR_ROOT/dist/launcher.jar"

zip "../$name.zip" "metaomgraph4.jar" "launcher.jar"

ok "Zip archive exported! ($(realpath ../$name.zip))"
