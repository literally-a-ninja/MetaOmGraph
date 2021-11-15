#!/usr/bin/env bash

export BUILD_MODE=1

## CLI Framework
## ===================== ##
source "$(realpath "$(dirname "$0")")/bootstrap.sh"

## Apply hotfixes
## ===================== ##
## source "$DIR_SCRIPT/hotfixes/linux.sh"

## Build app
## ===================== ##
## DOCKER_CONTAINER="localbuild-mog"
## DOCKER_IMG_TAG="iastate/mog-app:latest"
## if [[ ! -z "$(docker image inspect "$DOCKER_IMG_TAG" 2>&1 | grep "Error")" ]]; then
##   docker container rm -f "$DOCKER_CONTAINER" &>/dev/null
##   important "Starting docker build...."
##   docker build -t "$DOCKER_IMG_TAG" --compress $DIR_ROOT 
## 
##   # Exit on failure.
##   if [[ ! -z "$?" ]]; then
##     exit 1
##   fi
## 
##   ok "Finished build."
## fi


## bash ./build.sh

## Copy JAR
## ===================== ##
## docker cp "$DOCKER_CONTAINER":/opt/mog/mog.jar "$DIR_ROOT/metaomgraph-latest.jar"

cd "$DIR_ROOT"

BUILD_OUTPUT_SELECTOR="/target/metaomgraph4*.jar"

[[ -z "$CI_COMMIT_REF_SLUG" ]] \
  && version='latest' \
  || version="$CI_COMMIT_REF_SLUG"

name="metaomgraph4-jvm-$version"

cp "$DIR_ROOT/target/metaomgraph4*.jar" "$DIR_ROOT/$name.jar"
cp "$DIR_ROOT/target/metaomgraph4*.jar" "$DIR_ROOT/install/metaomgraph4.jar"
ok "Exported applet as $name.jar!"

tar czvf "$name.tar.gz" "$DIR_ROOT/install/*" 
ok "Exported install archive as $name.tar.gz."

zip -9 "$name.zip" "$DIR_ROOT/install/*" 
ok "Exported install archive as $name.zip."

## Launch app
## ===================== ##
## ok "Launching application!"
## if [[ ! -z "$(docker container inspect "$DOCKER_CONTAINER" 2>&1 | grep "Error")" ]]; then
##   docker container run -ti --name "$DOCKER_CONTAINER" \
##     --net=host \
##     --env="DISPLAY" \
##     --volume="$HOME:/$HOME:rw" \
##     --volume="/tmp/.X11-unix:/tmp/.X11-unix:rw" \
##     iastate/mog-app:latest
## else
##   docker container start localbuild-mog
## fi
