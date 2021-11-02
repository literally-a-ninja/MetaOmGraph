#!/usr/bin/env bash

export BUILD_MODE=1

# Bash framework
# -------------------------------
source "$(realpath "$(dirname "$0")")/bootstrap.sh"

# Hotfixes
# -------------------------------
source "$DIR_SCRIPT/hotfixes/linux.sh"

# Build docker image
# -------------------------------
if [[ ! -z "$(docker image inspect iastate/mog-app:latest 2>&1 | grep "Error")" ]]; then
	important "Starting docker build...."
	docker build --tag iastate/mog-app:latest --file "$DIR_ROOT/Dockerfile" $DIR_ROOT
	ok "Finished build."
fi

# Run docker container
# -------------------------------
ok "Launching application!"
docker inspect localbuild-mog &>/dev/null
if [[ "$?" > 0 ]]; then
	docker container run -ti --name localbuild-mog \
		--net=host \
		--env="DISPLAY" \
		--volume="$HOME/.Xauthority:/root/.Xauthority:rw" \
 		--volume /tmp/.X11-unix:/tmp/.X11-unix \
 		iastate/mog-app
else
	docker container start localbuild-mog
fi
