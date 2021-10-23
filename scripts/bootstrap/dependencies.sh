#!/usr/bin/env bash

if [[ $BUILD_MODE < $MODE_FULL ]]; then
docker-gui() {
	xhost +si:localuser:root;
	docker run -it --rm -e DISPLAY \
		-v /tmp/.X11-unix:/tmp/.X11-unix \
		--tmpfs /dev/shm \
		"${@}"
 
	return 0;
}
fi

apt() {
	if [[ -z "$1" ]]; then
		return 0;
	fi
	if [[ -f /usr/bin/apt ]]; then
		[[ -w "/var/lib/apt/lists/lock" ]] && ex="apt-get" || ex="sudo apt-get"
		echo "> $ex $1"
		info "Detected APT package manager! Proceeding with dependency installation..."

		$ex "install" "-y" $1

		ok "Installed all dependenices! Continuing..."
		return 0
	else
		error "Failed to locate required dependencies ($1). Aborting build."
		return 1
	fi
}

case "$BUILD_MODE" in
	"$MODE_MIN")
		deps=("docker")
		locs=("/usr/bin/docker")
	;;
	"$MODE_FULL")
		deps=("maven" "unzip" "wget" "file")
		locs=("/usr/bin/mvn" "/usr/bin/unzip" "/usr/bin/wget" "/usr/bin/file")
	;;
esac

missing=()
for (( i = 0; i < ${#deps[@]}; i++ )); do
	if [[ ! -f "${locs[$i]}" ]]; then
		missing+=("${deps[$i]}")
	fi
done

apt $missing
