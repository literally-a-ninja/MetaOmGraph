#!/usr/bin/env bash

# CLI Framework
# -------------------------------

function rel_pwd() {
	echo $(realpath "$(dirname "$0")$1")

	return 0;
}

# Build system modes:
# 1. Minimum mode (helpers + dependencies ONLY)
# 2. Building mode (everything)
if [[ ! -v BUILD_MODE ]]; then
	export BUILD_MODE=1
fi
export MODE_MIN=1
export MODE_FULL=2

export DIR_PROJECT="$(rel_pwd "/../src")"
export DIR_ROOT="$(rel_pwd "/../")"
export DIR_SCRIPT="$(rel_pwd "/.")"
export DIR_BOOT="$(rel_pwd "/bootstrap")"

source "${DIR_BOOT}/fun.sh"
source "${DIR_BOOT}/cli.sh"


# Pre-flight
# -------------------------------
if [[ $BUILD_MODE > $MODE_MIN ]]; then
	JAVA_VERSION=$(java --version | head -n1 | awk -F '[^0-9]*' '$0=$2')
	if [[ -z ${OJAVA_VERSION+z} && ${JAVA_VERSION} < 17 ]]; then
		error "Using $(link_man $(which java)) (JDK ${JAVA_VERSION}); build requires ${BOLD}JDK 17${NORMAL} or higher.";
		exit 1;
	fi
fi


# Dependencies
# -------------------------------
important "Checking build sys dependencies..."
source "${DIR_BOOT}/dependencies.sh"
ok "All dependencies are OK."


# Done.
# -------------------------------

# Make sure our pwd is in project root.
cd "$DIR_SCRIPT/.."


if [[ $BUILD_MODE > $MODE_MIN ]]; then
	source "${DIR_BOOT}/libraries.sh"
	ok "Build system is READY, building in three seconds."
	sleep 3
fi
