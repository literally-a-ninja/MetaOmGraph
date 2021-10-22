#!/usr/bin/env bash

# CLI Framework
# -------------------------------

function rel_pwd() {
	echo $(realpath "$(dirname "$0")$1")

	return 0;
}

export DIR_PROJECT="$(rel_pwd "/../src")"
export DIR_SCRIPT="$(rel_pwd "/.")"
export DIR_BOOT="$(rel_pwd "/bootstrap")"

source "${DIR_BOOT}/fun.sh"
source "${DIR_BOOT}/cli.sh"


# Pre-flight
# -------------------------------
JAVA_VERSION=$(java --version | head -n1 | awk -F '[^0-9]*' '$0=$2')
if [[ -z ${OJAVA_VERSION+z} && ${JAVA_VERSION} < 17 ]]; then
	error "Using $(link_man $(which java)) (JDK ${JAVA_VERSION}); build requires ${BOLD}JDK 17${NORMAL} or higher.";
	exit 1;
fi


# Dependencies
# -------------------------------
important "Checking build sys dependencies..."
source "${DIR_BOOT}/dependencies.sh"
ok "All dependencies are OK."

# Saved for later
# -------------------------------
# docker-gui() {
# 	xhost +si:localuser:root;
# 	docker run -it --rm -e DISPLAY \
# 		-v /tmp/.X11-unix:/tmp/.X11-unix \
# 		--tmpfs /dev/shm \
# 		"${@}"
# 
# 	return 0;
# }

# Libraries
# -------------------------------

source "${DIR_BOOT}/lib.sh"
important "Checking libraries..."

info "Processing l2fprod-common library..."
download_lib "http://www.java2s.com/Code/JarDownload/l2fprod/l2fprod-common-all.jar.zip" "l2fprod-common-all.jar" \
  && mvn install:install-file \
  -Dfile="$DIR_PROJECT/lib/l2fprod-common-all.jar" -Dpackaging=jar -DgeneratePom=true \
  -DgroupId='com.l2fprod' -DartifactId='l2fprod-common-all' -Dversion='0.1'

info "Processing CustomBrowserLauncher library..."
mvn install:install-file \
  -Dfile="$DIR_PROJECT/lib/CustomBrowserLauncher.jar" -Dpackaging=jar -DgeneratePom=true \
  -DgroupId='edu.iastate.metnet' -DartifactId='custombrowserlauncher' -Dversion='0.0.1'

ok "All libraries are OK."

# Done.
# -------------------------------

# Make sure our pwd is in project root.
cd "$DIR_PROJECT/.."


ok "Build system is READY, building in three seconds."
sleep 3
