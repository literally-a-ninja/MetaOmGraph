#!/bin/bash

export BUILD_MODE=2

# Bash framework
source "$(realpath "$(dirname "$0")")/bootstrap.sh"

# Libraries
# -------------------------------
if [[ $BUILD_MODE > $MODE_MIN ]]; then
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
fi


# Let's build the app!
# -------------------------------
important "Building app..."
mvn verify clean package
