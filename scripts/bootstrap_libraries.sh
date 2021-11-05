#!/usr/bin/env bash

# Bash framework
# -------------------------------
source "$(realpath "$(dirname "$0")")/bootstrap.sh"
source "${DIR_BOOT}/libraries.sh"


important "Checking libraries..."

pip3 install -r "$DIR_SCRIPT/build/mvnp/requirements" &>/dev/null

important "Processing l2fprod-common library..." \
[[ ! -f "$DIR_PROJECT/lib/l2fprod-common-all.jar" ]] \
  && download_libz "http://www.java2s.com/Code/JarDownload/l2fprod/l2fprod-common-all.jar.zip" "l2fprod-common-all.jar"
mvn install:install-file -Dfile="$DIR_PROJECT/lib/l2fprod-common-all.jar" -DgroupId='com.l2fprod' -DartifactId='l2fprod-common-all' -Dversion='0.1' -Dpackaging=jar -DgeneratePom=true 2>&1

important "Processing CustomBrowserLauncher library..."

mvn install:install-file -Dfile="$DIR_PROJECT/lib/CustomBrowserLauncher.jar" -DgroupId='edu.iastate.metnet' -DartifactId='custombrowserlauncher' -Dversion='0.0.1' -Dpackaging=jar -DgeneratePom=true 2>&1

important "Processing Hierarchial Clustering library..."
[[ ! -f "$DIR_PROJECT/lib/hierarchical-clustering-1.2.0.jar" ]] \
  && download_lib "https://github.com/lbehnke/hierarchical-clustering-java/releases/download/v1.2.0/hierarchical-clustering-1.2.0.jar" "hierarchical-clustering-1.2.0.jar"

mvn install:install-file -Dfile="$DIR_PROJECT/lib/hierarchical-clustering-1.2.0.jar" -DgroupId='com.apporiented' -DartifactId='hierarchical-clustering' -Dversion='1.2.0' -Dpackaging=jar -DgeneratePom=true 2>&1 
