#!/bin/bash

source scripts/helpers.sh

# Initialization
# -------------------------------

if [[ ! -f /usr/bin/mvn ]]; then
  if [[ -f /usr/bin/apt ]]; then
    important "Maven build system is a requirement of this project."
    echo "> apt-get install maven"
    apt-get update
    apt-get install maven

    ok "Installed Maven build system. Re-run the script."
  else
    error "Failed to locate maven build system (/usr/bin/mvn). Aborting build."
  fi

  exit 1;
fi

# Dependencies
# -------------------------------

debug "Installing l2fprod-common library..."
download_libz "http://www.java2s.com/Code/JarDownload/l2fprod/l2fprod-common-all.jar.zip" "l2fprod-common-all.jar" \
  && mvn install:install-file \
  -Dfile=src/lib/l2fprod-common-all.jar -Dpackaging=jar -DgeneratePom=true \
  -DgroupId='com.l2fprod' -DartifactId='l2fprod-common-all' -Dversion='0.1'

debug "Installing CustomBrowserLauncher library..."
mvn install:install-file \
  -Dfile=src/lib/CustomBrowserLauncher.jar -Dpackaging=jar -DgeneratePom=true \
  -DgroupId='edu.iastate.metnet' -DartifactId='custombrowserlauncher' -Dversion='0.0.1'

ok "Dependencies installed."

# Building
# -------------------------------

debug "MVN: Generating build system objects..."
mvn validate initialize generate-sources process-sources generate-resources process-resources

debug "MVN: Installing application..."
mvn clean package
