important "Checking libraries..."

info "Processing l2fprod-common library..."
download_lib "http://www.java2s.com/Code/JarDownload/l2fprod/l2fprod-common-all.jar.zip" "l2fprod-common-all.jar" \
  && mvn install:install-file \
  -Dfile="$DIR_PROJECT/lib/l2fprod-common-all.jar" -Dpackaging=jar -DgeneratePom=true \
  -DgroupId='com.l2fprod' -DartifactId='l2fprod-common-all' -Dversion='0.1' >/dev/null

info "Processing CustomBrowserLauncher library..."
if [[ -f "$DIR_PROJECT/lib/CustomBrowserLauncher.jar" ]]; then
	mvn install:install-file \
	  -Dfile="$DIR_PROJECT/lib/CustomBrowserLauncher.jar" -Dpackaging=jar -DgeneratePom=true \
	  -DgroupId='edu.iastate.metnet' -DartifactId='custombrowserlauncher' -Dversion='0.0.1' >/dev/null
	ok "Success. Library installed via Maven"
else
	ok "Skipping. Already installed at $(link_file $DIR_PROJECT'/lib/CustomBrowserLauncher.jar')"
fi

ok "All libraries are OK."