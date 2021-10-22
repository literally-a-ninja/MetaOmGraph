#!/usr/bin/env bash

apt() {
	if [[ -z "$1" ]]; then
		return 0;
	fi
	if [[ -f /usr/bin/apt ]]; then
		[[ -w "/var/lib/apt/lists/lock" ]] && ex="apt-get install" || ex="sudo apt-get install"
		echo "> $ex $1"
		info "Detected APT package manager! Proceeding with dependency installation..."

		$ex "$1"

		ok "Installed all dependenices! Continuing..."
		return 0
	else
		error "Failed to locate required dependencies ($1). Aborting build."
		return 1
	fi
}

deps=$(read -r -a array <<< "maven unzip wget")
locs=$(read -r -a array <<< "/usr/bin/mvn /usr/bin/unzip /usr/bin/wget")

missing=()
for (( i = 0; i < ${#deps[@]}; i++ )); do
	if [[ ! -f "${locs[$i]}" ]]; then
		missing+=("${deps[$i]}")
	fi
done

apt $missing
