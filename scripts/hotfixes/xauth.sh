#!/usr/bin/env bash

if [[ -z "$(xhost | tail -n+2 | grep 'localuser:root')" ]]; then
	warn "Cannot connect to X11 display server!"
	echo "The root user is not permitted to connect to the X11 display server."
	echo "As a result, Docker is conversely not permitted to use X11 unless you authorise this hotfix."; echo

	ANS=$(read -p"Permit root user to connect to X11? [Y/n] ")
	if [[ "$ANS" != "n" || "$ANS" != "N" ]]; then
		sudo xhost SI:localuser:root
	fi
fi
