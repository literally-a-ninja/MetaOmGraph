#!/usr/bin/env bash

dnc() {
	FILE_INFO="$(file $1)"
	FILE_LOCATION="$(dirname $1)"
	FILE_NAME="$(basename $1)"
	case "$FILE_INFO" in
		'Zip')
			mv "$FILE_LOCATION/$FILE_NAME{,.zip}" &&
				unzip "$FILE_LOCATION/$FILE_NAME.zip" -od "$DIR_PROJECT/lib" &>/dev/null &&
				rm -f "$FILE_LOCATION/$FILE_NAME.zip" &>/dev/null;
		;;

		'gzip')
			mv "$FILE_LOCATION/$FILE_NAME{,.gz}" &&
				gzip -d "$FILE_LOCATION/$FILE_NAME.gz" &>/dev/null
		;;

		*)
			error "Could not decode unknown archive type at \"$(link_file $FILE_DEST)\"."
		;;
	esac

	return 0;
}

download_lib() {
	FILE_DEST="$DIR_PROJECT/lib/$2"

	if [[ -f "$FILE_DEST" ]]; then
		ok "Skipping. Already installed at $(link_file $FILE_DEST)"
		return 0;
	fi

	wget "$1" -O"$FILE_DEST"

	if [[ "$(file $FILE_DEST)" =~ 'archive data' ]]; then
		dnc "$FILE_DEST"
	fi

	return 0;
}


