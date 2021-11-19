#!/usr/bin/env bash

# Assume bash framework is already loaded.
# source "$(realpath "$(dirname "$0")")/bootstrap.sh"
# source "${DIR_BOOT}/libraries.sh"

# Build archive and place in dir above.
dpkg-deb -v --build "dpkg/" "../metaomgraph4-all.deb"

ok "Debian package exported! ($(realpath ../metaomgraph4-all.deb))"
