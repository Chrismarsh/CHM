#
# Canadian Hydrological Model - The Canadian Hydrological Model (CHM) is a novel
# modular unstructured mesh based approach for hydrological modelling
# Copyright (C) 2018 Christopher Marsh
#
# This file is part of Canadian Hydrological Model.
#
# Canadian Hydrological Model is free software: you can redistribute it and/or
# modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Canadian Hydrological Model is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Canadian Hydrological Model.  If not, see
# <http://www.gnu.org/licenses/>.
#
# This script retrives the SVN version number of the tree on disk and displays
# a warning if it is unsuitable for an operational version (ie: composite version).
# cf http://subversion.tigris.org/faq.html#version-value-in-source
# Mathias Bavay - SLF - 02/2008

if [ $# -eq 1 ]; then
	FLAG="$1_"
else
	FLAG=""
fi

VERSION=$(svnversion -n ./)

NON_NUM=$(echo "${VERSION}" | grep "[^0-9]")
DATE=$(date +"%Y%m%d")

VER=$(printf "${FLAG}${DATE}.%s" "${VERSION}")
printf "***** Compiling Version ${VER}\n" > "/dev/stderr"

printf "${VER}"

