#!/usr/bin/env bash

loadProperties () {

# Copyright(C) 2019 Jason Fleming
# SPDX-FileCopyrightText: 2022 Renaissance Computing Institute. All rights reserved.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-License-Identifier: LicenseRef-RENCI
# SPDX-License-Identifier: MIT

propertiesFile=$1

   while read -r keyValuePair ; do 
      colonPosition=`expr index "$keyValuePair" ":" `
	  #colonPosition=`echo "$keyValuePair" | sed -n "s/[":"].*//p" | wc -c`
      key=${keyValuePair:0:$colonPosition-1}
      value=${keyValuePair:$colonPosition}
      # remove leading whitespace characters from key
      key="${key#"${key%%[![:space:]]*}"}" 
      # remove trailing whitespace characters from key
      key="${key%"${key##*[![:space:]]}"}"
      # remove leading whitespace characters from value
      value="${value#"${value%%[![:space:]]*}"}"
      # remove trailing whitespace characters from value
      value="${value%"${value##*[![:space:]]}"}"
      properties["$key"]="$value" 
      #echo $key, $value
   done < $propertiesFile
}
