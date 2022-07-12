# SPDX-FileCopyrightText: 2022 Renaissance Computing Institute. All rights reserved.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-License-Identifier: LicenseRef-RENCI
# SPDX-License-Identifier: MIT

for s in `cat url.list`; do 
	if [ ${s:0:1} == h ] ; then 
		bash compute_geotiffs.sh $s
	else
		echo "not running $s"
	fi
done
