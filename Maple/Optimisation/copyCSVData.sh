#!/bin/sh

DESTINATION="$1"

for n in {2,3,4}; do 
	[ -d "${DESTINATION}/n=${n}" ] || mkdir "${DESTINATION}/n=${n}"

	cp "n=${n}/data"/*.csv "${DESTINATION}/n=${n}"
done
