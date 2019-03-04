#!/bin/sh

DESTINATION="$1"

for n in {2,3,4}; do 
  cd "n=${n}"
  for file in plots/*.eps; do
    PDF_FILE="${DESTINATION}/$(basename ${file} .eps).pdf"
    ps2pdf -dEPSCrop "${file}" "../${PDF_FILE}"
  done
  cd ..
done
