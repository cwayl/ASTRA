#!/bin/bash

# show an inputbox
dialog --title "Current Stats" \
       --backtitle "ASTRA" \
       --infobox "Azimuth: $1\n
Elevation: $2\n
Range: $3 km\n
Range Rate: $4 km/s" 12 45
