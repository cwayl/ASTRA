#!/bin/bash

# show an inputbox
dialog --title "Station Coordinates" \
       --backtitle "ASTRA" \
       --nocancel \
       --inputbox "Enter the station Altitude (km)" 8 60 $1 2>/Users/cooperwayland/Desktop/text.txt
