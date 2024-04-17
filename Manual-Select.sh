#!/bin/bash

# show an menu
dialog --title "Mode Select" \
       --backtitle "ASTRA" \
       --nocancel \
       --menu "Select TLE Input Mode" 12 45 25 \
       1 "Using NORAD ID (Space-Track.org)" \
       2 "Manually Entering New TLE" \
       3 "Using current TLE" \
       2>/home/astra/ASTRA/text.txt
