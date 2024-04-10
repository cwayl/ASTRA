#!/bin/bash

# show an menu
dialog --title "Mode Select" \
       --backtitle "ASTRA" \
       --nocancel \
       --menu "Select TLE Input Mode" 12 45 25 \
       1 "Use NORAD ID (Space-Track.org)" \
       2 "Manually Enter New TLE" \
       3 "Use current TLE" \
       2>/home/astra/ASTRA/text.txt
