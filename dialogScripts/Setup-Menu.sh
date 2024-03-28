#!/bin/bash

# show an menu
dialog --title "Setup Menu" \
       --backtitle "ASTRA" \
       --nocancel \
       --menu "Select Setting" 12 45 25 \
       1 "Station Coordinates" \
       2 "Space-Track.com Credentials" \
       3 "" \
       4 " Back to Start Menu" \
       2>/home/astra/ASTRA/dialogScripts/text.txt