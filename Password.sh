#!/bin/bash

# show an inputbox
dialog --title "Space-Track.org Login Information" \
       --backtitle "ASTRA" \
       --nocancel \
       --inputbox "Enter your Space-Track.org password " 8 60 "$1" 2>/home/astra/ASTRA/text.txt

