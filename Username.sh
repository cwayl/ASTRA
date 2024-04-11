#!/bin/bash

# show an inputbox
dialog --title "Space-Track.com Login Information" \
       --backtitle "ASTRA" \
       --nocancel \
       --inputbox "Enter your Space-Track.org username " 8 60 "$1" 2>/home/astra/ASTRA/text.txt

