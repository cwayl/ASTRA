#!/bin/bash

# show an inputbox
dialog --title "Manual TLE" \
       --backtitle "ASTRA" \
       --nocancel \
       --inputbox "Enter line 2 of the TLE " 8 60 2>/home/astra/ASTRA/text.txt

