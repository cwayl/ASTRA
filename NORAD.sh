#!/bin/bash

# show an inputbox
dialog --title "NORAD ID" \
       --backtitle "ASTRA" \
       --nocancel \
       --inputbox "Enter the NORAD ID of the satellite " 8 60 2>/home/astra/ASTRA/text.txt

