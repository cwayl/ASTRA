#!/bin/bash

# show an inputbox
dialog --title "NORAD ID" \
       --backtitle "ASTRA" \
       --nocancel \
       --inputbox "Enter the NORAD ID of the sattelte " 8 60 2>/home/astra/ASTRA/dialogScripts/text.txt

