#!/bin/bash

# show an menu
dialog --title "Start Menu" \
       --backtitle "ASTRA" \
       --nocancel \
       --menu "Select" 12 45 25 \
       1 "Start ASTRA" \
       2 "Settings Menu" \
       2>/home/astra/ASTRA/text.txt