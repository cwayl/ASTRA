#!/usr/bin/bash
echo "P 180 0" | nc localhost 4533 &> /dev/null &
disown
exit

