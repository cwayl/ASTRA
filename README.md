# ASTRA
Active Sattelite Radio Tracking Antenna(ASTRA) is a software for automatic control of a rotator. 
## Table of Contents
* [Installation Instructions](#installation-instructions)
* [Running the Code](#running-the-code)


## Installation Instructions

0. Install [Raspberry Pi OS Legacy](https://www.raspberrypi.com/software/operating-systems/) (bullseye) 

1. Remove python3 
   ```
   sudo apt-get remove python3
   ```
   ```
   sudo apt-get remove autoremove
   ```
   ```
   sudo apt-get remove autoclean
   ```

2. Follow the installation instructions of the [gasctld driver](https://github.com/SmallSatGasTeam/greenctld)

3. Download ASTRA with the command `git clone -b Socket-Interface https://github.com/cwayl/ASTRA.git`


4. Run `gcc -o ASTRA_Program ASTRA/main.c -lm` to compile the code.

## Running the Code
1. Turn on the rotator driver by running `./gasctld/greenctld --az-device /dev/ttyUSB0 --el-device /dev/ttyUSB1` Remember that the USB adresses may be different in some cases. 
2. Turn on ASTRA by running `./ASTRA_Program` in a seperate terminal. 

