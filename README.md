# ASTRA
Active Sattelite Radio Tracking Antenna(ASTRA) is a software for automatic control of a rotator. 
## Table of Contents
* [Installation Instructions](#installation-instructions)


## Installation Instructions

0. Install [Raspberry Pi OS Legacy](https://www.raspberrypi.com/software/operating-systems/) (bullseye) 

1. Remove python3 
   ```python
   sudo apt-get remove python3
   ```
   ```python
   sudo apt-get remove autoremove
   ```
   ```python
   sudo apt-get remove autoclean
   ```

2. Follow the installation instructions of the [gasctld driver](https://github.com/SmallSatGasTeam/greenctld)
