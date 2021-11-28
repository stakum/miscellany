#!/usr/bin/env  python3
import numpy as np
import sys
import argparse
import ase.build
import ase.io 

parser = argparse.ArgumentParser(description = 'Post-process code for Green-Kubo calcualtion')

# arguments 
parser.add_argument('-c', '--coordinate', help = 'Input coordinate file')
parser.add_argument('--style', default = 'atomic', help = 'LAMMPS style')
parser.add_argument('-t', '--temperature', type = float, default = 300.0, help = 'Temperature')
parser.add_argument('-s', '--sampling', type = int, default = 1, help ='Sampling interval of autocorrelation function')
parser.add_argument('-dt', '--timestep', type = float, default = 0.01, help = 'Time step of MD simualtion')
parser.add_argument('-en', '--ensemble', type = int, default = 1, help = 'Number of ensembles')
parser.add_argument('-f', '--file_header_name', help = 'File header name of autocorrelation function')

args = parser.parse_args()

# Print arguments
print('Input coordinate file = ' + args.coordinate)
print('Temperature (K) = ' + str(args.temperature))
print('Sampling interval = ' + str(args.sampling))
print('Time step (ps) = ' + str(args.timestep))


# Physical constants
# 'metal unit' is used. 
A2m = 1.0e-10 # Angstrom
ps2s = 1.0e-12 # Picosecond
kB = 1.380649e-23 # Boltzmann constant (J/K)
eV2J = 1.602176634e-19 # Electron Volt (J)
conversion = eV2J * eV2J / ( ps2s * A2m )

# Get volume from lammps data file
atoms = ase.io.read(args.coordinate, format = 'lammps-data', style = args.style)
volume = atoms.get_volume()

print('Volume (Angstrom^3) = ' + str(volume))

# Collect autocorrelation functions of heat flux

for ind in range(args.ensemble):
    fname = args.file_header_name + str(ind + 1)
    f = open(fname, 'r')

    for jnd in range(3):
        line = f.readline()

    if ind == 0:
        line = f.readline()
        lines = line.split()
        ncor = int(lines[1])
        tcor = np.zeros(ncor)
        JJcorr = np.zeros((3, ncor))
    else:
        line = f.readline()

    for jnd in range(ncor):
        line = f.readline()
        lines = line.split()
        
        if ind == 0:
            tcor[jnd] = int(lines[1]) * args.timestep

        for xyz in range(3):
            JJcorr[xyz][jnd] = JJcorr[xyz][jnd] + float(lines[xyz + 3]) / float(args.ensemble)

    f.close()

scale = conversion / (kB * args.temperature * args.temperature * volume) * args.sampling * args.timestep
kappa=np.zeros((3,ncor))
for ind in range(ncor):
    for xyz in range(3):
        kappa[xyz][ind] = np.trapz(JJcorr[xyz][0:ind]) * scale
 
f = open('JJcorrKappa.txt','w')
for ind in range(ncor):
    print(tcor[ind], JJcorr[0][ind], JJcorr[1][ind], JJcorr[2][ind], \
          kappa[0][ind], kappa[1][ind], kappa[2][ind], end = "\n", file = f)

f.close()
