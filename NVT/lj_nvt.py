from random import randint
import os
import shutil
import subprocess

path = os.getcwd()
print(path)

Ran = int()
Dim = int(input('number of dimensions: '))
Temp = float(input('temperature: '))
GostPoints = int(input('number of density points: '))
GostStart = float(input('starting density: '))
GostStep = float(input('size of density step: '))

print(72*"-")
print("PARAMETERS INITIALIZED")
print(72*"-")

shutil.copy('lj_nvt.exe', '{:d}'.format(Dim))
os.chdir('{:d}'.format(Dim))
shutil.copy('lj_nvt.exe', '{:.3f}'.format(Temp))
os.chdir('{:.3f}'.format(Temp))
for i in range(GostPoints):
    Ran = randint(1, 32000)
    Gost = GostStart + i * GostStep
    os.mkdir('{:.3f}'.format(Gost))
    shutil.copy('lj_nvt.exe', '{:.3f}'.format(Gost))
    os.chdir('{:.3f}'.format(Gost))
    with open(file='param', mode='w') as f:
        f.write("'ensemble parameters'\n")
        f.write("'number of dimensions (D)'                  " + '{:d}'.format(Dim) + '\n')
        f.write("'number of particles (N)'                   100\n")
        f.write("'density'                                   " + '{:.3f}'.format(Gost) + 'd0' + '\n')
        f.write("'temperature'                               " + '{:.3f}'.format(Temp) + 'd0' + '\n')
        f.write("'simulation parameters'\n")
        f.write("'LJ potential parameters'                   1.0d0           1.0d0\n")
        f.write("'random seed'                               " + str(Ran) + '\n')
        f.write("'number of series of sampling (MAX 20)'     5\n")
        f.write("'equilibration length'                      100000\n")
        f.write("'MC parameters'\n")
        f.write("'number of cycles (N trial moves)'          10000\n")
        f.write("'trial moves per sample (same as N)'        100\n")
        f.write("'MAX random displacement by sigma'          0.1d0\n")
        f.write("'MD parameters'\n")
        f.write("'number of time steps'                      10000\n")
        f.write("'time step'                                 0.001d0\n")
        f.write("'temperature time constant'                 0.1d0\n")
        f.write("'time steps per frame'                      100\n")
        f.write("'correlation function parameters'\n")
        f.write("'g(r) radial interval'                      0.02d0\n")
        f.write("'cycles/steps per g(r) sampling'            10\n")
    print(72*"-")
    print("RUN {:d}".format(i))
    print(72*"-")
    if os.name == 'posix':
        Process = subprocess.Popen('./lj_nvt.exe', shell=False)
    elif os.name == 'nt':
        Process = subprocess.Popen('.\lj_nvt.exe', shell=False)
    Process.wait()
    print(72*"-")
    print("CHECKPOINT {:d}".format(i))
    print(72*"-")
    os.chdir('..')
os.chdir('..')
os.chdir('..')
print("DONE")
