import numpy as np
import matplotlib.pyplot as plt
import os

plt.rcParams.update({'font.size': 18})

path = os.getcwd()
print(path)

Dim = int(input('number of dimensions: '))
r_cutoff = float(input('RDF cutoff: '))

metode = ['mc', 'md']


def plot():
    for metoda in metode:
        r, gr, dgr = np.loadtxt(f'{metoda:s}_correlation').T
        r0, r1 = 0.0, r_cutoff
        plt.figure()
        plt.fill_between(r, gr+dgr, gr-dgr, alpha=0.25)
        plt.plot(r, gr, alpha=0.75)
        plt.title('Radial distribution function')
        plt.xlabel(r'$r^{*}$')
        plt.ylabel(r'$g(r^{*})$')
        plt.xlim(r0, r1)
        plt.tight_layout()
        plt.savefig(f'{metoda:s}_correlation.png', dpi=300)
        plt.close()


def correlation():
    os.chdir(f'{Dim:d}')
    temperature = [ f.path for f in os.scandir('.') if f.is_dir() ]
    print(temperature)
    for temperatura in temperature:
        os.chdir(temperatura)
        gostote = [ f.path for f in os.scandir('.') if f.is_dir() ]
        print(gostote)
        for gostota in gostote:
            os.chdir(gostota)
            plot()
            os.chdir('..')
        os.chdir('..')
    os.chdir('..')

correlation()
