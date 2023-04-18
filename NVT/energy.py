import numpy as np
import matplotlib.pyplot as plt
import os

plt.rcParams.update({'font.size': 18})

path = os.getcwd()
print(path)

Dim = int(input('number of dimensions: '))

mc_procesi = ['mc_minimization', 'mc_equilibration', 'mc_sam001', 'mc_sam005']
md_procesi = ['md_equilibration', 'md_sam001', 'md_sam005']


def plot():
    for proces in mc_procesi:
        t, u, p = np.loadtxt(f'{proces:s}').T
        t0, t1 = np.min(t), np.max(t)
        fig, ax = plt.subplots(2, 1, figsize=(15, 10))
        ax[0].plot(t, u, alpha=0.75)
        ax[0].set_title('Potencialna energija')
        ax[0].set_xlabel(r'$i$-ti korak')
        ax[0].set_ylabel(r'$U^{*}$')
        ax[0].set_xbound(t0, t1)
        ax[1].plot(t, p, alpha=0.75)
        ax[1].set_title('Tlak')
        ax[1].set_xlabel(r'$i$-ti korak')
        ax[1].set_ylabel(r'$F^{*}r^{*}$')
        ax[1].set_xbound(t0, t1)
        plt.tight_layout()
        fig.savefig(f'{proces:s}.png', dpi=300)
        plt.close(fig)
    for proces in md_procesi:
        t, u, ekin, etot, p, tem = np.loadtxt(f'{proces:s}').T
        t0, t1 = np.min(t), np.max(t)
        fig, ax = plt.subplots(5, 1, figsize=(15, 25))
        ax[0].plot(t, u, alpha=0.75)
        ax[0].set_title('Potencialna energija')
        ax[0].set_xlabel(r'$t^{*}$')
        ax[0].set_ylabel(r'$U^{*}$')
        ax[0].set_xbound(t0, t1)
        ax[1].plot(t, ekin, alpha=0.75)
        ax[1].set_title('Kinetična energija')
        ax[1].set_xlabel(r'$t^{*}$')
        ax[1].set_ylabel(r'$E_{kin}^{*}$')
        ax[1].set_xbound(t0, t1)
        ax[2].plot(t, etot, alpha=0.75)
        ax[2].set_title('Celotna notranja energija')
        ax[2].set_xlabel(r'$t^{*}$')
        ax[2].set_ylabel(r'$E_{tot}^{*}$')
        ax[2].set_xbound(t0, t1)
        ax[3].plot(t, p, alpha=0.75)
        ax[3].set_title('Tlak')
        ax[3].set_xlabel(r'$t^{*}$')
        ax[3].set_ylabel(r'$F^{*}r^{*}$')
        ax[3].set_xbound(t0, t1)
        ax[4].plot(t, tem, alpha=0.75)
        ax[4].set_title('Temperatura')
        ax[4].set_xlabel(r'$t^{*}$')
        ax[4].set_ylabel(r'$T^{*}$')
        ax[4].set_xbound(t0, t1)
        plt.tight_layout()
        fig.savefig(f'{proces:s}.png', dpi=300)
        plt.close(fig)


def energy():
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

energy()
