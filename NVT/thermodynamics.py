import pandas as pd
import matplotlib.pyplot as plt
import os

plt.rcParams.update({'font.size': 18})

path = os.getcwd()
print(path)

Dim = int(input('number of dimensions: '))

metode = ['mc', 'md']

kolicine = ['T', 'rho', 'U/N', 'Err U/N', 'p', 'Err p', 'C_V/N', 'Err C_V/N']


def thermodynamics(metoda):
    os.chdir(f'{Dim:d}')
    with open(f'{metoda:s}_thermodynamics_annotated', 'w') as fw:
        for kolicina in kolicine:
            fw.write(f'{kolicina:s}'.rjust(16))
        fw.write('\n')
        temperature = [ f.path for f in os.scandir('.') if f.is_dir() ]
        print(temperature)
        for temperatura in temperature:
            os.chdir(temperatura)
            gostote = [ f.path for f in os.scandir('.') if f.is_dir() ]
            print(gostote)
            for gostota in gostote:
                os.chdir(gostota)
                with open(f'{metoda:s}_thermodynamics', 'r') as fr:
                    for line in fr:
                        line = line.strip('\n')
                        fw.write(line)
                    fw.write('\n')
                os.chdir('..')
            os.chdir('..')    
    data = pd.read_fwf(f'{metoda:s}_thermodynamics_annotated')
    temperature = set(data['T'])
    fig, ax = plt.subplots(3, 1, figsize=(15, 15))
    for temperatura in temperature:
        temp, dens, u, du, p, dp, cv, dcv = data.loc[data['T'] == temperatura].to_dict('list').values()
        ax[0].errorbar(dens, u, du, marker='o', markeredgecolor='black', alpha=0.75, label=r'$T^{*}=$'+f'{temperatura:.3f}')
        ax[0].set_title('Gostota presežne notranje energije')
        ax[0].set_xlabel(r'$\rho^{*}$')
        ax[0].set_ylabel(r'$u^{*}=U^{*}/N$')
        ax[0].legend()
        ax[1].errorbar(dens, p, dp, marker='o', markeredgecolor='black', alpha=0.75, label=r'$T^{*}=$'+f'{temperatura:.3f}')
        ax[1].set_title('Tlak')
        ax[1].set_xlabel(r'$\rho^{*}$')
        ax[1].set_ylabel(r'$p^{*}$')
        ax[1].legend()
        ax[2].errorbar(dens, cv, dcv, marker='o', markeredgecolor='black', alpha=0.75, label=r'$T^{*}=$'+f'{temperatura:.3f}')
        ax[2].set_title('Gostota presežne toplotne kapacitete')
        ax[2].set_xlabel(r'$\rho^{*}$')
        ax[2].set_ylabel(r'$c_{V}^{*}=C_{V}^{*}/N$')
        ax[2].legend()
    plt.tight_layout()
    fig.savefig(f'{metoda:s}_thermodynamics.png', dpi=300)
    plt.close(fig)
    os.chdir('..')

for metoda in metode:
    thermodynamics(metoda)
