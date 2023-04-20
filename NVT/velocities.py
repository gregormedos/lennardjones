import numpy as np
import matplotlib.pyplot as plt
import os

plt.rcParams.update({'font.size': 18})

path = os.getcwd()
print(path)

Dim = int(input('number of dimensions: '))

components = [r'$v_{x}$', r'$v_{y}$', r'$v_{z}$']
distributions = [r'$P(v_{x})$', r'$P(v_{y})$', r'$P(v_{z})$']


def plot():
    velocities = np.loadtxt('md_velocities0001')
    for i in range(1,20):
        velocities = np.vstack((velocities, np.loadtxt(f'md_velocities{i:04d}')))
        magnitudes = np.sqrt(np.sum(velocities**2, axis=1))
    fig, ax = plt.subplots(Dim+1, 1, figsize=((Dim+1)*5, 10))
    for i in range(Dim):
        ax[i].hist(velocities[:, i], 50, density=True, alpha=0.75)
        ax[i].set_xlabel(f'{components[i]:s}')
        ax[i].set_ylabel(f'{distributions[i]:s}')
    ax[Dim].hist(magnitudes, 50, density=True, alpha=0.75)
    ax[Dim].set_xlabel(r'$v$')
    ax[Dim].set_xlabel(r'$p(v)$')
    plt.title('Velocity distribution')
    plt.tight_layout()
    fig.savefig('md_velocities.png', dpi=300)
    plt.close(fig)


def velocity():
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

velocity()
