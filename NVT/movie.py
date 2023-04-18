import os
import shutil
import subprocess

path = os.getcwd()
print(path)

Dim = int(input('number of dimensions: '))


def movie() :
    snapshot_script = f'snapshot{Dim:d}.plt'
    frame_script = f'frame{Dim:d}.plt'
    shutil.copy(snapshot_script, f'{Dim:d}')
    shutil.copy(frame_script, f'{Dim:d}')
    os.chdir(f'{Dim:d}')
    temperature = [ f.path for f in os.scandir('.') if f.is_dir() ]
    print(temperature)
    for temperatura in temperature:
        shutil.copy(snapshot_script, temperatura)
        shutil.copy(frame_script, temperatura)
        os.chdir(temperatura)
        gostote = [ f.path for f in os.scandir('.') if f.is_dir() ]
        print(gostote)
        for gostota in gostote:
            shutil.copy(snapshot_script, gostota)
            shutil.copy(frame_script, gostota)
            os.chdir(gostota)
            Process = subprocess.Popen(f'gnuplot {snapshot_script:s}', shell=True)
            Process.wait()
            Process = subprocess.Popen(f'gnuplot {frame_script:s}', shell=True)
            Process.wait()
            os.chdir('..')
        os.chdir('..')
    os.chdir('..')

movie()
