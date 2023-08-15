from multiprocessing import Pool
import subprocess
import os
import time

cpus = 80

structure_path = '/storage/evagsm/nobackup/crystal_dataset/structures/AF2_multi/L/'
output_path = '/storage/evagsm/nobackup/crystal_dataset/structures/AF2_multi_L_h/'
structures = os.listdir(structure_path)

def run(structure):
    print(structure)
    cmd = './reduce -FLIP {} > {}'.format(structure_path + structure, output_path + structure)
    subprocess.call(cmd, shell = True) 
      
if __name__ == '__main__':
    pool = Pool(processes = cpus)
    results = pool.map(run, structures)
    pool.close()
    pool.join()
