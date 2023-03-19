import numpy as np
import h5py
import argparse
import os
import re
from tqdm import tqdm
#========= Configuration ===========


parser = argparse.ArgumentParser(description='Population DATA compressor PICSP')
parser.add_argument('-f', '--folder', default="../../", type=str, help='path to data file')
parser.add_argument('-s', '--species', default=2, type=int, help='number of species')
args = parser.parse_args()


folder = args.folder
species = args.species


file_name = os.path.join('data', 'pop.pop.h5')  # rhoNeutral" #"P"
h5 = h5py.File(os.path.join(folder, file_name), 'r')


def sort_by_number(name):
    return float(name.split('=')[1])


for p in tqdm(range(species), desc='Loading species data', unit='species'):
    pos = []
    vel = []
    pos_group = h5[f'/pos/specie {p}/']
    vel_group = h5[f'/vel/specie {p}/']
    # data_len = len(pos_group.keys())

    sorted_pos_names = sorted(pos_group, key=sort_by_number)
    sorted_vel_names = sorted(vel_group, key=sort_by_number)
    time_list = [sort_by_number(name) for name in sorted_pos_names]

    for d_name in tqdm(sorted_pos_names, desc='Processing position data', leave=False):
        # Access the dataset and extract its data as a NumPy array
        pos_dataset = pos_group[d_name]
        # Append the array to the list of data arrays
        pos.append(np.array(pos_dataset))
        tqdm.write('Processed species position data {}'.format(d_name), end='\r')
    for d_name in tqdm(sorted_vel_names, desc='Processing velocity data', leave=False):
        # Access the dataset and extract its data as a NumPy array
        vel_dataset = vel_group[d_name]
        # Append the array to the list of data arrays
        vel.append(np.array(vel_dataset))
        tqdm.write('Processed species velocity data {}'.format(d_name), end='\r')
    print(np.array(pos).shape)
    np.savez_compressed(os.path.join(folder, f'pop_species{p}.npz'), time=np.array(time_list), pos=np.array(pos), vel=np.array(vel))
    # Save the data to the CSV file
    np.savetxt("time.csv", np.array(time_list), delimiter=",")
