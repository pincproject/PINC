import numpy as np
import h5py
import matplotlib as mp
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits import mplot3d
import argparse
import os
import glob
# import plotly.graph_objects as go
#========= Configuration ===========


parser = argparse.ArgumentParser(description='Phase - Space Data Animator PICSP')
parser.add_argument('-f', '--folder', default="../../", type=str, help='path to data file')
parser.add_argument('-a', '--animate', default=True, type=bool, help='Show Animation (True/False)')
parser.add_argument('-s', '--save', action='store_true', help='Add if want to Save Animation (True/False)')
parser.add_argument('-c', '--config', action='store_true', help='Add if want to plot config space')
parser.add_argument('-b', '--both', action='store_true', help='Plot both species together')
args = parser.parse_args()

show_anim = args.animate
save_anim = args.save
folder = args.folder
config = args.config
both = args.both

def sort_by_name(file_path):
    return file_path.split('/')[-1]

file_name = os.path.join('data', 'pop.pop.h5')  # rhoNeutral" #"P"

# ========== Figure Directory Setup =============
figPath = os.path.join(folder, "figures")  # DO NOT CHANGE THE PATH
if os.path.exists(figPath):
    print("figure directory exists. Existing figures will be replaced.")
else:
    os.mkdir(figPath)


h5 = h5py.File(os.path.join(folder, file_name), 'r')

npz_files = sorted(glob.glob(os.path.join(folder, '*.npz')), key=sort_by_name)
print(npz_files)
time = np.genfromtxt("time.csv", delimiter=",")
num_species = len(npz_files)

species_pos = []
species_vel = []


for i, file in enumerate(npz_files):
    data = np.load(file)
    species_pos.append(data["pos"])
    species_vel.append(data["vel"])

species_pos = np.array(species_pos)
species_vel = np.array(species_vel)
# Print the keys in the file
# print(data.keys())


if (show_anim == True):
    def animate(i):
        ax1.cla()
        if config:
            ax1.scatter(species_pos[0,i,:,0], species_pos[0,i,:,1], marker='.', color='b', alpha=1.0, s=1)
            # ax1.scatter(x2, y2, marker='.', color='r', alpha=1.0, s=1)
        else:
            ax1.scatter(species_pos[0,i,:,0], species_vel[0,i,:,0], marker='.', color='b', alpha=1.0, s=1)
            # ax1.scatter(x2, vx2, marker='.', color='r', alpha=1.0, s=1)
        ax1.set_title(f'Species 1 Phase Space ( cold electrons)(TimeSteps = {time[i]})')
        ax1.set_xlabel("$x$")
        ax1.set_ylabel("$v_x$")
        if both:
            ax1.set_ylim([np.min(np.array([species_vel[0,i,:,0],species_vel[1,i,:,0]])), np.max(np.array([species_vel[0,i,:,0],species_vel[1,i,:,0]]))])
        if not both:
            ax2.cla()
            if config:
                ax2.scatter(species_pos[1,i,:,0], species_pos[1,i,:,1], marker='.' ,color='r' ,alpha=1.0, s=1)
            else:
                ax2.scatter(species_pos[1,i,:,0], species_vel[1,i,:,0], marker='.' ,color='r' ,alpha=1.0, s=1)
            ax2.set_title(f'Species 2 Phase Space (Precipitating electrons)(TimeSteps = {time[i]})')
            ax2.set_xlabel("$x$")
            ax2.set_ylabel("$v_x$")
        # ax3.cla()
        # if config:
        #     ax3.scatter(x3, y3, marker='.', color='b', alpha=1.0, s=1)
        #     ax3.scatter(x2, y2, marker='.', color='r', alpha=1.0, s=1)
        # else:
        #     ax3.scatter(x3, vx3, marker='.', color='b', alpha=1.0, s=1)
        #     # ax1.scatter(x2, vx2, marker='.', color='r', alpha=1.0, s=1)
        #     ax3.set_title('Species 3 Phase Space (ions ) (TimeSteps = %d' % (i*dp)+')')
        #     ax3.set_xlabel("$x$")
        #     ax3.set_ylabel("$v_x$")
        # ax4.cla()
        # if config:
        #     ax4.scatter(x4, y4, marker='.', color='b', alpha=1.0, s=1)
        #     ax4.scatter(x2, y2, marker='.', color='r', alpha=1.0, s=1)
        # else:
        #     ax4.scatter(x4, vx4, marker='.', color='b', alpha=1.0, s=1)
        #     # ax1.scatter(x2, vx2, marker='.', color='r', alpha=1.0, s=1)
        #     ax4.set_title('Species 4 Phase Space (Ions ) (TimeSteps = %d' % (i*dp)+')')
        #     ax4.set_xlabel("$x$")
        #     ax4.set_ylabel("$v_x$")
        # ax2.set_xlim([0, Lx])
        # ax2.set_xlim([0, Ly])

        # ax1.set_ylim([-1, 1])
        # ax1.set_zlim([-Lz, Lz])


##### FIG SIZE CALC ############
figsize = np.array([200,200/1.618]) #Figure size in mm
dpi = 300                         #Print resolution
ppi = np.sqrt(1920**2+1200**2)/24 #Screen resolution

mp.rc('text', usetex=False)
mp.rc('font', family='sans-serif', size=10, serif='Computer Modern Roman')
mp.rc('axes', titlesize=10)
mp.rc('axes', labelsize=10)
mp.rc('xtick', labelsize=10)
mp.rc('ytick', labelsize=10)
mp.rc('legend', fontsize=10)


if both:
    fig, ax1 = plt.subplots(figsize=figsize/25.4, constrained_layout=True, dpi=ppi)
else:
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figsize/25.4, constrained_layout = True, dpi=ppi)

ani = animation.FuncAnimation(fig, animate, frames=len(time), interval=1, blit=False)

if (show_anim == True):
    plt.show()
if (save_anim ==True):
    try:
        Writer = animation.writers['ffmpeg']
        if (n <= growthPeriod and n % interval1 == 0):
            writer = Writer(fps=(1/interval1), metadata=dict(artist='Me'), bitrate=1800)
        elif (n % interval2 == 0):
            writer = Writer(fps=(1/interval2), metadata=dict(artist='Me'), bitrate=1800)
    except RuntimeError:
        print("ffmpeg not available trying ImageMagickWriter")
        if (n <= growthPeriod and n % interval1 == 0):
            writer = animation.ImageMagickWriter(fps=(1/interval1))
        elif (n % interval2 == 0):
            writer = animation.ImageMagickWriter(fps=(1/interval2))
        
    print("Saving movie to "+figPath+"/. Please wait .....")
    ani.save(figPath+'/phasespace_animation_PICSP.mp4')