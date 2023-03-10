import numpy as np
import h5py
import matplotlib as mp
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits import mplot3d
import argparse
import os
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

# interval = 100  # in seconds

# DIR ="../../data/"

file_name = os.path.join('data', 'pop.pop.h5')  # rhoNeutral" #"P"

# ========== Figure Directory Setup =============
figPath = os.path.join(folder, "figures")  # DO NOT CHANGE THE PATH
if os.path.exists(figPath):
    print("figure directory exists. Existing figures will be replaced.")
else:
    os.mkdir(figPath)


h5 = h5py.File(os.path.join(folder, file_name), 'r')

# file = h5py.File('../../probe/mag_5e-5_64_128_832/data_10/pop.pop.h5','r')

n = 1500000
dp = 1000
# dp1 = 500

growthPeriod = 50000
interval1 = 100
interval2  = 200

if (n <= growthPeriod and n % interval1 == 0):
    data_num = np.arange(start=100.0, stop=n, step=dp)
elif (n % interval2 == 0):
    data_num = np.arange(start=200.0, stop=n, step=dp)



if (show_anim == True):
    def animate(i):
        # ======Species 1 Data=========
        pos1 = h5['/pos/specie 0/n=%3.1f' % (data_num[i])]
        vel1 = h5['/vel/specie 0/n=%3.1f' % (data_num[i]+0.5)]
        x1 = pos1[:, 0]
        y1 = pos1[:, 1]
        z1 = pos1[:, 2]
        vx1 = vel1[:, 0]
        vy1 = vel1[:, 1]
        vz1 = vel1[:, 2]

        # ======Species 2 Data=========
        pos2 = h5['/pos/specie 1/n=%3.1f' % (data_num[i])]
        vel2 = h5['/vel/specie 1/n=%3.1f' % (data_num[i]+0.5)]
        x2 = pos2[:, 0]
        y2 = pos2[:, 1]
        z2 = pos2[:, 2]
        vx2 = vel2[:, 0]
        vy2 = vel2[:, 1]
        vz2 = vel2[:, 2]

        # # ================== species 3 data =============
        pos3 = h5['/pos/specie 2/n=%3.1f' % (data_num[i])]
        vel3 = h5['/vel/specie 2/n=%3.1f' % (data_num[i]+0.5)]
        x3 = pos3[:, 0]
        y3 = pos3[:, 1]
        z3 = pos3[:, 2]
        vx3 = vel3[:, 0]
        vy3 = vel3[:, 1]
        vz3 = vel3[:, 2]

        # #  ================== species 3 data =============
        # pos4 = h5['/pos/specie 3/n=%3.1f' % (data_num[i])]
        # vel4 = h5['/vel/specie 3/n=%3.1f' % (data_num[i]+0.5)]
        # x4 = pos4[:, 0]
        # y4 = pos4[:, 1]
        # z4 = pos4[:, 2]
        # vx4 = vel4[:, 0]
        # vy4 = vel4[:, 1]
        # vz4 = vel4[:, 2]



        ax1.cla()
        if config:
            ax1.scatter(x1, y1, marker='.', color='b', alpha=1.0, s=1)
            ax1.scatter(x2, y2, marker='.', color='r', alpha=1.0, s=1)
        else:
            ax1.scatter(x1, vx1, marker='.', color='b', alpha=1.0, s=1)
            # ax1.scatter(x2, vx2, marker='.', color='r', alpha=1.0, s=1)
        ax1.set_title('Species 1 Phase Space ( cold electrons)(TimeSteps = %d' % (i*dp)+')')
        ax1.set_xlabel("$x$")
        ax1.set_ylabel("$v_x$")
        if both:
            ax1.set_ylim([np.min(np.array([vx1,vx2])), np.max(np.array([vx1,vx2]))])
        if not both:
            ax2.cla()
            if config:
                ax2.scatter(x2, y2, marker='.' ,color='r' ,alpha=1.0, s=1)
            else:
                ax2.scatter(x2, vx2, marker='.' ,color='r' ,alpha=1.0, s=1)
            ax2.set_title('Species 2 Phase Space (Precipitating electrons)(TimeSteps = %d'%(i*dp)+')')
            ax2.set_xlabel("$x$")
            ax2.set_ylabel("$v_x$")
        ax3.cla()
        if config:
            ax3.scatter(x3, y3, marker='.', color='b', alpha=1.0, s=1)
            ax3.scatter(x2, y2, marker='.', color='r', alpha=1.0, s=1)
        else:
            ax3.scatter(x3, vx3, marker='.', color='b', alpha=1.0, s=1)
            # ax1.scatter(x2, vx2, marker='.', color='r', alpha=1.0, s=1)
            ax3.set_title('Species 3 Phase Space (ions ) (TimeSteps = %d' % (i*dp)+')')
            ax3.set_xlabel("$x$")
            ax3.set_ylabel("$v_x$")
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
    fig,(ax1, ax2,ax3) = plt.subplots(3, 1, figsize=figsize/25.4, constrained_layout = True, dpi=ppi)

if (n <= growthPeriod and n % interval1 == 0):
    ani = animation.FuncAnimation(fig, animate, frames=len(data_num), interval=interval1, blit=False)
elif (n % interval2 == 0):
    ani = animation.FuncAnimation(fig, animate, frames=len(data_num), interval=interval2, blit=False)

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