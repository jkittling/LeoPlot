import h6py
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.image as image
import numpy as np
import glob
import os

##################
#Import File Names
##################
snaps_fn = sorted(glob.glob("../output/snap_*.hdf5"))
#gives list of all output files according to order (000->###)

##########
#Functions
##########
def plotsnap(fn, ax, rr, axis):
    """
    Plots and saves image of a single snap file
    """
    part = h5py.File(fn,'r')

    BoxSize = part['Header'].attrs['BoxSize']
    center=np.array([0.5*BoxSize,0.5*BoxSize,0.5*BoxSize])

   # pos = part['PartType1/Coordinates'][:]

    pos_dm = part['PartType1/Coordinates'][:]
    pos_bstar = part["PartType2/Coordinates"][:]
    pos_dstar = part["PartType3/Coordinates"][:]
    pos_star = np.append(pos_bstar, pos_dstar, axis=0)

    print(pos_dm)
    print(pos_bstar)
    print(pos_star)

   # ax.hexbin(pos[:,axis[0]] ,pos[:,axis[1]], bins='log', gridsize=200, alpha=0.4,
   #         extent = (center[axis[0]]-rr, center[axis[1]]+rr, center[axis[0]]-rr, center[axis[1]]+rr),
   #         cmap = "Greys_r", mincnt=2)

    ax[0].hexbin(pos_dm[:,axis[0]], pos_dm[:,axis[1]], bins='log', gridsize=200, alpha=0.4,
            extent = (center[axis[0]]-rr, center[axis[1]]+rr, center[axis[0]]-rr, center[axis[1]]+rr),
            cmap = "Greys_r", mincnt=2)

    ax[1].hexbin(pos_star[:,axis[0]] ,pos_star[:,axis[1]], bins='log', gridsize=200, alpha=0.4,
            extent = (center[axis[0]]-rr, center[axis[1]]+rr, center[axis[0]]-rr, center[axis[1]]+rr),
            cmap = "Oranges_r", mincnt=2)


    im_name = f"{fn[:-5]}.png"
    fig.savefig(im_name)
    #return None
    return im_name

def init():
    """
    Initial state of axes image object
    """
    axim.set_data(np.zeros((1000,1000)))

    return axim

def animate(i):
    """
    Function that updates the axes image
    with the ith snap plot
    """
    fn = snaps_fn[i][:-5]+".png"
    im_file = image.imread(fn)[-1::-1]
    axim.set_data(im_file)

def animate_list(i):
    im_file = image.imread(snaps_imn[i])
    axim.set_data(im_file)
##########
#Animating
##########
#First thing is to make a figure object
fig, axs = plt.subplots(1,2)
plt.axis("off")

#Make an image file for each snap file
for snap in snaps_fn:
    plotsnap(snap, axs, 800, [0,1])

#Make an axis image object
axim_dm = axs[0].imshow(np.zeros((100, 100)), origin='lower', alpha=1.0, zorder=1, aspect=1)
#axim_star = axs[0].imshow(np.zeros((100, 100)), origin='lower', alpha=1.0, zorder=1, aspect=1)

#Number of frames
nframes = len(snaps_fn)

#Set up video saver file format
writer = animation.ImageMagickWriter()

#Animate those guys
ani_dm = animation.FuncAnimation(fig, animate, frames = nframes,  interval = 30)
ani.save("snapanimation_dm.gif", writer=writer)
os.system("display snapanimation_dm.gif")

#ani = animation.FuncAnimation(fig, animate, frames = nframes,  interval = 30)
#ani.save("snapanimation_star.gif", writer=writer)
#os.system("display snapanimation_star.gif")

os.system("rm ../output/*.png")

