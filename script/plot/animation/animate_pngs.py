from numpy import *
from pylab import *
import matplotlib.pyplot as plt 
import matplotlib.image as mgimg
from matplotlib import animation

## Read in graphs


startindex = 100
stopindex = 40000
step = 100
p = startindex

myimages = []
fig = plt.figure()

for k in range(startindex, stopindex, step):

    ## Read in picture
    print(k)
    fname = "dt%i.png" %k
    img = mgimg.imread(fname)
    imgplot = plt.imshow(img)

    # append AxesImage object to the list
    myimages.append([imgplot])

## create an instance of animation
my_anim = animation.ArtistAnimation(fig, myimages, interval=100, blit=True, repeat_delay=1000)

## NB: The 'save' method here belongs to the object you created above
#my_anim.save("animation.mp4")

## Showtime!
plt.show()
