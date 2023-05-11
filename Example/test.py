from PyTrackBead.bindings import Tracker
from matplotlib import pyplot as plt
from PIL import Image
import numpy as np
import time
from pyqtgraph import QtCore
import pyqtgraph as pg
from CrossItem import CrossItem


im_buffer_size = 128


#Just define images for dummy demo. In reality they will be passed by a camera
imfolder = "14_16_04/"
im1 = np.asarray(Image.open(imfolder+"file0.tiff"))
height = im1.shape[0]
width = im1.shape[1]
print(im1.shape)

Tr=Tracker(height,width) #Tracker class fro PyTrackBead

#Create 128 arrays to store the images and share with the module.
#NOTE THE TYPE SPECIFICATION
#NEVER DEREFERENCE THE BUFFER IN PYTHON OR THEY WILL BE DESTROYED BY THE GARBAGE COLLECTOR
buffers = [np.zeros([height,width],dtype=np.ubyte) for i in range(im_buffer_size)]

#CREATE THE BUFFER containing the number of frames put and to set_frames_to_analyze
#NOTE THE TYPE SPECIFICATION
nb_frame_put = np.array([0],dtype=np.intc)

#CREATE THE BUFFER containing time stamps of each frame
#NOTE THE TYPE SPECIFICATION
time_buf = np.zeros([im_buffer_size],dtype=np.uint64)

#REGISTER THE BUFFERS
for i in range(im_buffer_size):
    Tr.addBuffer(buffers[i],buffers[i].size)
    #UNCOMMENT THE FOLLOWING LINE IF A CAMERA WRITE TO THAT Buffer
    #IT WILL MAKE THE NUMPY ARRAY READ ONLY FOR NUMPY
    #buffers[i].flags.writeable = False

Tr.define_timebuffer(time_buf)
Tr.define_image_counter(nb_frame_put,nb_frame_put.size)


cl = 32
cw = 8
#Add Bead With cross lengths and widths above
Tr.addBead(cl,cw,48,48)
Tr.addBead(cl,cw,62,22)


#Start tracking and saving
Tr.Start()
Tr.set_saving("output.txt") #Define Output for saving of the trajectories

#The tracking thread has started but is now waiting for the buffer to be feed
#and waiting for the image counter to increase
# Typically this will be performed by a camera thread to which the buffers are also passed









###Now feed the buffers for demo purposes
# THE FOLLOWING IS NOT INTERESTING PER SE
# It is just a DUMMY GRAPHICAL INTERFACE
ci = 0 #index in buffers
imnb = 0 #image file loaded


app = pg.mkQApp("ImageItem Example")

## Create window with GraphicsView widget
win = pg.GraphicsLayoutWidget()
win.show()  ## show widget alone in its own window
win.setWindowTitle('PyTrackBead Demo')
view = win.addViewBox()

## lock the aspect ratio so pixels are always square
view.setAspectLocked(True)

## Create image item
img = pg.ImageItem(border='w')
view.addItem(img)


ci = 0
imnb = 0
crosses = [CrossItem(cl,cw,48,48),CrossItem(cl,cw,22,62)]
for c in crosses:
    view.addItem(c)
first = 1
st = time.time()

def updateData():
    global timer,ci,imnb, first
    if first != 1:
        for i,c in enumerate(crosses):
            c.set_prop(cl,cw,Tr.get_y(i),Tr.get_x(i))

    time_buf[ci] = (time.time()-st)*1000
    buffers[ci][:,:] = np.asarray(Image.open(imfolder+"file"+str(imnb)+".tiff"))
    nb_frame_put[0] +=1
    ci+=1
    if (ci == im_buffer_size):
        ci = 0
    imnb+=1
    if (imnb == 2000): #Just because there are 2000 files in the demo
        imnb = 0
    timer.start(50)
    img.setImage(buffers[ci-1])
    first = 0




timer = QtCore.QTimer()
timer.setSingleShot(True)
timer.timeout.connect(updateData)
updateData()


if __name__ == '__main__':
    pg.exec()
