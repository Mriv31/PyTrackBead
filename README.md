#PyTrackBead

This is a python wrapper of the XY tracking functions written by Vincent Croquette.

The tracking will run in a separate C++ Thread and analyze images as they are fed.



##Definitions of Buffer

You need to register one or several shared image buffers (unsigned char*) to the module as well as a buffer containing the number of frames received (int). When the number of images tracked by the module is smaller than the number of frames received, it will track the image.

This module will NOT modify any of the buffer given to him, so there cannot be a competition of writing between the writing module (for example the camera) and that buffer.

You can also define a timebuffer (uint64_t), to store the timestamps of the different images.

TYPES ARE IMPORTANT. If you create a numpy array as a buffer, if the type is wrong, COPIES will be made instead of reference passing.


##Saving
You can set up an output file and a saving mode with the methods "define_timebuffer" and the methods "set_saving". The current version will write in a text file but can be easily modified to write in a binary file.

You can get the last x tracked with the method get_x and the last y tracked with get_y. get_copy_x_array will return a copy of the whole buffer of x positions (which size is hard coded as TRACKING_BUFFER_SIZE)

## Installation of pre-built assembly

### Python3.10 Windows 64 bits
```pip install git+https://github.com/Mriv31/PyTrackBead/PyTrackBead-1.0-cp310-cp310-win_amd64.whl```

### Python3.10 Windows 32 bits
```pip install git+https://github.com/Mriv31/PyTrackBead/PyTrackBead-1.0-cp310-cp310-win32.whl```


## Build yourself
Shouldn't be too hard but you need to define a correct build environment and you will need pybind11. 
```pip install git+https://github.com/Mriv31/PyTrackBead```


## example
Install pyqtgraph ```pip install pyqtgraph```
Clone the repo : ```git clone https://github.com/Mriv31/PyTrackBead```
Make the Example directory as your current directory ```cd Examples```
Run the Example : ```python test.py```
