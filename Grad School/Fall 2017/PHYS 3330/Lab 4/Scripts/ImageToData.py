from scipy import ndimage
import numpy as np
from scipy import misc
import matplotlib.pyplot as plt

file_name = '80umSS.jpg'

# Import image
# The mode is an import, see:
# https://docs.scipy.org/doc/scipy-0.19.1/reference/generated/scipy.misc.imread.html
# Basically, we just care about the intensities, so the 'L' mode reads each pixel as a
# single value (0-255 8 bit - black & white) as opposed to three values (RGB - color)
image = misc.imread(file_name,mode='L')

# Rotate image
# This seems a little bit easier than making a diagonal slice. 
# Ideally, the image would be taken in the correct orientation...
rotated = ndimage.rotate(image, -6)

# The image is a ny X nx matrix, just pick the values.
# Plot the image and read off pixel values for the corners of
# a rectangle of interest.

# The rows correspond to y values and the columns correspond to x values
# this is counter intuitive to me, since I want to put the x's in the first place
cropped = rotated[630:692, 380:2260]

plt.figure()
plt.imshow(cropped)

# You can average across a slice to improve the signal to noise.
# axis = 0 => average columns individually
# axis = 1 => average rows individually

plt.figure()
plt.plot(cropped.mean(axis=0))
plt.show()

