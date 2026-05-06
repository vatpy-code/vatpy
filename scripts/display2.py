# Packages --------------------------------------------------------------------
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Button
from PIL import Image


# Functions -------------------------------------------------------------------
def refresh(event):
    # Image update:
    img = Image.open(file)
    ax.imshow(img)
    ax.axis('off')
    plt.draw()


# Script ----------------------------------------------------------------------
print('\nWelcome to Vatpy Image Display')

# Image file from command line argument:
file = sys.argv[1]
if not file:
    print('  * No image provided')
    print('  * Exiting...')
    sys.exit()

# Load image and get ratio:
img = Image.open(file)
img_array = np.asarray(img)
ratio = np.shape(img_array)[0] / np.shape(img_array)[1]

# Create a figure:
fig, ax = plt.subplots(figsize=(6/ratio, 6), facecolor='gray')
fig.subplots_adjust(left=0, right=1, bottom=0, top=0.9)

# Display the image:
ax.imshow(img)
ax.axis('off')

# Create a refresh button:
button_ax = plt.axes([0.4, 0.93, 0.2, 0.05])
button = Button(button_ax, 'REFRESH', color='darkgray', hovercolor='lightgray')

# Connect the button to the refresh function (see above):
button.on_clicked(refresh)

# Create the GUI:
print(f'  * Displaying {file}')
try:
    show = True
    while show is True:
        plt.show()
        show = False
except KeyboardInterrupt:
    print('\n  * Program interupted by the user!')

# End statement:
print('  * Done!\n')

# End of file -----------------------------------------------------------------
