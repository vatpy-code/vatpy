import sys
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image

print('\n Welcome to VATPY Image Display (by Jonathan Petersson)')

# File:
file = str(sys.argv[1])
print(f'  * Displaying {file}')

# Display:
img = Image.open(file)
img_array = np.asarray(img)
ratio = np.shape(img_array)[0] / np.shape(img_array)[1]
plt.figure(figsize=(6/ratio, 6), layout='constrained')
plt.imshow(img)
plt.axis('off')
plt.show()

print('  * Done!\n')
