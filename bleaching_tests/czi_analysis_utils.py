#%%
import czifile
import napari
import re

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

from skimage import data
from skimage.filters import threshold_otsu
from skimage.segmentation import clear_border
from skimage.measure import label, regionprops
from skimage.morphology import closing, square
from skimage.color import label2rgb

import numpy as np
import pandas as pd

from glob import glob 
import os


#%%
def create_masks(image):
    
    thresh = threshold_otsu(image)
    bw = closing(image > thresh, square(3))
    
    # remove artifacts connected to image border
    cleared = clear_border(bw)
    # label image regions
    label_image = label(cleared)
    
    # Iterate over each region in the labeled image
    for region in regionprops(label_image):
        # Check if the region's area is below the threshold (500 pixels in this case)
        if region.area < 500:
            # Set the pixels corresponding to the current small region to 0 in the label_image
            label_image[label_image == region.label] = 0
    
    # Now, label_image will only contain regions with an area above 500 pixels
    # You may need to relabel to ensure continuous labels after removal
    label_image = label(label_image > 0)
    
    # to make the background transparent, pass the value of `bg_label`,
    # and leave `bg_color` as `None` and `kind` as `overlay`
    image_label_overlay = label2rgb(label_image, image=image, bg_label=0, kind='overlay')
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.imshow(image, cmap='gray')
    ax.imshow(image_label_overlay, alpha=0.5)
    
    ax.set_axis_off()
    plt.tight_layout()
    plt.show()
    
    return label_image

