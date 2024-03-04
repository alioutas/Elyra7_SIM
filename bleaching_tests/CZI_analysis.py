#%%
import czifile
import napari

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


# %%
# find all images in the folder
img_files = glob('/Volumes/Genetics/Wu_Lab-Vutara/Experiments/Eunice/Elyra_Eunice/WGI/bleaching_test/2xSSCT_4channels_WGI2.0_20240223/*.czi')

# from the file name select any character before the first underscore
# this will be the name of the image
img_names = [img_file.split('/')[-1].split('_')[0] for img_file in img_files]
# find unique image names and order them
img_names = sorted(list(set(img_names)))

channels = sorted(list(set([img_file.split('/')[-1].split('_')[1] for img_file in img_files])))




#%%
img = czifile.imread('/Volumes/Genetics/Wu_Lab-Vutara/Experiments/Eunice/Elyra_Eunice/WGI/bleaching_test/2xSSCT_4channels_WGI2.0_20240223/001_405_2xSSCT_buffer_bleachingtest_2024_02_23__14_06_09.czi')
print(img.shape)
# %%
img_ = img[0,0,0,:,0,:,:,0]
img_.reset_index(drop=True, inplace=True)
img_.shape

#%%
img_100 = img_[:100,:,:]


# %%
viewer, image_layer = napari.imshow(img_100, rgb=False)

# %%
# theshold my first image to get a mask with otsu

#%%
image = img[0,0,0,0,0,:,:,0]
thresh = threshold_otsu(image)
bw = closing(image> thresh, square(3))

# %%
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
image_label_overlay = label2rgb(label_image, image=image, bg_label=0)
fig, ax = plt.subplots(figsize=(10, 6))
ax.imshow(image_label_overlay)

for region in regionprops(label_image):
    # take regions with large enough areas
    if region.area >= 100:
        # draw rectangle around segmented coins
        minr, minc, maxr, maxc = region.bbox
        rect = mpatches.Rectangle((minc, minr), maxc - minc, maxr - minr,
                                  fill=False, edgecolor='red', linewidth=2)
        ax.add_patch(rect)

ax.set_axis_off()
plt.tight_layout()
plt.show()
# %%

#Calculate fluorescence intensity for each mask for every time point in img_
time_point = []
mask = []
avg_int = []
median_int = []
channel = []



for i in range(img_.shape[0]):  # Assuming you're iterating over the first image for now; replace 1 with img_.shape[0] for all images
    intensity_image = img_[i, :, :]
    
    # It's important to pass intensity_image to regionprops to calculate intensity metrics
    for region in regionprops(label_image, intensity_image=intensity_image):
        # Extract the region's intensity values using its mask
        region_intensity_values = intensity_image[region.coords[:, 0], region.coords[:, 1]]
        
        # Append the data to the lists
        time_point.append(i)
        mask.append(region.label)
        avg_int.append(np.mean(region_intensity_values))
        median_int.append(np.median(region_intensity_values))

# %%
# normalize the intensity values min max for each mask
df_int_time = pd.DataFrame({'time_point': time_point, 
              'mask': mask, 
              'avg_int': avg_int, 
              'median_int': median_int})
# df_int_time['norm_avg_int'] = (df_int_time['avg_int'] - df_int_time.groupby('mask').avg_int.transform('min')) / (df_int_time.groupby('mask').avg_int.transform('max') - df_int_time.groupby('mask').avg_int.transform('min'))
# df_int_time['norm_median_int'] = (df_int_time['median_int'] - df_int_time.groupby('mask').median_int.transform('min')) / (df_int_time.groupby('mask').median_int.transform('max') - df_int_time.groupby('mask').median_int.transform('min'))
df_int_time['norm_avg_int'] = (df_int_time['avg_int'] - df_int_time.groupby('mask').avg_int.transform('min')) / (df_int_time.groupby('mask').avg_int.transform('max') - df_int_time.groupby('mask').avg_int.transform('min'))
df_int_time['norm_median_int'] = (df_int_time['median_int'] - df_int_time.groupby('mask').median_int.transform('min')) / (df_int_time.groupby('mask').median_int.transform('max') - df_int_time.groupby('mask').median_int.transform('min'))

# %%
# plot the normalized intensity values
fig, ax = plt.subplots(figsize=(10, 6))
for mask, group in df_int_time.groupby('mask'):
    ax.plot(group.time_point, group.norm_avg_int, label=mask)
    ax.legend(title = 'Mask')
plt.show()
# %%
