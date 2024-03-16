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

from czi_analysis_utils import create_masks
from tqdm import tqdm




# %%
# find all images in the folder
condition = 'SSCT2X'
path = "/Volumes/Genetics/Wu_Lab-Vutara/Experiments/Eunice/Elyra_Eunice/WGI/bleaching_test/2xSSCT_4channels_WGI2.0_20240223/"
img_files = glob(os.path.join(path, '*.czi'))

# from the file name select any character before the first underscore
# this will be the name of the image
img_names = [img_file.split('/')[-1].split('_')[0] for img_file in img_files]
# find unique image names and order them
img_names = sorted(list(set(img_names)))

channels = sorted(list(set([img_file.split('/')[-1].split('_')[1] for img_file in img_files])))



#%%
#Calculate fluorescence intensity
time_point = []
mask = []
avg_int = []
median_int = []
ch = []
im_name = []

# loop over all images and channels
for img_name in tqdm(img_names):
    # loop over all channels
    for channel in channels:
        print(img_name, channel)
        # find the img_name+'_'+channel+'_' within the img_files
        img_file = [img_file for img_file in img_files if img_name+'_'+channel+'_' in img_file]
        if len(img_file) == 0:
            next
        else:
            # find the 405 channel, create masks, and quantify the mean intensity
            if re.search('405', img_file[0]):
                print(img_file)
                # read image
                img = czifile.imread(img_file[0])
                img_ = img[0,0,0,:,0,:,:,:] # last dimension is channel
                image = img[0, 0, 0, 0, 0, :, :, 0]
                thresh = threshold_otsu(image)
                bw = closing(image > thresh, square(3))
                
                # remove artifacts connected to image border
                cleared = clear_border(bw)
                # label image regions
                label_image = label(cleared)
                
                # Find regions < 500 pixels and remove them
                for region in regionprops(label_image):
                    if region.area < 500:
                        label_image[label_image == region.label] = 0
                
                label_image = label(label_image > 0)
                
                # plot image and masks
                image_label_overlay = label2rgb(label_image, image=image, bg_label=0, kind='overlay')
                fig, ax = plt.subplots(figsize=(10, 6))
                ax.imshow(image, cmap='gray')
                ax.imshow(image_label_overlay, alpha=0.5)
                ax.set_axis_off()
                plt.tight_layout()
                plt.show()

                # loop over all time points
                for t in range(img_.shape[0]): 
                    intensity_image = img_[t, :, :]
                    
                    # It's important to pass intensity_image to regionprops to calculate intensity metrics
                    for region in regionprops(label_image, intensity_image=intensity_image):
                        # Extract the region's intensity values using its mask
                        region_intensity_values = intensity_image[region.coords[:, 0], region.coords[:, 1]]
                        
                        # Append the data to the lists
                        time_point.append(t)
                        mask.append(region.label)
                        avg_int.append(np.mean(region_intensity_values))
                        median_int.append(np.median(region_intensity_values))
                        ch.append(channel)
                        im_name.append(img_name)
            else:
                #use the mask created from the 405 channel to quantify the mean intensity of the other channels                # read image
                img = czifile.imread(img_file[0])
                img_ = img[0,0,0,:,0,:,:,:] # last dimension is channel
                # loop over all time points
                for t in range(img_.shape[0]): 
                    intensity_image = img_[t, :, :]
                    
                    # It's important to pass intensity_image to regionprops to calculate intensity metrics
                    for region in regionprops(label_image, intensity_image=intensity_image):
                        # Extract the region's intensity values using its mask
                        region_intensity_values = intensity_image[region.coords[:, 0], region.coords[:, 1]]
                        
                        # Append the data to the lists
                        time_point.append(t)
                        mask.append(region.label)
                        avg_int.append(np.mean(region_intensity_values))
                        median_int.append(np.median(region_intensity_values))
                        ch.append(channel)
                        im_name.append(img_name)
        # if img_file == 405 then run the following code
        # if '405' in img_file:
        #     print(img_file)
            # img = czifile.imread(img_file[0])
            # print(img.shape)
            # reshape the image to the correct dimensions (time, x, y)
            # img_ = img[0,0,0,:,0,:,:,:]
            # img_.shape

            # plt.imshow(img_[0,0,:,:])

######################## TODO:
# 1. finish up the above loop so that when channel is 405 segment and quantify the mask mean intentisites
# 2. Fiinish up the normalization of the plots below

#%%

# %%





# for i in range(img_.shape[0]):  # Assuming you're iterating over the first image for now; replace 1 with img_.shape[0] for all images
#     intensity_image = img_[i, :, :]
    
#     # It's important to pass intensity_image to regionprops to calculate intensity metrics
#     for region in regionprops(label_image, intensity_image=intensity_image):
#         # Extract the region's intensity values using its mask
#         region_intensity_values = intensity_image[region.coords[:, 0], region.coords[:, 1]]
        
#         # Append the data to the lists
#         time_point.append(i)
#         mask.append(region.label)
#         avg_int.append(np.mean(region_intensity_values))
#         median_int.append(np.median(region_intensity_values))

# %%
# normalize the intensity values min max for each mask
df_int_time = pd.DataFrame({
              'image': im_name,
              'condition': condition,
              'time_point': time_point, 
              'channel': ch,
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
    ax.plot(group.time_point, group.median_int, label=mask)
    ax.legend(title = 'Mask')
plt.show()
# %%

