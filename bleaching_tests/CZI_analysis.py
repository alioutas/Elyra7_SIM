#%%
import czifile
import napari
import re
import cellpose
from cellpose import models, io

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

from skimage import data
from skimage.segmentation import clear_border
from skimage.measure import regionprops
from skimage.color import label2rgb
from skimage.measure import find_contours

from scipy.optimize import curve_fit




import numpy as np
import pandas as pd

from glob import glob 
import os

from czi_analysis_utils import create_masks
from tqdm import tqdm


# %%

# find all folders in the following path
path = "/Volumes/Genetics/Wu_Lab-Vutara/Experiments/Eunice/Elyra_Eunice/WGI/bleaching_test_new"
folders = [os.path.join(path, f) for f in os.listdir(path) if os.path.isdir(os.path.join(path, f))]

# Set the model parameters
model = models.Cellpose(gpu=False, model_type='nuclei')
diameter = 200

# Create a conditions list based on the folder names
conditions = []
for folder in folders:
    if 'SSC' in folder:
        conditions.append('SSCT2X')
    elif 'Sorb' in folder:
        conditions.append('Sorb70VE')
    elif 'PCA' in folder:
        conditions.append('PCA_PCD')
    else:
        conditions.append('unknown')

#%%
#Calculate fluorescence intensity
time_point = []
mask_label = []
avg_int = []
median_int = []
sum_int = []
sd_int = []
ch = []
im_name = []
area = []
treatment = []


#%%
for folder in folders:

    img_files = glob(os.path.join(folder, '*.czi'))
    condition = conditions[folders.index(folder)]
    # from the file name select any character before the first underscore
    # this will be the name of the image
    img_names = [img_file.split('/')[-1].split('_')[0] for img_file in img_files]
    # find unique image names and order them
    img_names = sorted(list(set(img_names)))

    channels = sorted(list(set([img_file.split('/')[-1].split('_')[1] for img_file in img_files])))
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
                    # Run the segmentation with cellpose
                    label_image, flows, styles, diams = model.eval(image, diameter=diameter, channels=[0, 0])

                    # remove artifacts connected to image border
                    label_image = clear_border(label_image)

                    # plot image and masks
                    image_label_overlay = label2rgb(label_image, image=image, bg_label=0, kind='overlay')
                    fig, ax = plt.subplots(figsize=(10, 6))
                    ax.imshow(image, cmap='gray')
                    ax.imshow(image_label_overlay, alpha=0.5)
                    for region in regionprops(label_image):
                        # Get contour coordinates
                        contours = find_contours(label_image == region.label, 0.5)
                        for contour in contours:
                            ax.plot(contour[:, 1], contour[:, 0], linewidth=2, c='red')
                    ax.set_axis_off()
                    plt.tight_layout()
                    plt.show()

                    # plot the max difference in time with the mask contour in red
                    img_t_diff = img_[1:,:,:,0] - img_[:-1,:,:,0] 
                    # max the first dimension of img_t_diff and plot
                    fig, ax = plt.subplots(figsize=(10, 6))
                    im = ax.imshow(img_t_diff.max(axis=(0)), cmap='gray') 
                    plt.colorbar(im) 
                    for region in regionprops(label_image):
                        # Get contour coordinates
                        contours = find_contours(label_image == region.label, 0.5)
                        for contour in contours:
                            ax.plot(contour[:, 1], contour[:, 0], linewidth=2, c='red')
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
                            mask_label.append(region.label)
                            avg_int.append(np.mean(region_intensity_values))
                            median_int.append(np.median(region_intensity_values))
                            sum_int.append(np.sum(region_intensity_values))
                            sd_int.append(np.std(region_intensity_values))
                            ch.append(channel)
                            im_name.append(img_name)
                            area.append(region.area)
                            treatment.append(condition)
                else:
                    #use the mask created from the 405 channel to quantify the mean intensity of the other channels                # read image
                    img = czifile.imread(img_file[0])
                    img_ = img[0,0,0,:,0,:,:,0] # last dimension is channel
                    
                    
                    # plot image and masks
                    fig, ax = plt.subplots(figsize=(10, 6))
                    ax.imshow(img_[0,:,:]*10, cmap='gray')
                    # draw the masks on the image as a read contour line
                    for region in regionprops(label_image):
                        # Get contour coordinates
                        contours = find_contours(label_image == region.label, 0.5)
                        for contour in contours:
                            ax.plot(contour[:, 1], contour[:, 0], linewidth=2, c='red')
                    ax.set_axis_off()
                    plt.tight_layout()
                    plt.show()

                    # plot the max difference in time with the mask contour in red
                    img_t_diff = img_[1:,:,:] - img_[:-1,:,:] 
                    # max the first dimension of img_t_diff and plot
                    fig, ax = plt.subplots(figsize=(10, 6))
                    im = ax.imshow(img_t_diff.max(axis=(0)), cmap='gray') 
                    plt.colorbar(im) 
                    for region in regionprops(label_image):
                        # Get contour coordinates
                        contours = find_contours(label_image == region.label, 0.5)
                        for contour in contours:
                            ax.plot(contour[:, 1], contour[:, 0], linewidth=2, c='red')
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
                            mask_label.append(region.label)
                            avg_int.append(np.mean(region_intensity_values))
                            median_int.append(np.median(region_intensity_values))
                            sum_int.append(np.sum(region_intensity_values))
                            sd_int.append(np.std(region_intensity_values))
                            ch.append(channel)
                            im_name.append(img_name)
                            area.append(region.area)
                            treatment.append(condition)


# %%
# Create a df of all values
df_int_time = pd.DataFrame({
              'image': im_name,
              'condition': treatment,
              'time_point': time_point, 
              'channel': ch,
              'mask': mask_label, 
              'avg_int': avg_int, 
              'median_int': median_int,
              'sum_int': sum_int,
              'sd_int': sd_int,
              'area': area
              })

# save the df
# df_int_time.to_csv('/Volumes/Genetics/Wu_Lab-Vutara/Experiments/Eunice/Elyra_Eunice/WGI/analysis_output/20240320_df_int_time_raw.csv')


#%%
# read in the results df
df_int_time = pd.read_csv('/Volumes/Genetics/Wu_Lab-Vutara/Experiments/Eunice/Elyra_Eunice/WGI/analysis_output/20240319_df_int_time_raw.csv')

#%%



#%%
# normalize the intensity values min max for each mask
groups = ['image','mask','channel', 'condition']
# df_int_time['norm_avg_int'] = (df_int_time['avg_int'] - df_int_time.groupby('mask').avg_int.transform('min')) / (df_int_time.groupby('mask').avg_int.transform('max') - df_int_time.groupby('mask').avg_int.transform('min'))
# df_int_time['norm_median_int'] = (df_int_time['median_int'] - df_int_time.groupby('mask').median_int.transform('min')) / (df_int_time.groupby('mask').median_int.transform('max') - df_int_time.groupby('mask').median_int.transform('min'))
df_int_time['norm_avg_int'] = (df_int_time['avg_int'] - df_int_time.groupby(groups).avg_int.transform('min')) / (df_int_time.groupby(groups).avg_int.transform('max') - df_int_time.groupby(groups).avg_int.transform('min'))
df_int_time['norm_median_int'] = (df_int_time['median_int'] - df_int_time.groupby(groups).median_int.transform('min')) / (df_int_time.groupby(groups).median_int.transform('max') - df_int_time.groupby(groups).median_int.transform('min'))

# now bring the values for all conditions between 0 and 1
df_int_time['norm_avg_int'] = (df_int_time['norm_avg_int'] - df_int_time.groupby('condition').norm_avg_int.transform('min')) / (df_int_time.groupby('condition').norm_avg_int.transform('max') - df_int_time.groupby('condition').norm_avg_int.transform('min'))
df_int_time['norm_median_int'] = (df_int_time['norm_median_int'] - df_int_time.groupby('condition').norm_median_int.transform('min')) / (df_int_time.groupby('condition').norm_median_int.transform('max') - df_int_time.groupby('condition').norm_median_int.transform('min'))


#%%
# sns.kdeplot(data=df_int_time, x="avg_int", hue="condition", log_scale=True, col = 'channel')
g = sns.FacetGrid(df_int_time, col="channel", hue='condition')
# g.map(sns.lineplot, 'time_point', 'norm_avg_int') #hue='channel' #(data=df_int_time,
g.map(sns.kdeplot, "time_point", "norm_avg_int", ci=None)
g.add_legend()

# %%
g = sns.FacetGrid(df_int_time, col="channel", hue='condition')
g.map(sns.lineplot, 'time_point', 'norm_avg_int') #hue='channel' #(data=df_int_time,
# g.map(sns.regplot, "time_point", "norm_avg_int", ci=None)
g.add_legend()
g.fig.suptitle(f'Normalized average intensity for {groups}', y = 1.2)

#%%
# Calculate hald life of measured intensities for each condition and channel
half_time_df = []
I0_df = []
k_df = []
treatment_df = []
ch_df = []

for condition in df_int_time['condition'].unique():
    for channel in df_int_time['channel'].unique():
        # get the data for the condition and channel
        data = df_int_time[(df_int_time['condition'] == condition) & (df_int_time['channel'] == channel)]
        # Sample data - Replace this with your actual fluorescence intensity data
        time = data['time_point'] 
        intensity = data['norm_avg_int']

        # Define exponential decay function
        def exponential_decay(t, I0, k):
            return I0 * np.exp(-k * t)

        # Fit the exponential decay model to the data
        popt, pcov = curve_fit(exponential_decay, time, intensity)

        # Extract fitted parameters
        I0_fit, k_fit = popt
        I0_df.append(I0_fit)
        k_df.append(k_fit)

        # Calculate half-life from fitted decay rate constant (k)
        half_life = np.log(2) / k_fit
        half_time_df.append(half_life)

        treatment_df.append(condition)
        ch_df.append(channel)

half_times = pd.DataFrame({'condition': treatment_df, 
                           'channel': ch_df, 
                           'half life': half_time_df, 
                           'I0': I0_df, 
                           'k': k_df})


#%%
sns.stripplot(data = half_times, 
              x= 'channel', 
              y = 'half life', 
              hue = 'condition', 
              size = 8, 
              dodge = True, 
              log_scale=True)

#change font size
plt.xticks(fontsize=12)

# add text annotation to the max value per channel an condition
unique_chan = half_times.channel.unique()
max_half_time = half_times.groupby(['channel', 'condition']).max('half life').reset_index()
for i in range(len(unique_chan)):
    plt.text(i+0.05, max_half_time.loc[max_half_time['channel'] == unique_chan[i], 'half life'].iloc[0], 
             f"{int(max_half_time.loc[max_half_time['channel'] == unique_chan[i], 'half life'].iloc[0])}", 
             ha = 'left', va = 'bottom')


