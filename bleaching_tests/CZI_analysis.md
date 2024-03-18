# Analysis of fluoresence for images taken in different imaging buffers on an Elyra7
Authors: Eunice Fabian Morales, Antonios Lioutas

Updated on : 20240304

## Introduction
We want to select the optimal imaging buffer for Elyra 7 so we have imaged FISH samples with the following buffers:
1. SSCT2X (Control)
2. PCA/PCD
3. SORB70VE (Sorbitol 63.3% and Vectashield 20%)
4. Eternity buffer

The imaging protocol is 500time points of a single Z, and 4 channels (405, 488, 561, 642)
Each channel was saved in a separate file for 1, and 3. For 2 the images were saves in the same file. 

### Image analysis pipeline
1. Find images within a folder
2. Split filenames
3. Select for 405
4. Make a mask with OTSU threshold for each nucleus
5. Filter anything that is below 500pixels
6. Measure the mean intensity of channels 405, 488, 561, 642 for each nuleus 
7. Save the results in a csv file
8. Plot the results for comparison to choose the best buffer for imaging


#### Pipeline progess
1. [x] Find images within a folder
2. [x] Split filenames
3. [x] Select for 405
4. [x] Make a mask with OTSU threshold for each nucleus
5. [x] Filter anything that is below 500pixels
6. [x] Measure the mean intensity of channels 405, 488, 561, 642 for each nuleus
7. [ ] Save the results in a csv file
8. [ ] Plot the results for comparison to choose the best buffer for imaging


#### TO DO 
- Normalization of data plotting (Min max per nucleus or min max of all nuclei?)
- Plot all cells individually or plot the mean of all cells?
- Make the different types of images readable for the pipeline (1. channerl are split, 2. channels are in the same file)



