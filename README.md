# PICO-500 Particulate Counting Matlab Code - Version Overview

This document provides an overview of the different versions of the MATLAB scripts used for particulate counting in the PICO experiments.

## Version 1.0
particulate_counting_using_matlab_V1.m

- **Description**: Original PICO-40L Particulate counting code.
- **Known Issues**:
  - The error calculation is incorrect, the microscope parameters are outdated

## Version 1.1
PICO500_Particulate_Counting_Matlab_Code_V1_1.m

- **Description**: This version introduces updated microscope calibration parameters to enhance the accuracy of measurements.
- **Improvements**:
  - Updated microscope calibration parameters.
  - General comments within the code have been updated for clarity.
  - Some small changes for code efficiency

## Version 1.2
PICO500_Particulate_Counting_Matlab_Code_V1_2.m

- **Description**: Significant structural changes have been made to improve user interaction and error calculation.
- **Major Changes**:
  - **User Interaction**: Users must now specify the number of images to analyze in a prompt. Entering `0` will command the script to automatically analyze all images in the folder of the first image selected in the previous user prompt.
  - **Output Option**: Users will be prompted to decide if all analyzed images should be outputted. This is essential for reviewing if the e-threshold is consistent across all images. Selecting 'n' will speed up computation time and only produce figures from the first loop iteration (first image analyzed).
  - **Close Figure Prompt**: A prompt now appears at the end to simultaneously close all figures.
  - **Error Correction**: The error in counts has been corrected and scaled to the 90% confidence interval. Error bar reflect this change. Small visual changes to the  plotting for better readability and organization.
  - **Code Organization and Cleanup**: General code cleanup and additional comments have been added for better readability and understanding of the new changes. The script also organizes the output in a more user-friendly format.
  - The name of Version 1.1 and Version 1.2 was changed to reflect its application towards PICO-500.

### Additional Notes

-One may want to review the area correction factor calculation. This and the calibration parameters where assumed correct and up-to-date during the changes in Version 1.2.
-The SampleID variable is hard coded. Due to inconsistency in the naming of samples in 2023, this was left as is. Since a consistent system will be in place for the 
real cleaning of PICO-500, this is only expected to be adjusted once (if need be) to display the sample ID properly on the histograms. SampleID is defined on line 195 in Version 1.2. 
