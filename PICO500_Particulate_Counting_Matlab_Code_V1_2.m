function PICO500_Particulate_Counting_Matlab_Code_V1_2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% PICO-500 Particulate Counting Matlab Code Version 1.2
%%%%%%%%%% Last Updated: March 5th, 2024
%%%%%%%%%% Updated by: Shawn Miller-Chikowski
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Version 1.0
%% Modified on 29 March 2018. Mostly comments
%% Usman Chowdhury- at Canadian Nuclear laboratories from April 2018
%% Contact at usmanphy@gmail.com if any help is required or if you have any questions

%% Version 1.1
%% Last Modified on 15th June 2021. Comments, code changed for efficiency, and change of constants
%% Quinn Malin - University of Alberta
%% Contact at qmalin@ualberta.ca for questions

%% Version 1.2
%% Last Modified in March 2024. 
%% Additional Comments, code changed for efficiency -> only the first iteration shows images if user requested, added user inputs (vlose all plots, show all images, etc), added alternative plot display, code formatting for readability 
%% The main fix was to the error calculation. Error is sqrt(n_bins) * FilterAreaCorrectionFactor * Z_90CI -> error bars are of the 90% Confidence interval 
%% sampleID will need to be adjusted if different sequencing formats are used. Naming consistency will allow this to be a constant value
%% Shaw Miller-Chikowski - University of Alberta
%% Contact Quinn Malin if there are questions: qmalin@ualberta.ca

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Provide all microscope and image related constants here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%sampleName: This is very important. Make sure that the
% image names inside this folder are in the format XY-DDMMYY-S-N, where S is the sample 
% number for a particular day. N is the image number from that particular
% sample. Usually we take 10 random pictures, all have the same size,
% zooming(i.e. applification) but may have different focus (as some samples 
% might have bumpy regions)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FirstImageNumber= 1 ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is normally "1". However, sometimes the first image
% is biased when we want to focus the microscope. So you may want to skip
% the first image and start from 2, 3, 4, ... whatever you want.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
e_threshold = 0.42 ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Should be between 0.5 and 0.64>> 0.58 worked fine. 0.53 
% is good for old halogen light
%For PICO-500, 0.64 worked well
%This constant is used to differentiate between particulates based on light
%levels. Light levels should stay the same for each photo taken
%This constant is very important for differentiating between particulates
%and not particulates

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filterSize = 25 ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 0.25 for 25 mm filters used with dishwasher works good
% filter size in millemeter. Effective radius is 11.5 mm. Rim of the filter is excluded
% The previous code (by Carsten and Pitam) considers only the area covered
% by the retainer metal mesh. I don't think that is correct as the water
% flow is not constrained within the holes
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sampleFluidVolme = 1 ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the amount of water passed through the filter in litre. 

%%camera specs:
    %pixel side length = 2.4; %um
    %magnification = 9; 4.5x plus 2x barlow lense
    %bin size = 1 for largest resolution
    %pix per um = magnification/pixels; %pixels/um
    %resolution = 5440x3648

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pix_per_um = 4.082641 ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%pixels/um found using calibration slide
%you can also calculate it as follows:
% (magnification/pixel side length) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
area_pixels = 5440*3648 ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%pixels^2, resolution was changed in ROI of program for focusing purposes, I recommend doing the same if using the same microscope

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% counts the total time spent on processing the images. "tic" is the beginning of time count.
tic

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Check if system has image processign toolbox %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check that user has the Image Processing Toolbox installed.
hasIPT = license('test', 'image_toolbox');
if ~hasIPT
    % User does not have the toolbox installed.
    message = sprintf('Sorry, but you do not seem to have the Image Processing Toolbox.\nDo you want to try to continue anyway?');
    reply = questdlg(message, 'Toolbox missing', 'Yes', 'No', 'Yes');
    if strcmpi(reply, 'No')
        % User said No, so exit.
        return;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% end of this section

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%   This setion reads the file     %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% User selects the first image
[fileName, pathName] = uigetfile('*.jpg;*.jpeg;*.tif;*.png;*.gif', 'Select the first image of the sample');
if isequal(fileName,0)
   disp('User selected Cancel');
   return;
end

disp(' ');
% Ask the user for the number of images to process or 0 to process all sequentially
numImagesToProcess = input('Enter the number of images to analyze (or 0 to analyze all): ');

% Extract the base name (excluding extension) and sequential number from the filename
[~, baseName, ext] = fileparts(fileName);
numericPartMatch = regexp(baseName, '-(\d+)$', 'tokens');
if isempty(numericPartMatch)
    disp('No numeric part found in the filename.');
    return;
else
    numericPart = str2double(numericPartMatch{1}{1});
end

% Prepare the prefix for subsequent filenames, which excludes the numeric part
prefix = baseName(1:end-length(numericPartMatch{1}{1})-1);

% Initialize the counter with the extracted numeric part
counter = numericPart;

% Initialize a cell array to hold the file names
validFileNames = {};

% Determine the initial number of digits in the numeric part
initialNumDigits = length(numericPartMatch{1}{1});

% Start the loop from the extracted numeric part
while true
    % Determine the current number of digits needed for the counter
    currentNumDigits = max(initialNumDigits, length(num2str(counter)));
    
    % Construct the filename using dynamic formatting to maintain initial leading zeros
    formatStr = sprintf('%%s-%%0%dd%%s', currentNumDigits);
    currentFileName = sprintf(formatStr, prefix, counter, ext);
    currentFilePath = fullfile(pathName, currentFileName);
    
    % Check if the file exists
    if exist(currentFilePath, 'file')
        % File exists, add it to the list
        validFileNames{end+1} = currentFileName;
        
        % Increment the counter for the next file
        counter = counter + 1;
        
        % If we've reached the user-defined limit, if any, exit the loop
        if numImagesToProcess > 0 && length(validFileNames) >= numImagesToProcess
            break; % Exit the loop if the limit is reached
        end
    else
        % No more files found, exit the loop
        disp(' ');
        disp(['No more files after ', currentFileName]);
        break;
    end
end

% Now, process each file collected in validFileNames
%for i = 1:length(validFileNames)
%    currentFilePath = fullfile(pathName, validFileNames{i});
%    img = imread(currentFilePath);
%    disp(['Processed ', validFileNames{i}]);
%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sampleID = fileName(1:length(fileName) - 9  );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The images we take are usually very large, and Matlab sends you a friendly warning
% this unnecessarily fills the command window. So let's get rid of this
warning('off', 'Images:initSize:adjustingMag');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Prepare the image for analysis

disp(' ');
userChoice = input('Show All images? If no, only the first image is processed for demonstration purposes (y/n): ', 's');

% Ensure userResponse is processed correctly
if isempty(userChoice) || (~strcmpi(userChoice, 'y') && ~strcmpi(userChoice, 'n'))
    error('Invalid input. Please enter Y or N.');
end

showAllImages = lower(userChoice) == 'y';

% Initialize a variable to track whether it's the first iteration
isFirstIteration = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate constants based on number of images and initial defined parameters

numberOfImage = length(validFileNames); % Number of unbiased images of the filter taken randomly

area_picture = area_pixels/(pix_per_um^2); %um^2
area_filter = pi*((filterSize/2)*10^3)^2; %um^2
image_per_filter = (numberOfImage-FirstImageNumber+1)*area_picture/area_filter; % multiplication factor from # of pictures
FilterAreaCorrection = 1/image_per_filter;

% Print out the Percent of the fitler surveyed
AreaCorrectionPercent = 100 * 1/FilterAreaCorrection;

callibrationFactorPixels2microMeter=pix_per_um;
% 3.75 pixels per micrometer written on microscope, calibration slide
% showed 3.74666667
% Up to date callibration is achived by using the callibration slide

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Loop sequentially from the first image picked to the number of images the
% user input
for i = 1:length(validFileNames)
    
    clf; % Clear the current figure
    
    disp(' ');
    fprintf('Image #%d being processed now, ', i);
    
    currentFileName = validFileNames{i}; % Correctly reference the current file
    
    % Corrected file reading using 'currentFileName'
    myImageOriginal = double(imread(fullfile(pathName, currentFileName)))/255.0;
    %myImage = rgb2gray(myImageOriginal);
    myImage = im2gray(myImageOriginal);
    myImage = imsharpen(myImage, 'Radius', 10, 'Amount', 2, 'Threshold', 0.01);
    
    % Display processed grayscale image based on user choice
    if showAllImages || isFirstIteration
        figure; imshow(myImage); title(['Processed Image #', num2str(i)]);
    end

    [background] = imopen(myImage,strel('disk',5)) > e_threshold; 
    if showAllImages || isFirstIteration
        figure; imshow(background); title('The background'); % Show background extraction
    end
    
    myImage = myImage - background; 
    if showAllImages || isFirstIteration
        figure; imshow(myImage); title('Background Removed'); % Show image after background subtraction
    end
    
    myImage = imbinarize(myImage); 
    if showAllImages || isFirstIteration
        figure; imshow(myImage); title('Binarized Image'); % Show binarized image
    end

    se = strel('disk',2); 
    myImage = imclose(myImage,se); 
    
    myImage = imfill(myImage,'holes'); 
    if showAllImages || isFirstIteration
        figure; imshow(myImage); title('Filled Holes'); % Show image after filling holes
    end
    
    myImage = bwareaopen(myImage, 5); 
    if showAllImages || isFirstIteration
        figure; imshow(myImage); title('Removed Small Objects'); % Show image after removing small objects
    end

    showBoundaries = bwboundaries(myImage);
    if showAllImages || isFirstIteration
        figure; imshow(myImage); hold on; visboundaries(showBoundaries); title('Boundaries'); % Show boundaries
    end

    stats = regionprops('table',myImage,'Centroid','MajorAxisLength','MinorAxisLength');
    centers = stats.Centroid;
    radii = stats.MajorAxisLength / 2;

    % Now we convert the pixel unit to MKS unit.
    % Here the unit is Micro Meter, or µm
    particulateSizeInPixels=stats.MajorAxisLength;
    
    particulateSizeInMicroMeterInCell{i} = particulateSizeInPixels./callibrationFactorPixels2microMeter;
    
    if showAllImages || isFirstIteration
        figure; imshow(myImage); hold on; viscircles(centers, radii); title('Detected Particles'); % Show detected particles
    end
    
    if showAllImages || isFirstIteration
        figure; imshow(myImageOriginal); hold on; viscircles(centers, radii); title('Original Image with Detected Particles'); % Show original image with detected particles
    end
    
    % After the first iteration, if 'n' was selected, ensure no more images are displayed
    isFirstIteration = false;
end

particulateSizeInMicroMeter=cell2mat(particulateSizeInMicroMeterInCell(:));
numItems = length(particulateSizeInMicroMeter);
%disp(['Number of items: ' num2str(numItems)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sorting the particulate in different bins with specific titles

lessThan5=particulateSizeInMicroMeter(particulateSizeInMicroMeter>=0 ...
    & particulateSizeInMicroMeter<5);
% numItems_lessThan5 = length(lessThan5);
% mean_lessThan(1) = numItems_lessThan5/45
% % Calculate the sum of squared differences
% sumOfSquaredDifferences_lessThan5 = sum((lessThan5 - mean_lessThan(1)).^2);
% % Calculate the standard deviation
% std_lessThan5 = sqrt(sumOfSquaredDifferences_lessThan5 / 45);

lessThan15=particulateSizeInMicroMeter(particulateSizeInMicroMeter>=5 ...
    & particulateSizeInMicroMeter<15);
% numItems_lessThan15 = length(lessThan15);
% mean_lessThan(2) = numItems_lessThan15/45
% % Calculate the sum of squared differences
% sumOfSquaredDifferences_lessThan15 = sum((lessThan15 - mean_lessThan(2)).^2);
% % Calculate the standard deviation
% std_lessThan15 = sqrt(sumOfSquaredDifferences_lessThan15 / 45);

lessThan25=particulateSizeInMicroMeter(particulateSizeInMicroMeter>=15 ...
    & particulateSizeInMicroMeter<25);
% numItems_lessThan25 = length(lessThan25);
% mean_lessThan(3) = numItems_lessThan25/45
% % Calculate the sum of squared differences
% sumOfSquaredDifferences_lessThan25 = sum((lessThan25 - mean_lessThan(3)).^2);
% % Calculate the standard deviation
% std_lessThan25 = sqrt(sumOfSquaredDifferences_lessThan25 / 45);

lessThan50=particulateSizeInMicroMeter(particulateSizeInMicroMeter>=25 ...
    & particulateSizeInMicroMeter<50);
% numItems_lessThan50 = length(lessThan50);
% mean_lessThan(4) = numItems_lessThan50/45
% % Calculate the sum of squared differences
% sumOfSquaredDifferences_lessThan50 = sum((lessThan50 - mean_lessThan(4)).^2);
% % Calculate the standard deviation
% std_lessThan50 = sqrt(sumOfSquaredDifferences_lessThan50 / 45);

lessThan100=particulateSizeInMicroMeter(particulateSizeInMicroMeter>=50 ...
    & particulateSizeInMicroMeter<100);
% numItems_lessThan100 = length(lessThan100);
% mean_lessThan(5) = numItems_lessThan100/45
% % Calculate the sum of squared differences
% sumOfSquaredDifferences_lessThan100 = sum((lessThan100 - mean_lessThan(5)).^2);
% % Calculate the standard deviation
% std_lessThan100 = sqrt(sumOfSquaredDifferences_lessThan100 / 45);

lessThan5000=particulateSizeInMicroMeter(particulateSizeInMicroMeter>=100 ...
    & particulateSizeInMicroMeter<5000);

% Find the number of particulates sampled. Needed to get the error
totalLessThan5microMeter=size(lessThan5,1);
totalLessThan15microMeter=size(lessThan15,1);
totalLessThan25microMeter=size(lessThan25,1);
totalLessThan50microMeter=size(lessThan50,1);
totalLessThan100microMeter=size(lessThan100,1);
totalLessThan5000microMeter=size(lessThan5000,1);

particulatesAtEachBin_raw=[totalLessThan5microMeter, totalLessThan15microMeter,...
    totalLessThan25microMeter, totalLessThan50microMeter,...
    totalLessThan100microMeter, totalLessThan5000microMeter]./sampleFluidVolme;

% 1 sigma error before area correction
error=sqrt(particulatesAtEachBin_raw);

% Find the total number of particulates by applying the filter area
% correction facotr
totalLessThan5microMeter=size(lessThan5,1)*FilterAreaCorrection;
totalLessThan15microMeter=size(lessThan15,1)*FilterAreaCorrection;
totalLessThan25microMeter=size(lessThan25,1)*FilterAreaCorrection;
totalLessThan50microMeter=size(lessThan50,1)*FilterAreaCorrection;
totalLessThan100microMeter=size(lessThan100,1)*FilterAreaCorrection;
totalLessThan5000microMeter=size(lessThan5000,1)*FilterAreaCorrection;

particulatesAtEachBin=[totalLessThan5microMeter, totalLessThan15microMeter,...
    totalLessThan25microMeter, totalLessThan50microMeter,...
    totalLessThan100microMeter, totalLessThan5000microMeter]./sampleFluidVolme;

% Convert and area scale error to the two-sided 90% CI
% For two-sided 90% CI, the z-value is 1.645 for a normal distribution;
% Counts are high, so ~ normal. 
z_value = 1.645;
CI_90_error = error * z_value * FilterAreaCorrection;

dataBins=linspace(1.5,6.5,6); % move the bin location by 0.5 so that the data point is plotted at the center. Quick fix, google to find better eay

% Used to draw horizontal lines
xerror = ones(size(dataBins, 2), 1) / 2; 

% Print out some useful information
Total_Particulates_Area_Corrected = round(sum(particulatesAtEachBin));
Total_Particulates_Area_Corrected_error = round(sqrt(sum(particulatesAtEachBin_raw)) * FilterAreaCorrection * z_value);

disp(' ');
disp(' ');
disp(['Sample ID: ', string(sampleID)]);

disp(['e-Threshold Used:                                     ', num2str(e_threshold)]);
disp(['Number of Images:                                     ', num2str(numberOfImage)]);
disp(['Total Particulates (Area Corrected):                  ', num2str(Total_Particulates_Area_Corrected)]);
disp(['Total Particulates (Area Corrected) Error 90% CI:     ', num2str(Total_Particulates_Area_Corrected_error)]);
disp(['Filter Area Surveyed (%):                             ', num2str(AreaCorrectionPercent)]);
disp(['Filter Area Correction:                               ', num2str(FilterAreaCorrection)]);
disp(' ');

% Convert array to a formatted string, specifying number of decimal places
particulatesAtEachBin_raw_2string = sprintf('%.2f ', particulatesAtEachBin_raw); % Creates a string of array elements formatted to 2 decimal places
particulatesAtEachBin_raw_error_2string = sprintf('%.2f ', error * z_value); % Creates a string of array elements formatted to 2 decimal places

particulatesAtEachBin_2string = sprintf('%.2f ', particulatesAtEachBin); % Creates a string of array elements formatted to 2 decimal places
particulatesAtEachBin_error_2string = sprintf('%.2f ', CI_90_error); % Creates a string of array elements formatted to 2 decimal places

% Raw Counts
disp(['Total Binned Particulates:                            ', particulatesAtEachBin_raw_2string]);
disp(['Total Binned Particulates Error (90% CI):             ', particulatesAtEachBin_raw_error_2string]);
disp(' ');

% Corrected Counts
disp(['Total Binned Particulates (Corrected):                ', particulatesAtEachBin_2string]);
disp(['Total Binned Particulates (Corrected) Error (90% CI): ', particulatesAtEachBin_error_2string]);
disp(' ');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%  AREA reserved for IEST-STD-CC1246D standard histogram calculations STARTS %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

binNumber4AllLevels=[1 2 3 4 5 6 7];
particulateCount4Level25=[758 33 10 0 0 0 0];
particulateCount4Level50=[3990 246 72 10 0 0 0];
particulateCount4Level100=[40710 2640 784 107 10 0 0];
particulateCount4Level25(particulateCount4Level25==0)=1;
particulateCount4Level50(particulateCount4Level50==0)=1;
particulateCount4Level100(particulateCount4Level100==0)=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%     This is the histogram DRAWING section   %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(10)

% Draw the histograms of the standards
stairs(binNumber4AllLevels, particulateCount4Level100,'LineWidth', 2.75,'color',[0.6,0, 0]);
hold on;

stairs(binNumber4AllLevels, particulateCount4Level50,'LineWidth', 2.75,'color',[0,0,0.6]);
hold on;

stairs(binNumber4AllLevels, particulateCount4Level25,'LineWidth', 2.75,'color',[0,0.6,0]);
hold on;

 errorbar(9, particulatesAtEachBin(1), CI_90_error(5)/2, CI_90_error(5)/2, xerror, xerror, 'LineStyle', 'none', ...
    'Color', 'black', 'Marker', '.', 'MarkerSize', 35, 'LineWidth', 3, 'CapSize', 12); % Dummy Plot. For legend symbol only

% Plot vertical Y-axis error bars with caps
hError = errorbar(dataBins, particulatesAtEachBin, CI_90_error/2, CI_90_error/2, 'LineStyle', 'none', ...
    'Color', 'black', 'Marker', '.', 'MarkerSize', 35, 'LineWidth', 3, 'CapSize', 12);

hold on;

% Manually draw horizontal lines for X-axis errors without caps
for i = 1:length(dataBins)
    % Draw horizontal line for each data point
    line([dataBins(i) - xerror(i), dataBins(i) + xerror(i)], [particulatesAtEachBin(i), particulatesAtEachBin(i)], ...
         'Color', 'black', 'LineWidth', 3);
end

% Do some plot labelling
xAxisLabel={'\bf\leq5µm' '\leq15µm' '\leq25µm','\leq50µm' '\leq100µm' '>100µm',};
yAxisLabel={'0','10','10^2','10^3','10^4','10^5', '10^6', '10^7','10^8'};
xtickPositions=1.5:1:8.5;
yTickRange=1:10^7;
xlim([1,7]);
set(gca,'xtick',xtickPositions,'xticklabel',xAxisLabel);
set(gca,'fontSize',20,'fontWeight','bold','lineWidth', 1.5,'TickLength',[0.03, 0.01]);

set(gca, 'YScale', 'log', 'yticklabel', yAxisLabel,'YMinorGrid', 'off')
%
xlabel('Particulate size','fontSize',24,'fontWeight','bold');
ylabel('Particulates / litre','fontSize',24,'fontWeight','bold');

%ytickformat('10^{%g}')
titleOfHistogram=['Sample ID: ', sampleID];
title(titleOfHistogram)
legend('IEST-STD-CC1246D-100','IEST-STD-CC1246D-50','IEST-STD-CC1246D-25',...
    sampleID, 'Fontsize', 20, 'Location', 'best');

% Enable the grid
grid on;

% Get the current axes handle
ax = gca;

% Set major grid lines style
ax.GridLineStyle = '-';  % Solid line
ax.GridAlpha = 0.7;     % Transparency
ax.GridColor = [0.8, 0.8, 0.8];  % Light grey color

% Enable minor grid lines only in Y
ax.XMinorGrid = 'off';   % Turn off minor grid lines on X-axis
ax.YMinorGrid = 'on';    % Turn on minor grid lines on Y-axis

% Set minor grid lines style (same as major for consistency)
ax.MinorGridLineStyle = '-';  % Solid line
ax.MinorGridColor = [0.8, 0.8, 0.8];  % Light grey color (same as major grid lines)
ax.MinorGridAlpha = 0.1;  

hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% counts the total time spent on processing the images. "toc" is the end of time count.
toc
disp(' ');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%     Ask User if they want to close all figures   %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%prompt the user to close all figures
closeFiguresDialog();

function closeFiguresDialog
    % Get screen size
    screenSize = get(0, 'ScreenSize'); % Returns [left, bottom, width, height]
    screenWidth = screenSize(3);
    screenHeight = screenSize(4);
    
    % Define the dialog size
    dialogWidth = 300;
    dialogHeight = 100;
    
    % Calculate position for center-right alignment
    posX = screenWidth - dialogWidth - 50; % 50 pixels from the right edge of the screen
    posY = (screenHeight / 2) - (dialogHeight / 2); % Vertically centered
    
    % Create a figure for the dialog
    dlg = figure('Name', 'Close All Figures?', ...
                 'NumberTitle', 'off', ...
                 'MenuBar', 'none', ...
                 'ToolBar', 'none', ...
                 'Position', [posX, posY, dialogWidth, dialogHeight], ...
                 'CloseRequestFcn', @closeDialog); % Handle closing the window

    % Add "OK" button
    uicontrol('Style', 'pushbutton', ...
              'Position', [50 30 100 40], ...
              'String', 'OK', ...
              'Callback', @okCallback);

    % Function to close all figures and the dialog itself
    function okCallback(~, ~)
        close all; % Close all figures
        delete(dlg); % Close the dialog window
    end
    
    % Function to simply close the dialog if the window'sCan  close button
    % is clicked 
    function closeDialog(~, ~)
        delete(dlg);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end






