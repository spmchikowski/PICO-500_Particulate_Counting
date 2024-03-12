function particulate_counting_using_matlab_V_14

%% Modified on 29 March 2018. Mostly comments
%% Usman Chowdhury- at Canadian Nuclear laboratories from April 2018
%% Contact at usmanphy@gmail.com if any help is required or if you have any questions

%% Last Modified on 15th June 2021. Comments, code changed for efficiency, and change of constants
%% Quinn Malin - University of Alberta
%% Contact at qmalin@ualberta.ca for questions


%% Provide all the required information here

%%sampleName: This is very important. Make sure that the
% image names inside this folder are in the format XY-DDMMYY-S-N, where S is the sample 
% number for a particular day. N is the image number from that particular
% sample. Usually we take 10 random pictures, all have the same size,
% zooming(i.e. applification) but may have different focus (as some samples 
% might have bumpy regions)

numberOfImage=45% Number of unbiased images of the filter taken randomly
% that you want to analyze. You can give any number here, but we usually 
% use 10 images, as mentioned above.

FirstImageNumber=2 % This is normally "1". However, sometimes the first image
% is biased when we want to focus the microscope. So you may want to skip
% the first image and start from 2, 3, 4, ... whatever you want.


e_threshold=0.60;% Should be between 0.5 and 0.64>> 0.58 worked fine. 0.53 
% is good for old halogen light
%For PICO-500, 0.60 worked well
%This constant is used to differentiate between particulates based on light
%levels. Light levels should stay the same for each photo taken



% 0.25 for 25 mm filters used with dishwasher works good
filterSize=25; % filter size in millemeter. Effective radius is 11.5 mm. Rim of the filter is ecluded
% The previous code (by Carsten and Pitam) considers only the area covered
% by the retainer metal mesh. I don't think that is correct as the water
% flow is not constrained within the holes. Contact Usman Chowdhury for
% more details.

sampleFluidVolme=1; % the amount of water passed through the filter in litre. 

%%camera specs:
    %pixel side length = 2.4; %um
    %magnification = 4.5; questionable because the calibartion slide showed otherwise
    %bin size = 1 for largest resolution
    %pix per um = magnification/pixels; %pixels/um
    %resolution = 5440x3648
  
pix_per_um = 3.74666667; %pixels/um found using calibration slide

area_pixels = 3400*2280; %pixels^2, resolution was changed in ROI of program for focusing purposes, I recommend doing the same

area_picture = area_pixels/(pix_per_um^2); %um^2
area_filter = pi*((filterSize/2)*10^3)^2; %um^2
image_per_filter = (numberOfImage-FirstImageNumber+1)*area_picture/area_filter; % multiplication factor from # of pictures
FilterAreaCorrection = 1/image_per_filter

callibrationFactorPixels2microMeter=pix_per_um
% 1.875 pixels per micrometer written on microscope, calibration slide
% showed 3.74666667
% Up to date callibration is achived by using the callibration slide


tic % counts the total time spent on processing the images. "tic" is the beginning of time count.

% 2.4 um per pixel
% so 2.4/4.5 = 0.5333, then 1/0.5333 = 1.875 pixels/um at 4.5x magnification, then the
% resolution is 2736 x 1824, so (2736/1.875)(1824/1.875) = 1419509.76
% take 1419509.76/415475628.43725 = 0.0034165897,
% times 10 (for 10 pictures) and invert: 29.26895187


%% this is another section

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Check if system has image processign toolbox %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Check of image processign toolbox finished %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% end of this section

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%   This setion reads the file     %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 [fileName, pathName] = uigetfile('*.jpg;*.tif;*.png;*.gif',...
     'Select the first image of the sample');
  ImageAddress = fullfile(pathName, fileName);
  sampleName = fileName(1:length(fileName)-6);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Find the file section ends   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% The images we take are usually very large, and Matlab sends you a friendly warning
% this unnecessarily fills the command window. So let's get rid of this
warning('off', 'Images:initSize:adjustingMag');


%% 1st section. Prepare the image for analysis


for i=FirstImageNumber:numberOfImage
    
    clf
    
    fprintf('Image # %d being processed now',i)
    
    fileName=[sampleName, '-',num2str(i), '.jpg']
    folder='sampleName';
    myImageOriginal=double(imread([pathName fileName]))/255.0;
    
    myImage=rgb2gray(myImageOriginal);
    
    myImage=imsharpen(myImage,'Radius',10,'Amount',2, 'Threshold', 0.01); % This is critical.
    % If we make the iamge too sharp, it takes for ever to analyze and also can merge
    % two nearby particulates. Play around a bit if you need to adjust any
    % of the parameters by reading the Matlab website.
    
    figure(1);
    imshow(myImage);
  
 
    size_of_image=size(myImage); % Just to display the size of the image in pixels
       
    [background] = imopen(myImage,strel('disk',5))>e_threshold; % I don't see a very strong effect of n on ('disk',n), however, this Structure element
    % operation helps define the backgroudn better. You may want to play
    % around a bit with the size of the element if necessary.
    
    figure(2); % This figure will show how the background looks.
    imshow(background);
    title('The background')
    backgroundOnly=['background' num2str(i) '.jpg']
    
    size_of_background=size(background); % Just a check of the image size again
    
    % figure(23)
    % surf(double(background(1:8:end,1:8:end))),zlim([0 255]);
    % set(gca,'ydir','reverse');
    
    myImage = myImage - background; % subtract the background from the original image
    figure(3) % Image after background subraction.
    imshow(myImage);
    title('Background Removed')
    backgroundRemoved = ['backgroundRemoved' num2str(i) '.jpg']
    imwrite(myImage,backgroundRemoved)
    
    % myImage{i} = adapthisteq(myImage{i});
    % figure(3)
    % imshow(myImage{i});
    % title('Enhanced');
    % adapthisteqApplied = ['adapthisteqApplied' num2str(i) '.jpg']
    % imwrite(myImage{i},adapthisteqApplied)
    
    % myimage2=rgb2gray(myimage2);
    % figure(5)
    % imshow(myimage2)
    
    % myimage3 = imadjust(myimage2);
    % figure(5)
    % imshow(myimage3)
    
    
    myImage = imbinarize(myImage); % imbinarize causing problems
    figure(4)
    imshow(myImage)
    
    binaryImage = ['binary' num2str(i) '.jpg']
    imwrite(myImage,binaryImage)
    
    
    se = strel('disk',2); % This operation compensates for the size of the particulates that got
    % reduced because of the bacground reduction. This operation might be bad for some applications, but for now I am ...
%    marging very closely located particulates. 10 is bad, changed to 2.
%    The higher the value of you put the further apart pixel groups will be
%    connected. Choose it carefully.

    myImage = imclose(myImage,se); % The gap between the very close pixel group is closed using the above mentioned
% rule
    
%     figure(666)
%     imshow(myImage)
%     
    
    myImage = imfill(myImage,'holes'); % During backgroudn cancellation some particulates will 
    % show up with holes in the middle. This operation fills those holes.
    
    figure(5);
    imshow(myImage);
    
    myImage=bwareaopen(myImage, 5); % remove small objects
    %myimage3 = imadjust(myimage2);
     figure(6)
    imshow(myImage);
    tinyPixelsRemoved = ['tinyPixelsRemoved' num2str(i) '.jpg']
    imwrite(myImage,tinyPixelsRemoved)
    
    
    
    %% end 1st section
    
    % Binary image is ready after this. DON'T touch anything before this
    % point
    % for now
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    showBoundaries=bwboundaries(myImage);
    
    figure(7);
    imshow(myImage);
    hold on;
    visboundaries(showBoundaries);
    boundaryMarked = ['boundaryMarked' num2str(i) '.jpg'];
%    imwrite(myImage,boundaryMarked)
    
    
    stats = regionprops('table',myImage,'Centroid',...
        'MajorAxisLength','MinorAxisLength');
    
    
    centers = stats.Centroid;
    %diameters = mean([stats.MajorAxisLength stats.MinorAxisLength],2);
    % The original code plots a circle by taking the averate of the Major and
    % Minor axis. I define the major axis length as the length of the object,
    % wich is slightly bigger than what it is in reality. But this is a very
    % good conservative estimation
    
    
    % measurementsBoxing = regionprops(myimage4, 'BoundingBox');
    %
    % xOfEachBox=measurementsBoxing.BoundingBox(:,1);
    % yOfEachBox=measurementsBoxing.BoundingBox(:,2);
    % widths = measurementsBoxing.BoundingBox(:,3);
    % heights = measurementsBoxing.BoundingBox(:,4);
    %
    %
    % radii=sqrt((widths).^2+(heights).^2)/2;
    
    radii=stats.MajorAxisLength/2;
    
    % Now we convert the pixel unit to MKS unit.
    % Here the unit is Micro Meter, or µm
    particulateSizeInPixels=stats.MajorAxisLength;
    
    particulateSizeInMicroMeterInCell{i}...
        =particulateSizeInPixels./callibrationFactorPixels2microMeter
    
    
    figure(8)
    imshow(myImage);
    hold on
    viscircles(centers,radii);
    hold off
    
    
    
    figure(9);
    imshow(myImageOriginal);
    hold on
    viscircles(centers, radii)
    
    clear myImage myimageOrigina
    
        
end


particulateSizeInMicroMeter=cell2mat(particulateSizeInMicroMeterInCell(:))


%% Sorting the particulate in different bins with specific titles



lessThan5=particulateSizeInMicroMeter(particulateSizeInMicroMeter>=0 ...
    & particulateSizeInMicroMeter<5)

lessThan15=particulateSizeInMicroMeter(particulateSizeInMicroMeter>=5 ...
    & particulateSizeInMicroMeter<15)

lessThan25=particulateSizeInMicroMeter(particulateSizeInMicroMeter>=15 ...
    & particulateSizeInMicroMeter<25)

lessThan50=particulateSizeInMicroMeter(particulateSizeInMicroMeter>=25 ...
    & particulateSizeInMicroMeter<50)

lessThan100=particulateSizeInMicroMeter(particulateSizeInMicroMeter>=50 ...
    & particulateSizeInMicroMeter<100)

lessThan5000=particulateSizeInMicroMeter(particulateSizeInMicroMeter>=100 ...
    & particulateSizeInMicroMeter<5000)

totalLessThan5microMeter=size(lessThan5,1)*FilterAreaCorrection;
totalLessThan15microMeter=size(lessThan15,1)*FilterAreaCorrection;
totalLessThan25microMeter=size(lessThan25,1)*FilterAreaCorrection;
totalLessThan50microMeter=size(lessThan50,1)*FilterAreaCorrection;
totalLessThan100microMeter=size(lessThan100,1)*FilterAreaCorrection;
totalLessThan5000microMeter=size(lessThan5000,1)*FilterAreaCorrection;

particulatesAtEachBin=[totalLessThan5microMeter, totalLessThan15microMeter,...
    totalLessThan25microMeter, totalLessThan50microMeter,...
    totalLessThan100microMeter, totalLessThan5000microMeter]./sampleFluidVolme


%sumOfWeight=sum(particulatesAtEachBin)
error=sqrt(particulatesAtEachBin);%sumOfWeight*ones(1, 6) %

particulatesAtEachBin(particulatesAtEachBin==0)=1;

dataBins=linspace(1.5,6.5,6) % move the bin location by 0.5 so that the data point is plotted at the center. Quick fix, google to find better eay

xerror=ones(size(dataBins,2),1)/2;





%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%  AREA reserved for IETS-STD-CC1246D standard histogram calculations STARTS %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

binNumber4AllLevels=[1 2 3 4 5 6 7];

particulateCount4Level25=[758 33 10 0 0 0 0];


particulateCount4Level50=[3990 246 72 10 0 0 0];


particulateCount4Level100=[40710 2640 784 107 10 0 0];


particulateCount4Level25(particulateCount4Level25==0)=1;
particulateCount4Level50(particulateCount4Level50==0)=1;
particulateCount4Level100(particulateCount4Level100==0)=1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%  AREA reserved for IETS-STD-CC1246D standard histogram ends %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%     This is the histogram DRAWING section   %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(10)

stairs(binNumber4AllLevels, particulateCount4Level25,'LineWidth',2,'color','green');
hold on;
stairs(binNumber4AllLevels, particulateCount4Level50,'LineWidth',2,'color','blue');
hold on;
stairs(binNumber4AllLevels, particulateCount4Level100,'LineWidth',2,'color','red');
hold on;
errorbar(dataBins, particulatesAtEachBin, error/2, error/2, xerror, xerror, ...
    'black.','MarkerSize',30,'LineWidth',2.5, 'CapSize',3);

xAxisLabel={'\bf\leq5µm' '\leq15µm' '\leq25µm','\leq50µm' '\leq100µm' '>100µm',};
yAxisLabel={'0','10','10^2','10^3','10^4','10^5', '10^6', '10^7','10^8'};
xtickPositions=1.5:1:8.5;
yTickRange=1:10^7;
set(gca,'xtick',xtickPositions,'xticklabel',xAxisLabel);
set(gca,'fontSize',12,'fontWeight','bold','lineWidth', 1.5,'TickLength',[0.03, 0.01]);

set(gca, 'YScale', 'log', 'yticklabel', yAxisLabel,'YMinorGrid', 'off')
%
xlabel('Particulate size','fontSize',18,'fontWeight','bold');
ylabel('Particulates/litre','fontSize',18,'fontWeight','bold');

%ytickformat('10^{%g}')
titleOfHistogram=['Sample ID:', sampleName]
title(titleOfHistogram)
legend('IETS-STD-CC1246D-25','IETS-STD-CC1246D-50','IETS-STD-CC1246D-100',...
    sampleName);

grid on;

hold off

toc % counts the total time spent on processing the images. "toc" is the end of time count.

'';






