%
% NAME:
%               Particle_boundaries_tracking_v2.0
%
% PURPOSE:
%               For analyzing image videos of fluorescent spots. Identify particles in specific ROIs, define the boundries and change of the 
%               boundries among frames, connect the identified particles to trajectories in subsequent image frames, and the numbered 
%               particles are to be further analyzed by Single_particle_analysis.m
%
%               Require the script of track.m and bpass.m
%               and Matlab R2015b
% 
%
%
%               Written by Dr Aleks Ponjavic and Jason C Sang, University of Cambridge, 
%               2015-2016
%
%               Last updated on 2019/05/11
%       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [images, imDataOri, mag, results, boundListN, roiListN, iListN, sbrListN]=Particle_boundaries_tracking_v2(parameter1)

clearvars -except parameter1

% Append the SPT codes
addpath('C:\Pile Higher and Deeper (PhD)\Area for Data\Analysis methods\Jason coding\main_tracking source code\')

datafile=parameter1.datafile;
savefigure=parameter1.savefigure;
figuredir=parameter1.figuredir;
localized_tracking=parameter1.localized_tracking;
mag = parameter1.mag;
ig_size = parameter1.ig_size;
pixel_size = parameter1.pixel_size;

%% Parameter setting

%datafile = 'C:\ProgramData\IObit\Driver Booster\Download\20170710 Linda\20170710_5uM p53 seeding-1.tif';      % This is the destination for loading images

%savefigure = 1;                                                                                                  % Set 1 for saving all analyzed frames; 0 for not saving

%figuredir = 'C:\Pile Higher and Deeper (PhD)\Area for Data\Analysis 20170720_Linda p53\';                                           % Destination for figure saving

%localized_tracking = 0;                                                                                    % Set 1 for tracking specific area; 0 for entire images

%mag = 10;                                                                                                           % Magnification of images in calculation

%% Stepper parameters for localized tracking

if localized_tracking==1
    
    %stepper = 50;       % Size of ROI

    %yc = 52;                % Coordinate of X
    %xc = 330;                % Coordinate of Y
       
    %y1 = yc-stepper; % Coordinate of X, from
    %y2 = yc+stepper; % Coordinate of X, to
    %x1 = xc-stepper; % Coordinate of Y, from
    %x2 = xc+stepper; % Coordinate of Y, to
    
    y1 = 110; % Coordinate of X, from
    y2 = 370; % Coordinate of X, to
    x1 = 110; % Coordinate of Y, from
    x2 = 380; % Coordinate of Y, to
end


%% Image input and processing

% Check existence of folder
if savefigure==1
    if exist(figuredir, 'dir')==0
        warning('off','MATLAB:MKDIR:DirectoryExists')
        mkdir(figuredir);
        display('Folder not existed. Generating a new folder...');
    end
end

fname = datafile;
info = imfinfo(fname);
images = numel(info);

for i = 1:images
    % Bin Up
    imTemp(:,:) = imread(fname, i);
    imTemp(imTemp == 0) = mean2(imTemp);

    if localized_tracking==1
        imDataOri(:,:,i) = imTemp(x1:x2,y1:y2); % if using stepper
    elseif localized_tracking==0
        imDataOri(:,:,i) = imTemp(:,:); % if NOT using stepper
    end
    clear imTemp
end

clear i

list = [];
boundList = [];
iList = [];
sbrList = [];
roiList = [];

for i = 1:images
    
    display(['Loading frame No ', num2str(i), ' of ', num2str(images), '. Filtering image...'])
    
    se = strel('disk', 10);
    
    avg_background(i) = sum(sum(double(imDataOri(:,:,i))))/(512*512);
    
    imDataB = imtophat(imDataOri(:,:,i), se); % top-hat filtered image

    imDataG(:,:) = imDataOri(:,:,i) - imDataB; % background per unit

    imDataBG = imresize(imDataG, mag,'nearest'); %10x background

    imDataF = imresize(imDataB, mag,'nearest'); % 10x image

    imData = imresize(imDataOri(:,:,i), mag,'nearest'); % 10x original image

    
    imDataGauss = imgaussfilt(imData,3); % Gaussian filter for images

    imDataF(:,:)=bpass(imgaussfilt(imDataF,3),3,50); % Filter for image thresholding; usually 30-50, Syn=50, PrP=30
    
    clear imDataB imDataG imData
    
    
    
%% Identify and track spots

    display('Thresholding images...')

    BW = im2bw(imDataF, avg_background(i)/1294*0.001); % Transform to BW for boundaries; set the threshold; usually 0.005-0.015 0.0015
    BW = imfill(BW,'holes');

    clear out intensity sbr image blank background

    % parameters for peak finding
    %th = 200; % threshold; usually 200
    %sz = 30; % size of the ROI; usually 30

    %surf(imDataOri(:,:,i))
    %hold on
    figure('visible','on');
    imagesc(imDataOri(:,:,i)); %Plot
    colormap(parula)
    colorbar
    daspect([1 1 1])
    hold on %Plot

    image = imDataGauss;
    blank = imDataBG;
    
    %for j=1:length(boundaries')
               
        %boundaries{j} = boundaries{j}/10+0.5;
       % b = boundaries{j};

        %out(j,1) = mean(b(:,2));
        %out(j,2) = mean(b(:,1));
        %plot(b(:,2),b(:,1),'g'); %Plot
        %plot(out(j,1),out(j,2),'rx') %Plot
                
        %B = (b-0.5)*10;
        %roi = [roi; B];
        
        %clear b B
    %end
    
    %clear j
    
    clear imT
    
    % Calculate intensity
    [imT roiN] = bwlabel(BW);
    boundaries = bwboundaries(BW, 'noholes');
    roi = [];
    display(['Found ', num2str(roiN), ' particles.'])
     
    for j = 1:roiN
        ind = find(imT==j);
        [m n] = ind2sub(size(BW), ind);
        roi{j,1} = [m n];

        intensity(j) = sum(sum(double(image(ind)))); % sum of intensity of pixels
        background(j) = sum(sum(double(blank(ind)))); % sum of intensity of pixels 
        %sbr(j) = sum(sum(double(image(ind))./double(blank(ind)))); % sum of signal-to-background ratio of each pixels
        sbr(j) = sum(sum((double(image(ind))-double(blank(ind)))./double(blank(ind))));
        roiOri{j,1} = roi{j}/10+0.5;
        
        b = boundaries{j,1}/10+0.5;
        
        out(j,1) = mean(roiOri{j,1}(:,1));
        out(j,2) = mean(roiOri{j,1}(:,2));
        
        plot(b(:,2),b(:,1),'r'); %Plot
        %plot(out(j,2),out(j,1),'rx') %Plot
        
        clear m n b ind
    end
    
    
    %pause(0.01); % for plot; time gap between frames

    if(~isempty(roi))
        out = [out ones(size(out,1),1)*i];

        list = [list; out];
        boundList = [boundList; boundaries];
        roiList = [roiList; roi];
        iList = [iList; intensity'];
        sbrList = [sbrList; sbr'];
     end

    if savefigure==1
        saveas(gcf, [figuredir 'Figure ', num2str(i), '.png']);
    end
    clear imDataBG imDataG imDataF imDataGauss BW roi boundaries image blank out intensity sbr
    
end

%% Tracking particles

param.mem = 3; % steps particle can be lost for             4
param.dim = 2; % dimensions
param.good = 0; % delete short tracks shorter than n     5
param.quiet = 1; % set this keyword to 1 if you don't want any text

if isempty(list) == 0
    if (max(list(:,3))-min(list(:,3))) ~= 0
        results = track(list(:,1:3),4,param);
    else
        results = [0 0 0 0];
    end
else
    results = [0 0 0 0];
end

figure %Plot
imagesc(imDataOri(:,:,end)) %Plot
hold on %Plot

if any(results) ~= 0
for k = 1:max(results(:,4))
    frames = find(results(:,4)==k);
    if isempty(frames)==0
    plot(results(frames,2),results(frames,1),'r-') %Plot
    colormap(parula)
    colorbar
    text(results(frames(1),2)+1,results(frames(1),1)+1,['\color{red} ' num2str(k)], 'FontSize', 5) %Plot
    daspect([1 1 1])
    end
end
end

clear k


if savefigure==1
    saveas(gcf, [figuredir 'Particles with numbering', '.png']);
    if any(results)~= 0
    display(['Tracking completed! Identified ', num2str(max(results(:,4))), ' particles.']);
    else
    display(['Tracking completed! Identified 0 particles.']);
    end
    display(['Figures saved in: "',  num2str(figuredir), '".']);
    
elseif savefigure==0
    if any(results)~= 0
    display(['Tracking completed! Identified ', num2str(max(results(:,4))), ' particles.']);
    else
    display(['Tracking completed! Identified 0 particles.']);
    end
    display('No figures saved.');
    close all
end

if any(results) ~= 0
% Organise the lists
for k = 1:size(results,1)
    xx = list(:,1)==results(k,1);
    yy = list(:,2)==results(k,2);
    pos = xx&yy;
    
    posVal = find(pos==1);
    boundListN(k) = boundList(posVal(1));
    roiListN(k) = roiList(posVal(1));
    iListN(k) = iList(posVal(1));
    sbrListN(k) = sbrList(posVal(1));
end

else
    boundListN = [];
    roiListN = [];
    iListN = [];
    sbrListN = [];
end
end