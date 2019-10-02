tic
clear all
addpath('C:\Pile Higher and Deeper (PhD)\Area for Data\Analysis methods\Jason coding\main_tracking source code\');

dir = 'J:\20190923\Image stack\test\';

[dname, pathname] = uigetfile('*.tif',  'MultiSelect', 'on');

%%
 for t = 1:size(dname,2)
    warning('off','MATLAB:MKDIR:DirectoryExists')
    if iscell(dname)==1 % multiple files exist
        mkdir([dir 'ROI ' strrep(dname{t},'.tif','') '\']);
        display(['Processing video ', dname{t}])
 
        parameter1.datafile = [pathname dname{t}];
        parameter1.savefigure = 1;
        parameter1.figuredir = [dir 'ROI ' [strrep(dname{t},'.tif','')] '\'];
        parameter1.localized_tracking = 1;
        parameter1.mag = 10;
        parameter1.ig_size = 512;
        parameter1.pixel_size = 235;                                 % Diffraction limit for length calculation (nm)
        
        [images, imDataOri, mag, results, boundListN, roiListN, iListN, sbrListN]=Particle_boundaries_tracking_v2(parameter1);
        %%
        if any(results) ~= 0
        display(['Analysing video ', dname{t}])
        parameter2.analysis_all = 1;
        parameter2.savedata = 1;
        parameter2.datadir = [dir 'ROI ' [strrep(dname{t},'.tif','')] '\'];
        parameter2.pixel_size = 235;
        
        [Area, Intensity, Length, SBR] = Single_particle_analysis_v2(parameter2, images, imDataOri, mag, results, boundListN, roiListN, iListN, sbrListN);
        else
        display('No particles have been tracked.')
        end
        close all
        
    else % only one file to be analysed
        mkdir([dir 'ROI ' [strrep(dname,'.tif','')] '\']);
        display(['Processing SINGLE video ', dname])
        
        parameter1.datafile = [pathname dname];
        parameter1.savefigure = 1;
        parameter1.figuredir = [dir 'ROI ' [strrep(dname,'.tif','')] '\'];
        parameter1.localized_tracking = 0;
        parameter1.mag = 10;
        parameter1.ig_size = 512;
        parameter1.pixel_size = 350;                                 % Diffraction limit for length calculation (nm)
        
        [images, imDataOri, mag, results, boundListN, roiListN, iListN, sbrListN]=Particle_boundaries_tracking_v2(parameter1);
        
        if any(results) ~= 0
        display(['Analysing single video ', dname])
        parameter2.analysis_all = 1;
        parameter2.savedata = 1;
        parameter2.datadir = [dir 'ROI ' [strrep(dname,'.tif','')] '\'];
        parameter2.pixel_size = 350;
        
        [Area, Intensity, Length, SBR] = Single_particle_analysis_v2(parameter2, images, imDataOri, mag, results, boundListN, roiListN, iListN, sbrListN);
        else
        display('No particles have been tracked.')
        end
        close all
        break
    end
 end
display('Execution completed!')
toc