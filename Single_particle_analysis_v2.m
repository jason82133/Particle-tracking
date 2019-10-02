%
% NAME:
%               Single_particle_analysis_v2
%
% PURPOSE:
%               Calculate area, intensity, and length of the particles within defined ROIs.
%               
%               Require Particle_boundaries_tracking.m and Length_calculation.m
%               
% 
%               Written by Dr Aleks Ponjavic and Jason C Sang, University of Cambridge, 
%               2015-2016
%
%               Last updated on 2018/05/05
%       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Area, Intensity, Length, SBR] = Single_particle_analysis_v2(parameter2, images, imDataOri, mag, results, boundListN, roiListN, iListN, sbrListN)
%clear a I l A Area Intensity Length SBR track_total trackNo;
clearvars -except parameter2 images imDataOri mag results boundListN roiListN iListN sbrListN

analysis_all=parameter2.analysis_all;
savedata=parameter2.savedata;
datadir=parameter2.datadir;
pixel_size=parameter2.pixel_size;

%% Parameter setting
%analysis_all = 1;                % Set 1 for all-particle tracking; 0 for single tracking

%trackNo = 6;                   % Key in a number of particle. No need to change if tracking all particles

%savedata = 1;                    % Set 1 for saving analysis results in the same directory; 0 for not saving

%datadir = 'C:\Pile Higher and Deeper (PhD)\Area for Data\Analysis 20170720_Linda p53\';                                 % Destination for figure saving

%pixel_size = 235;                    % Diffraction limit for length calculation (nm)

%%

% Check existence of folder
check_folder=exist(datadir, 'dir');
if check_folder==0
    display('Folder not existed! Script terminated. Have you saved the figures?');
    pause
end


% For all-particle tracking

if analysis_all==1
    
    track_total=max(results(:,4));
    
    % Define format of final results
    Area=zeros(images+1, track_total+1);
    Intensity=zeros(images+1, track_total+1);
    Length=zeros(images+1, track_total+1);
    SBR=zeros(images+1, track_total+1);

    Area(1,2:end)=1:track_total;
    Intensity(1,2:end)=1:track_total;
    Length(1,2:end)=1:track_total;
    SBR(1,2:end)=1:track_total;
    
    %Area(2:end,1)=1:images;
    %Intensity(2:end,1)=1:images;
    %Length(2:end,1)=1:images;
    %SBR(2:end,1)=1:images;
    
    %% Get Length
    NList = [];
    lengthList = [];
    lengthListN = [];

    for i = 1:max(results(:,3)) % looking for frames
        display(['Analysing frame ', num2str(i), ' of ', num2str(max(results(:,3)))]);
        [i n] = size(imDataOri(:,:,i));
        R = [];
        N = [];
    
        posF = find(results(:,3)==i); % positions of the given frame in the results

        for j = 1:length(posF)
            R = [R; roiListN{posF(j)}]; % particle ROIs in the results in the ascending order of numbering
            N = [N; ones(size(roiListN{posF(j)},1),1).*results(posF(j),4)]; % corresponding particle numbers in the results
        end
        clear j posF
        
        if isempty(R)==0
            [ll, skel] = Length_calculation_v2(R, i, n, mag);
        
            [label, number] = bwlabel(skel); % thinned particles has their own numbering; look for original numbers in the following loop
        
            for num = 1:number  % find the numbering of thinning data
               ind = find(label==num);
               Rind = sub2ind(size(label), R(:,1), R(:,2));
               pos = find(ismember(Rind, ind));
               Nnum = N(pos);
               findN(num,1) = mean(Nnum);
               clear pos ind Rind Nnum
            end
            clear num skel label number R N m n
        

            lengthList = [lengthList; ll];
            NList = [NList; findN];
            %frameList = [frameList; ones(size(ll,1),1)*i];
            clear ll findN
        end
    end
    clear i

    for k = 1:max(NList)
        posN = find(NList==k);
        
        for j = 1:length(posN)
            lengthListN=[lengthListN; lengthList(posN(j))];
        end
        clear posN j
    end
    clear k
    

    for trackNo = 1:track_total
        frames = find(results(:,4)==trackNo); % frames of the given particle existing
        
        for i = 1:length(frames)
            time = results(frames(i),3);
            [m n] = size(imDataOri(:,:,time));

            %R = roiListN{frames(i)};
            %I(i) = iListN(frames(i));
            
        %% Get area
            b = boundListN{frames(i)}/mag+0.5;
            
            xb=b(:,2);
            yb=b(:,1);
            a= polyarea(xb,yb);
            A(i)=a;
            
        %% Get length
            %[ll skel] = Length_calculation_v2(R, m, n, mag);
            %l(i) = ll;

        %% Final results reformation
            Area(time+1, trackNo+1) = A(i); % area in pixel squared
            %Intensity(time+1, trackNo+1) = iListN(frames(i))/100; % real int
            %Length(time+1, trackNo+1) = pixel_size*(lengthListN(frames(i))/10+0.5); % length in nanometre
            SBR(time+1, trackNo+1) = sbrListN(frames(i))/100; % signal-to-background ratio
        
        end
    end
    
    
% for single particle tracking
elseif analysis_all==0
    
    Area=zeros(images+1, 2);
    Intensity=zeros(images+1,2);
    Length=zeros(images+1, 2);
    SBR=zeros(images+1, 2);

    Area(1,2)=trackNo;
    Intensity(1,2)=trackNo;
    Length(1,2)=trackNo;
    SBR(1,2)=trackNo;
    Area(2:end,1)=1:images;
    Intensity(2:end,1)=1:images;
    Length(2:end,1)=1:images;
    SBR(2:end,1)=1:images;

        %% Get Length
    NList = [];
    lengthList = [];
    lengthListN = [];

    for i = 1: max(results(:,3)) % looking for frame
        display(['Analysing frame ', num2str(i), ' of ', num2str(max(results(:,3)))]);
        [i n] = size(imDataOri(:,:,i));
        R = [];
        N = [];
    
        posF = find(results(:,3)==i); % positions of the given frame in the results

        for j = 1:length(posF)
            R = [R; roiListN{posF(j)}]; % particle ROIs in the results in the descending order of numbering
            N = [N; ones(size(roiListN{posF(j)},1),1).*results(posF(j),4)]; % corresponding particle numbers in the results
        end
        clear j posF
       
        [ll, skel] = Length_calculation_v2(R, i, n, mag);
        
        [label, number] = bwlabel(skel); % thinned particles has their own numbering; look for original numbers in the following loop
        
        for num = 1:number  % find the numbering of thinning data
               ind = find(label==num);
               Rind = sub2ind(size(label), R(:,1), R(:,2));
               pos = find(ismember(Rind, ind));
               Nnum = N(pos);
               findN(num,1) = mean(Nnum);
               clear pos ind Rind Nnum
        end
           clear num skel label number R N m n
        

        lengthList = [lengthList; ll];
        NList = [NList; findN];
        %frameList = [frameList; ones(size(ll,1),1)*i];
        clear ll findN
    end
    clear i

    for k = 1:max(NList)
        posN = find(NList==k);
        
        for j = 1:length(posN)
            lengthListN=[lengthListN; lengthList(posN(j))];
        end
        clear posN j
    end
    clear k
    
    
    frames = find(results(:,4)==trackNo);

    for i = 1:length(frames)
        time = results(frames(i),3);
        [m n] = size(imDataOri(:,:,time));
        %R = roiListN{frames(i)};
        
        %I(i) = iListN(frames(i));
        
        xp =R(:,2)/mag+0.5;
        yp = R(:,1)/mag+0.5;
        
        imagesc(imDataOri(:,:,time)) %Recover this part if requiring plots
        hold on %Recover this part if requiring plots
        b = boundListN{frames(i)}/mag+0.5;
        plot(b(:,2),b(:,1),'r-') %Recover this part if requiring plots
    
    %% Get area
        xb=b(:,2);
        yb=b(:,1);
        a= polyarea(xb,yb);
        A(i)=a*(pixel_size/1000)^2;
        
    %% Get length
        %[ll skel] = Length_calculation_v2(R, m, n, mag);
        %l(i) = ll;
        
        %% Plot skeleton
        %[row col] = find(skel==1);
        %hold on
        %plot(col/10+0.5,row/10+0.5,'r.')
        %display(['Frame ', num2str(i),' of particle No ', num2str(trackNo)]);
        
        %pause(0.15)
    
        %% Final results reformation
        Area(time+1, 2) = A(i); % area in pixels squared
        %Intensity(time+1, 2) = iListN(frames(i))/100; % real int
        %Length(time+1, 2) = pixel_size*(lengthListN(frames(i))/10+0.5); % in nanometre
        SBR(time+1, 2) = sbrListN(frames(i))/100; % ratio
        
    end

    figure

    subplot(1,3,1)
    plot(A)
    subplot(1,3,2)
    plot(I)
    subplot(1,3,3)
    plot(l)
end

display('All particles analysed !');


%% Export to subfolder of where the figures saved
if savedata==1
    save ([datadir 'Area.txt'], 'Area', '-ascii');
    %save ([datadir 'Intensity.txt'], 'Intensity', '-ascii');
    save ([datadir 'Length.txt'], 'Length', '-ascii');
    save ([datadir 'SBR.txt'], 'SBR', '-ascii');
    
    display(['Data saved in: "', datadir, '".']);

elseif savedata==0
    display('No data saved.');

end
close all
end