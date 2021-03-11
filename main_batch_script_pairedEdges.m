% Always to check appropriate parameters before running this script!
%
% This script is very specific for dot assay time lapse images where two
% positions for each cell sheet were imaged.  Other data should be 
% processed using the code subfolders folders.
%
% Note: this script was optimized for MATLAB 2019a. Some of the features,
% especially the plots, may be incompatible with other MATLAB versions.
%
% 2015/03/17 Created by RML
% 2018/10/20 Adapted for data from spinning disk microscope
% 2019/08/19 RML: naming convention assumes two locations per well next to
%                 each other in the file names
% 2021/03/11 RML: updated commenting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%% USER PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Folders of .tif images to be analyzed
% Can have multiple using directories = {'dir1\','dir2\','dir3\','etc'};
directories = {'E:\2021-03_multiScalePIVanalysis\ExampleImages\'};

% Repeat already completed analysis?
over_write = 0; %When over_write = 1, the edge ramen image and PIV analytics
%(basic and time) will be run again and will REPLACE existing files
overwriteEdge = 0; % 1 = update edge detection/REPLACE exisiting files

% These two frames will be used in edge_ramen plots and PIV_analytics
% (Raw PIV and edge detection will run on all available frames)
user_firstframe = 1;
user_lastframe = 240; % Can be set to inf to analyze all frames

% Imaging resolutions
r_scale = 0.582; % microns/pixel
t_scale = 3; % minutes/frame

% Filtering parameters <--robust and should not be changed within a project
size_exclude = 5000;
size_erode = 2;
size_gap = 2*size_exclude;

%For FTLE analysis:
T=40; % deformation time
trange=1:200; % frames of the movie to consider, must be <= (user_lastframe - T)

% Relative image positions
% Can be extracted from microscope stage information
% See for example convertMVD2.m for extracting PerkinElmer information
pos_type = 1; % 0 = get position from zvi files, 1 = get postion from excel sheet, 2 = analytics previously run, use those positions

% Parallel processing
parpoolName = 'local'; % Set the parallel profile to use.  
% Set to 'local' if you have not created custom profiles in MATLAB.

% Plotting parameter            
frameskip = 25; %Every frameskip-th frame will be ploted during edge detection and filter check

                 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MATLAB PREP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section makes sure the path contains the necessary subfunctions.
p = path; %Save the original path to reset later
p2 = genpath('PIV');
p3 = genpath('EdgeDetect');
p4 = genpath('Analytics');
p5 = genpath('FTLE');
p6 = genpath('CoarseGrain');
p7 = genpath('Correlations');
path(p,[p2,p3,p4,p5,p6 p7]);
clear p2 p3 p4 p5 p6

% This sets up the ability to run parfor loops.  Change to appropriate
% parpool profile as necessary.
if isempty(gcp('nocreate'))
    parpool(parpoolName);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RUN EACH STEP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The steps are run in 3 separate loops over all files:
% (1) Raw PIV, edge detection, and FTLE calculations.  These are some of
% the most time consuming steps.
% (2) A user sanity check step where edge detection is plotted.
% (3) Downstream calculations that combine images from the same well,
% calculate trends over time and space, etc.

clc, close all
disp(datestr(now)) % Display what time the batch began

% Loop #1: Raw PIV, edge detection, FTLE calculations
for kk = 1:numel(directories)  % Run each directory, then each analysis

    % Format the directory string correctly
    if ~strcmp(directories{kk}(end),'\')
        directory = [directories{kk} '\'];
    else
        directory = directories{kk};
    end
    disp(['Working on the folder: ' directory ' (' num2str(kk) ' of ' num2str(numel(directories)) ')'])

    % Find relevant files
    list = dir([directory '*.tif']);

    %%%%% PIV
    % Based on MatPIV 1.6.1
    % Sveen, Johan Kristian. “An Introduction to MatPIV v. 1.6.1,” 2004. https://www.duo.uio.no/handle/10852/10196.
    for ii = 1:numel(list) % For each tif file
        imname = list(ii).name(1:(end-4));
        if ~exist([directory imname 'dot_piv.mat'],'file') % Don't run analysis if output already exists
            disp(['     ' datestr(now) ' - Running PIV on file ' imname ' (' num2str(ii) ' of ' num2str(numel(list)) ')'])
            dot_matpiv5(directory, imname,'multiff') % Syntax: dot_matpiv5(directory, imname,image_type,firstframe,lastframe,...)
        else
            disp(['     PIV already exists for ' imname ' (' num2str(ii) ' of ' num2str(numel(list)) ')'])
        end
    end
    disp(' ')    

    %%%%% FILTERING
    % Requires MATLAB's image processing toolbox
    list = dir([directory '*dot_piv.mat']); % Find piv files in the folder
    saveimpose = [directory 'quiverSpeed\'];
    for ii = 1:numel(list) % For each piv file
        imname = list(ii).name(1:(end-4));
        if ~exist([directory imname '_filtered.mat'],'file')
            disp(['     ' datestr(now) ' - Filtering PIV on ' imname ' (' num2str(ii) ' of ' num2str(numel(list)) ')'])
            load([directory imname '.mat'])
            % The next 2 lines make sure the data will save in a nice format even if the files have been moved from their original location
            dot_piv.imname = imname(1:end-7);
            dot_piv.directory = directory;
            %  mask_and_filter removed non-cell areas from the PIV data
            dot_piv = mask_and_filter(dot_piv,size_exclude,size_erode,directory,size_gap);
            %Plot the filtered PIV
            savedir = [saveimpose 'quiverSpeed_' dot_piv.imname filesep];
            plot_piv_quiverSpeed_v2(dot_piv, 0, 1.5, r_scale, t_scale,[],[char(181) 'm'],'min',savedir)
            close all

        else
            disp(['     PIV filtering already complete on ' imname ' (' num2str(ii) ' of ' num2str(numel(list)) ')'])
        end       
    end
    disp(' ')

    %%%%% EDGE DETECTION
    % Before running the edge detection for the first time on a system,
    % run mexme_dijkstra.m to compile the source code.
    saveimpose = [directory 'EdgeImposed\'];
    if ~exist(saveimpose,'file')
        mkdir(saveimpose)
    end
    list = dir([directory '*.tif']);
    NedgeFrames = NaN*ones(length(list),1);
    for ii = 1:numel(list)
        imname = list(ii).name(1:(end-4));
        if ~exist([directory imname '_edgedat.mat'],'file') || overwriteEdge > 0
            disp(['     ' datestr(now) ' - Running edge dectection on ' imname ' (' num2str(ii) ' of ' num2str(numel(list)) ')'])
                edgedat = sheet_edger_combo_v2(imname,directory,'multiff',[],[],[],[],[],[],[],frameskip,saveimpose);
                NedgeFrames(ii) = size(edgedat.points,3); % This checks that the edge detection does not end prematurely               
        else
            disp(['     Edge detection already exists for ' imname ' (' num2str(ii) ' of ' num2str(numel(list)) ')'])
        end
    end
    save([saveimpose 'NedgeFramesCheck.mat'],'NedgeFrames')
    disp(' ')

    %%%%% EDGE 'RAMEN NOODLE' PLOTS
    % Fun plots of edge dynamics
    list = dir([directory '*edgedat.mat']);
    % Make a folder to save images
    if ~exist([directory 'LengthStacks' num2str(user_lastframe) '\'],'file')
        mkdir(directory,['LengthStacks' num2str(user_lastframe)])
    end
    saveimpose = [directory 'LengthStacks' num2str(user_lastframe) '\'];
    % Loop through files and make images
    for ii = 1:numel(list)
        imname = list(ii).name(1:(end-4));
        if ~exist([saveimpose imname 'LengthColorStack.bmp'],'file') || over_write
            load([directory list(ii).name])
            disp(['    Creating edge image for ' imname ' (' num2str(ii) ' of ' num2str(numel(list)) ')'])
            edge_ramen_length_v2(directory,imname,edgedat,[], [], saveimpose,user_firstframe,user_lastframe,uint8(parula(2^8)*((2^8)-1)),2)
        end
    end
    disp(' ')

    %%%% FTLE
    % FTLE code was developed by Doug Kelley (University of Rochester)
    % See also: Kelley, Douglas H., and Nicholas T. Ouellette. “Separating Stretching from Folding in Fluid Mixing.” Nature Physics 7, no. 6 (March 6, 2011): 477–80. https://doi.org/10.1038/nphys1941.
    savename= [directory 'FTLEall.mat'];
    if ~exist(savename,'file')
        % Find the relevant piv files
        list = [dir([directory '*piv_filtered.mat']); dir([directory '*-C*dot_piv.mat'])];
        % Format dir output more conveniently
        pivfiles = cell(numel(list),1);
        for i = 1:numel(list), pivfiles{i} = list(i).name; end
        % Useful sizes
        Npiv=numel(pivfiles);
        Nt=numel(trange);
        % Preallocate
        ftle_all=cell(Npiv,1);
        % Loop through each PIV file
        for ii=1:Npiv
            disp(['    Calculating FTLE values for ' pivfiles{ii} ' (' num2str(ii) ' of ' num2str(Npiv) ').'])
            d=load([directory pivfiles{ii}]);
            fn=fieldnames(d);
            eval(['ftle_piv=d.' fn{1} ';']);
            test_length = size(d.dot_piv.u_fi,3);
            if test_length < Nt
                Nt_now = test_length - T;
                disp([num2str(Nt) ' frames not available, using ' num2str(Nt_now)])
            else
                Nt_now = Nt;
            end
            ftle_all{ii}=NaN([size(ftle_piv.x) Nt_now]);
            for jj=1:Nt_now
                % The main calculations are done by stretchg, which is
                % a function provided by Douglas Kelly
                ftle_all{ii}(:,:,jj)=log(stretchg( ...
                    ftle_piv.u_fi(:,:,trange(jj):trange(jj)+T), ...
                    ftle_piv.v_fi(:,:,trange(jj):trange(jj)+T), ...
                    ftle_piv.x,ftle_piv.y,0,0))/T;
            end
            % Remove edge effects
            ftle_all{ii}(1:2,:,:) = NaN;
            ftle_all{ii}(end-1:end,:,:) = NaN;
            ftle_all{ii}(:,1:2,:) = NaN;
            ftle_all{ii}(:,end-1:end,:) = NaN;
        end
        save(savename,'T','pivfiles','ftle_all');
    else
        disp(['    FTLE values already exist for ' directory])
    end
    disp(' ')

end

% Loop #2: Sanity check on edge detection and filtering
%%%% FILTER CHECK
% This creates a series of images that can be used to check whether the
% edge detection and PIV filtering agree -- this provides a sanity
% check on the previous image analysis and quality of the data
need_run = 0; %  First check if filter check has already run
for kk = 1:numel(directories)
    % Format the directory correctly
    if ~strcmp(directories{kk}(end),'\')
        directory = [directories{kk} '\'];
    else
        directory = directories{kk};
    end

    % If the folder doesn't exist, it needs to be run
    if ~exist([directory 'FilterCheck\'],'file')
        need_run = 1;
        break
    end

    % Find relevant files
    list = dir([directory '*edgedat.mat']);

    for ii = 1:numel(list)
        imname = list(ii).name(1:end-12);
        if ~exist([directory 'FilterCheck\' imname '\'],'file')
            need_run = 1;
            break
        end
    end
end
%If it needs to run, wait for user to watch the filter check images
if need_run
    wait_signal = input(['**************' 10 ...
        'When you are ready to watch the filter check images, hit enter to continue' 10 ...
        '**************']);
    disp(' ')
    for kk = 1:numel(directories)
        % Format the directory correctly
        if ~strcmp(directories{kk}(end),'\')
            directory = [directories{kk} '\'];
        else
            directory = directories{kk};
        end
        disp(['Filter check on the folder: ' directory ' (' num2str(kk) ' of ' num2str(numel(directories)) ')'])

        % Create folder for saving if necessary
        if ~exist([directory 'FilterCheck\'],'file')
            mkdir(directory,'FilterCheck')
        end

        % Find relevant files
        list = dir([directory '*edgedat.mat']);

        for ii = 1:numel(list)
            imname = list(ii).name(1:end-12);
            % Only run filter check if files don't already exist
            if ~exist([directory 'FilterCheck\' imname '\'],'file')
                filter_check(directory, imname, frameskip)
            end
        end
    end
    close all
end

% Loop #3: Downstream analysis
%%%% BEGIN ANALYTICS
for kk = 1:numel(directories)
    directory = directories{kk};

    if ~strcmp(directory(end),'\')
        directory = [directories{kk} '\'];
    end
    disp(['Analyzing PIV in the folder: ' directory ' (' num2str(kk) ' of ' num2str(numel(directories)) ')'])

    %%%First run basic analytics
    % This outputs basic quantities such as speed and also creates
    % figures and a summary excel file
    basic_analyticsPairedtif(directory,r_scale, t_scale, user_firstframe, user_lastframe,pos_type,over_write)
    close all
    disp(' ')

    %%%Then run time analytics - this reports trends over time
    % Look for possible files to analyze
    list = dir([directory '*tif']);
    names_list = cell(2,numel(list)/2);
    for i = 1:numel(list)
        names_list{i} = list(i).name(1:end-4);
    end
    names_list = names_list'; % Each row is a well now
    clear list
    % Loop through the files in this directory
    for jj = 1:size(names_list,1)
        imnames = names_list(jj,:);
        imname = [imnames{1} '_' imnames{2}];
        disp(['    Time Analysis of ' imname ' (' num2str(jj) ' of ' num2str(size(names_list,1)) ')'])
        if ~exist([directory 'Analytics Data\' imname '_time_analytics.mat'],'file') || over_write
            piv_analytics_timePairedtif(directory, imnames, r_scale, t_scale, user_firstframe,user_lastframe)
        end
    end
    disp(' ')
    close all

    %%%% Now make histograms over time
    for jj = 1:size(names_list,1)
        imnames = names_list(jj,:);
        imname = [imnames{1} '_' imnames{2}];
        disp(['    Histogram Analysis of ' imname ' (' num2str(jj) ' of ' num2str(size(names_list,1)) ')'])
        if ~exist([directory 'PIV Hist\' imname '_PIV_edge.mat'],'file') || over_write
            plot_piv_time_histPairedtif(directory, imnames, r_scale, t_scale, user_firstframe,user_lastframe)
            close all
        end
    end       
    disp(' ')

    %%%% Now make spatial plots over time
    for jj = 1:size(names_list,1)
        imnames = names_list(jj,:);
        imname = [imnames{1} '_' imnames{2}];
        disp(['    r/R Analysis of ' imname ' (' num2str(jj) ' of ' num2str(size(names_list,1)) ')'])
        if ~exist([directory 'PIV Time and R\' imname 'Radial Velocity.tif'],'file') || over_write
            plot_piv_time_rPairedTif(directory, imnames, r_scale, t_scale, user_firstframe,user_lastframe)
            close all
        end
    end       
    disp(' ')        

    %%%% FTLE Analytics
    % This provides basic analysis of the FTLE data
    disp(['Analyzing FTLEs in the folder: ' directory ' (' num2str(kk) ' of ' num2str(numel(directories)) ')'])
    run_ftle_analyticsPairedTif(directory,r_scale, t_scale, user_firstframe, user_lastframe,over_write)
    close all
    disp(' ')

    %%%% Now make FTLE histograms over time
    run_ftle_time_histPairedTif(directory,r_scale, t_scale, user_firstframe, user_lastframe,over_write)
    disp(' ')

    %%%% Now plot FTLEs spatially over time
    run_ftle_time_rPairedTif(directory,r_scale, t_scale, user_firstframe, user_lastframe,over_write)
    disp(' ')

    %%%% Coarse Graining
    % See Marel, et al., Flow and Diffusion in Channel-Guided Cell
    % Migration. Biophysical Journal, 107(5):1054-1064, 2014. for a
    % discussion of this analysis
    disp(['Running coarse graining on the folder: ' directory ' (' num2str(kk) ' of ' num2str(numel(directories)) ')'])
    run_coarse_grain_pairedTif(directory,user_lastframe,r_scale,t_scale,over_write)

    %%%% Spatial Auto-Correlations
    % This is a time consuming step. Multiple types of correlations (speed,
    % radial velocity) are calculated. piv_corr_r_NaNnorm_pairedTif.m can
    % be modified to change the correlations being performed.
    % Note, when the mean is subtracted, the correlations being calculated
    % are Pearson's correlations. Modified correlations with the mean
    % included, used by the Losert lab, are also calculated.
    for jj = 1:size(names_list,1)
        imnames = names_list(jj,:);
        imname = [imnames{1} '_' imnames{2}];
        disp(['    Correlation Analysis of ' imname ' (' num2str(jj) ' of ' num2str(size(names_list,1)) ')'])
        if ~exist([directory 'Correlation Data\' imname '_vth_radial_corrNaNnorm.mat'],'file') || over_write
            piv_corr_r_NaNnorm_pairedTif(directory, imnames, r_scale, t_scale, user_firstframe, user_lastframe)
            close all
        end
    end       
    disp(' ')            

end    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CLEANUP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reset the path
path(p)
% Clean up workspace
clear kk jj ii p directory imname myaddress mypassword need_run props saveimpose savename trange
close all
% If any parallel pools are still open, close them
if ~isempty(gcp('nocreate'))
    delete(gcp)
end
% Display and e-mail completion of batch
disp(' ')
disp('Batch Complete')
disp(' ')
disp(datestr(now)) % Display what time the batch ended

