% Function to collect PIV data on phase contrast images of cells.  See file
% for detailed comments on usage, inputs, and outputs.
%
% Basic Example of Usage:
% dot_matpiv5(directory, imname, [type], [firstframe], [lastframe], [channel],[num_channel],[imname2],[savedir])
%
% WARNING: This version of dot_matpiv automatically crops the top 50 pixels
% from all images (to account for Colibri errors).  See parameters to
% change this.
%
% Version 5 - adds new file formats and updates bfopen

function dot_matpiv5(directory, imname, type, firstframe, lastframe, channel,num_channel,imname2,savedir)

% INPUTS
%   directory   - The folder containing the images to be analyzed
%   imname      - The name of the images to be analyzed
%   type        - (Optional) Defaults to 'zvi' but will also accept 'czi',
%               multipage tif files*, and image sequences in MATLAB supported
%               formats (such as tif and jpg).
%               *Specify multipage tif files as 'multiff'
%   firstframe  - (Optional) The first frame to analyze. Defaults to 1.
%   lastframe   - (Optional) The last frame to analyze. Defaults to number
%                       of images found in .zvi or in image sequence
%   channel     - (Optional) Which set of images to analyze. Defaults to 1.
%   num_channel - (Optional) Number of channels total. Defaults to 1.
%   imname2     - (Optional) If using a non-zvi image sequence, any text
%                       after the numbers in the sequence (e.g., 'c3').
%                       Defaults to [].
%   savedir     - (Optional) Save to a different directory if desired.
%                       Defaults to the image directory.
%
% NOTES: - Check parameters (see below) before using dot_matpiv on new data
%        - Can use channel and num_channel to adjust delta t if desired
%        - When using zvi files, a folder called "temp_images" is created
%          and then DELETED.  Do not store files you want to keep in a
%          folder with this name.
%
% EXAMPLE USAGE: dot_matpiv5('C:\PIV_2013-04_RLee\','PIVexample','tif')
%
% OUTPUTS
%   No outputs are sent to the matlab workspace, however the following are
%   saved to the folder specified by directory:
%   - dot_piv.mat - This file is saved as imnamedot_piv.mat and contains
%   the following structure:
%       dot_piv.x = The x locations of the PIV velocity vectors (in pixels)
%       dot_piv.y = The y locations of the PIV velocity vectors (in pixels)
%       dot_piv.u = Unfiltered velocity in the x direction
%       dot_piv.v = Unfiltered velocity in the y direction
%       dot_piv.u_f = Filtered velocity in the x direction
%       dot_piv.v_f = Filtered velocity in the y direction
%       dot_piv.u_fi = Filtered and interpolated velocity in the x direction
%       dot_piv.v_fi = Filtered and interpolated velocity in the y direction
%       dot_piv.snr = The signal to noise ratio in unfiltered velocity (used for validation)
%       dot_piv.pkh = The correlation peak height (used for validation)
%       dot_piv.date = The current date at time of saving
%       dot_piv.imname = The name entered as an input
%       dot_piv.imname2 = Extra text in image name, if used
%       dot_piv.type = Image type
%       dot_piv.channel = Channel of zvi used
%       dot_piv.num_channel = Number of channels in zvi (or delta t)
%       dot_piv.directory = The directory added as an input
%       dot_piv.firstframe = The first frame analyzed
%       dot_piv.lastframe = The last frame analyzed
%       dot_piv.created_with = Which mfile was used to run PIV;
%
%
% 2012/06/07 Rachel Lee   Adapted for M1 dot assay migration
% 2012/08/21 Rachel Lee   Cleaned up for sharing
% 2012/08/25 Rachel Lee   Added ability to load non-zvi image sequences
%                         Added automatic calculation of storage sizes
% 2013/04/10 Rachel Lee   New version (dot_matpiv2):
%                         Changed to accept fluorescent data more
%                         effienctly. Separated plotting function.
%                         Automatic registration of number of frames.
%                         Additional saving of parameters. Added a factor
%                         of two to the waitbar to account for memory usage
% 2013/09/13 Rachel Lee   Added ability to save to different directory
%                         Changed factor of 2 in waitbar to 1.5
% 2014/01/08 Rachel Lee   New version (dot_matpiv3)
%                         Use parfor to speed up the code
%                          -ZVI files are now written to temporary tiffs
% 2015/02/04 Rachel Lee   New version (dot_matpiv4)
%                          -Now explicitly handles multi image tiff files
% 2015/03/17 Rachel Lee   New version (dot_matpiv5)
%                          -Now handles czi and has updated bfopen package
%                          -Make pixel top cropping parameter (microscope
%                          specific) an explicit variable
% 2015/03/26 Rachel Lee   Fixed error on multiple channel multi-tiffs
% 2020/01/17 Rachel Lee   Fixed error with multiple channel tiffs (only
%                          corrected regular tiffs with last edit??)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The top 50 pixels on the Colibri need to be cut off.  If that is no
% longer necessary, set top_cut = 0;
top_cut = 0;

% Parameters necessary for calculating velocities
win1 = 64; % Window size for large iteration
win2 = 32; % Window size for final iteration
overlap = 0.5; % Interogation window overlap (ranges from 0 to 1, 0.5 or 0.75 recommended)

% Default values
firstframe_def = 1;
lastframe_def = inf;
channel_def = 1;
num_channel_def = 1;
type_def = 'zvi';
savedir_def = directory;

%%%%%%%%%%%%%%%%%% Take care of inputs and defaults
if nargin<2
    error(['Usage: ' mfilename '(directory, imname, [type], [firstframe], [lastframe], [channel],[num_channel],[imname2],[savedir])'])
end

if ~exist('firstframe','var') || isempty(firstframe)
    firstframe = firstframe_def;
end
if ~exist('lastframe','var') || isempty(lastframe)
    lastframe = lastframe_def;
end
if ~exist('channel','var') || isempty(channel)
    channel = channel_def;
end
if ~exist('num_channel','var') || isempty(num_channel)
    num_channel = num_channel_def;
end
if ~exist('type','var') || isempty(type)
    type = type_def;
end
if ~exist('imname2','var') || isempty(imname2)
    imname2 = [];
end
if ~exist('savedir','var') || isempty(savedir)
    savedir = savedir_def; 
end

% Optimize memory usage by clearing unneccessary variables:
clear firstframe_def lastframe_def channel_def num_channel_def type_def savedir_def

%%%%%%%%%%%%%%%%%%%%%%% TEMPORARILY SET PATH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = path;
p2 = genpath('MatPIV161');
path(p,p2)
clear p2

%%%%%%%%%%%%%%%%%%%%%%% CALCULATE VELOCITIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
if strcmp(type,'zvi') || strcmp(type,'czi')
    
    % Load the images
    h = waitbar(0,'Opening Bio-Formats Reader...');
    images = bfopen4PIV([directory imname '.' type],h);
    images = images{1};
    delete(h); % Delete the waitbar
    
    % Sort through images to only keep desired set
    in_channel = channel:num_channel:size(images,1);
    images = images(in_channel,1);
    
    % Pick reasonable last frame
    lastframe = min([lastframe size(images,1)]);
    images = images(firstframe:lastframe);
    size3 = size(images,1) - 1;
    
    % Full zvi series uses too much memory -- write to temporary tiffs
    % Make a temporary folder -- WARNING: Any folder named "temp_images" 
    % will be deleted after running dot_matpiv3, regardless of folder contents.
    if ~exist('temp_images','file')
        mkdir('temp_images')
    end
    % Write the images
    for jj = 1:(size3+1)
        imwrite(images{jj},['temp_images\temp_' num2str(jj,'%04u') '.tif'],'tiff')
    end
    clear images
    list = dir('temp_images\*.tif');
    list2 = list(2:end); %Second list increases parfor efficiency
    
    % Preallocate storage
    u = cell(1,size3);
    v = cell(1,size3);
    u_f = cell(1,size3);
    v_f = cell(1,size3);
    u_fi = cell(1,size3);
    v_fi = cell(1,size3);
    snrs = cell(1,size3);
    pkhs = cell(1,size3);
    
    % Format the first images
    im1 = imread(['temp_images\' list(1).name],'tif');
    im1 = sobel(rescale_image(im1(top_cut+1:end,:),8));
    
    im2 = imread(['temp_images\' list(2).name],'tif');
    im2 = sobel(rescale_image(im2(top_cut+1:end,:),8));

    % Run PIV
    [xw,yw,uw,vw,snr,pkh]=matpiv(im1,im2,[win1 win1; win1 win1; win2 win2; win2 win2],1,overlap ,'multin');

    % Perform Signal to Noise Filter
    [su,sv]=snrfilt(xw,yw,uw,vw,snr,1.3);

    % Interpolate NaNs
    [fu,fv]=naninterp(su,sv,'linear');

    % Save Data
    u{1} = uw;
    v{1} = vw;
    u_f{1} = su;
    v_f{1} = sv;
    u_fi{1} = fu;
    v_fi{1} = fv;
    snrs{1} = snr;
    pkhs{1} = pkh;
    
    % Create output data structure to save
    dot_piv = struct('x',[],'y',[],'u',[],'v',[],'u_f', [],'v_f', [],'u_fi', [],'v_fi', []);
    % Save vector positions
    dot_piv.x = xw;
    dot_piv.y = yw;

    clear xw yw uw vw snr pkh su sv fu fv

    % Begin calculating velocities using matpiv
    parfor i = 2:size3

        % Load images
        im1 = imread(['temp_images\' list(i).name],'tif');
        im1 = sobel(rescale_image(im1(top_cut+1:end,:),8));
        %
        im2 = imread(['temp_images\' list2(i).name],'tif');
        im2 = sobel(rescale_image(im2(top_cut+1:end,:),8));
        
        % Run PIV
        [xw,yw,uw,vw,snr,pkh]=matpiv(im1,im2,[win1 win1; win1 win1; win2 win2; win2 win2],1,overlap ,'multin');
        
        % Perform Signal to Noise Filter
        [su,sv]=snrfilt(xw,yw,uw,vw,snr,1.3);
        
        % Interpolate NaNs
        [fu,fv]=naninterp(su,sv,'linear');
        
        % Save Data
        u{i} = uw;
        v{i} = vw;
        u_f{i} = su;
        v_f{i} = sv;
        u_fi{i} = fu;
        v_fi{i} = fv;
        snrs{i} = snr;
        pkhs{i} = pkh;
         
    end
    
    if exist('temp_images','file')
        rmdir('temp_images','s')
    end
    
elseif strcmp(type,'multiff')  %%% Multipage tif files
    
    % Find correct file names
    image_info = imfinfo([directory imname '.tif']);
       
    % Pick reasonable last frame
    lastframe = min([lastframe length(image_info)/num_channel]);  
    % Sort through images to only keep desired set
    in_channel = firstframe*channel:num_channel:lastframe*num_channel;
    size3 = numel(in_channel) - 1;
    in2 = in_channel(2:end); %makes parfor happy
    
    % Preallocate storage
    u = cell(1,size3);
    v = cell(1,size3);
    u_f = cell(1,size3);
    v_f = cell(1,size3);
    u_fi = cell(1,size3);
    v_fi = cell(1,size3);
    snrs = cell(1,size3);
    pkhs = cell(1,size3);

    % Format the first images
    im1 = imread([directory imname '.tif'],in_channel(1));
    im1 = sobel(rescale_image(im1(top_cut+1:end,:),8));
    
    im2 = imread([directory imname '.tif'],in_channel(2));
    im2 = sobel(rescale_image(im2(top_cut+1:end,:),8));

    % Run PIV
    [xw,yw,uw,vw,snr,pkh]=matpiv(im1,im2,[win1 win1; win1 win1; win2 win2; win2 win2],1,overlap ,'multin');

    % Perform Signal to Noise Filter
    [su,sv]=snrfilt(xw,yw,uw,vw,snr,1.3);

    % Interpolate NaNs
    [fu,fv]=naninterp(su,sv,'linear');

    % Save Data
    u{1} = uw;
    v{1} = vw;
    u_f{1} = su;
    v_f{1} = sv;
    u_fi{1} = fu;
    v_fi{1} = fv;
    snrs{1} = snr;
    pkhs{1} = pkh;
    
    % Create output data structure to save
    dot_piv = struct('x',[],'y',[],'u',[],'v',[],'u_f', [],'v_f', [],'u_fi', [],'v_fi', []);
    % Save vector positions
    dot_piv.x = xw;
    dot_piv.y = yw;

    clear xw yw uw vw snr pkh su sv fu fv

    % Prepare for next iteration
    clear im1 im2

    % Begin calculating velocities using matpiv
    parfor i = 2:size3
        
        % Load images
        im1 = imread([directory imname '.tif'],in_channel(i));
        im1 = sobel(rescale_image(im1(top_cut+1:end,:),8));
        
        im2 = imread([directory imname '.tif'],in2(i));
        im2 = sobel(rescale_image(im2(top_cut+1:end,:),8));

        % Run PIV
        [xw,yw,uw,vw,snr,pkh]=matpiv(im1,im2,[win1 win1; win1 win1; win2 win2; win2 win2],1,overlap ,'multin');

        % Perform Signal to Noise Filter
        [su,sv]=snrfilt(xw,yw,uw,vw,snr,1.3);

        % Interpolate NaNs
        [fu,fv]=naninterp(su,sv,'linear');

        % Save Data
        u{i} = uw;
        v{i} = vw;
        u_f{i} = su;
        v_f{i} = sv;
        u_fi{i} = fu;
        v_fi{i} = fv;
        snrs{i} = snr;
        pkhs{i} = pkh;
         
    end    
    
else  %%% Image sequences that are not in .zvi format, individual images
    
    % Find correct file names
    list = dir([directory imname '*' imname2 '.' type]);
    
    % Sort through images to only keep desired set
    in_channel = channel:num_channel:numel(list);
    list = list(in_channel);
    
    % Pick reasonable last frame
    lastframe = min([lastframe numel(list)]);
    list = list(firstframe:lastframe);
    size3 = numel(list) - 1;
    list2 = list(2:end); %Second list increases parfor efficiency
    
    % Preallocate storage
    u = cell(1,size3);
    v = cell(1,size3);
    u_f = cell(1,size3);
    v_f = cell(1,size3);
    u_fi = cell(1,size3);
    v_fi = cell(1,size3);
    snrs = cell(1,size3);
    pkhs = cell(1,size3);

    % Format the first images
    im1 = imread([directory list(1).name],type);
    im1 = sobel(rescale_image(im1(top_cut+1:end,:),8));
    
    im2 = imread([directory list(2).name],type);
    im2 = sobel(rescale_image(im2(top_cut+1:end,:),8));

    % Run PIV
    [xw,yw,uw,vw,snr,pkh]=matpiv(im1,im2,[win1 win1; win1 win1; win2 win2; win2 win2],1,overlap ,'multin');

    % Perform Signal to Noise Filter
    [su,sv]=snrfilt(xw,yw,uw,vw,snr,1.3);

    % Interpolate NaNs
    [fu,fv]=naninterp(su,sv,'linear');

    % Save Data
    u{1} = uw;
    v{1} = vw;
    u_f{1} = su;
    v_f{1} = sv;
    u_fi{1} = fu;
    v_fi{1} = fv;
    snrs{1} = snr;
    pkhs{1} = pkh;
    
    % Create output data structure to save
    dot_piv = struct('x',[],'y',[],'u',[],'v',[],'u_f', [],'v_f', [],'u_fi', [],'v_fi', []);
    % Save vector positions
    dot_piv.x = xw;
    dot_piv.y = yw;

    clear xw yw uw vw snr pkh su sv fu fv

    % Prepare for next iteration
    clear im1 im2

    % Begin calculating velocities using matpiv
    parfor i = 2:size3
        
        % Load images
        im1 = imread([directory list(i).name],type);
        im1 = sobel(rescale_image(im1(top_cut+1:end,:),8));
        
        im2 = imread([directory list2(i).name],type);
        im2 = sobel(rescale_image(im2(top_cut+1:end,:),8));

        % Run PIV
        [xw,yw,uw,vw,snr,pkh]=matpiv(im1,im2,[win1 win1; win1 win1; win2 win2; win2 win2],1,overlap ,'multin');

        % Perform Signal to Noise Filter
        [su,sv]=snrfilt(xw,yw,uw,vw,snr,1.3);

        % Interpolate NaNs
        [fu,fv]=naninterp(su,sv,'linear');

        % Save Data
        u{i} = uw;
        v{i} = vw;
        u_f{i} = su;
        v_f{i} = sv;
        u_fi{i} = fu;
        v_fi{i} = fv;
        snrs{i} = snr;
        pkhs{i} = pkh;
         
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%% Save Output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find sizes to reshape
[size1, size2] = size(u{1});

% Data
dot_piv.u = reshape(cell2mat(u),[size1 size2 size3]);
dot_piv.v = reshape(cell2mat(v),[size1 size2 size3]);
dot_piv.u_f = reshape(cell2mat(u_f),[size1 size2 size3]);
dot_piv.v_f = reshape(cell2mat(v_f),[size1 size2 size3]);
dot_piv.u_fi = reshape(cell2mat(u_fi),[size1 size2 size3]);
dot_piv.v_fi = reshape(cell2mat(v_fi),[size1 size2 size3]);
dot_piv.snr = reshape(cell2mat(snrs),[size1 size2 size3]);
dot_piv.pkh = reshape(cell2mat(pkhs),[size1 size2 size3]);
% Parameters
dot_piv.date = datestr(now);
dot_piv.imname = imname;
dot_piv.imname2 = imname2;
dot_piv.type = type;
dot_piv.channel = channel;
dot_piv.num_channel = num_channel;
dot_piv.directory = directory;
dot_piv.firstframe = firstframe;
dot_piv.lastframe = lastframe;
dot_piv.created_with = mfilename;
%dot_piv.creation_time = datevec(ending-start);

save(sprintf('%s%s%s',savedir,imname,'dot_piv.mat'),'dot_piv')

%%%%%%%%%%%%%%%%%%%%%%% RESET PATH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
path(p)