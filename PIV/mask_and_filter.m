% Creates a mask by taking the gradient of the original images, then uses
% this to filter PIV data.
% Usage: dot_piv_filtered = mask_and_filter(dot_piv,size_exclude,erode_size)

function dot_piv_filtered = mask_and_filter(dot_piv,size_exclude,erode_size, save_dir,gap_size)

%%%%% PARAMETERS:
top_cut = 0; %This many pixels will be cut off the top.  This was added to deal with a specific microscope error and can be set to zero if necessary.

% Change Log
% 2014/03/11 RML -- added cropping of top 50 pixels to images to account
% for Colibri microscope errors.
% 2015/03/17 RML -- made cropping of top 50 pixels easier to change
%                   changed file types available and updated bfopen
%                   in addition to zvi and image sequences, can now use czi
%                       and multi page tiffs
% 2015/03/26 RML -- Fixed error on multiple channel multi-tiffs
% 2019/08/28 RML -- added gap_size as a parameter in cell_area_filter

%%%% Set optional parameters to function defaults
if ~exist('size_exclude','var') || isempty(size_exclude)
    size_exclude = [];
end
if ~exist('erode_size','var') || isempty(erode_size)
    erode_size = [];
end
if ~exist('save_dir','var') || isempty(save_dir)
    save_dir = dot_piv.directory;
end
if ~exist('gap_size','var') || isempty(gap_size)
    gap_size = [];
end

%%% Make mask based on original images
if strcmp(dot_piv.type,'zvi')  || strcmp(dot_piv.type,'czi')
    
    % Load the images
    h = waitbar(0,'Opening Bio-Formats Reader...');
    images = bfopen4PIV([dot_piv.directory,dot_piv.imname,'.zvi'],h);
    images = images{1};
    delete(h); % Delete the waitbar
    disp(' ')

    in_channel = dot_piv.channel:dot_piv.num_channel:size(images,1);
    images = images(in_channel,1);
    images = images(dot_piv.firstframe:dot_piv.lastframe);

    % matlabpool local
    im_mask = cell(1,size(images,1));

    parfor i = 1:size(images,1)

        im_mask{i} = cell_area_filter(images{i}(top_cut+1:end,:),size_exclude,erode_size,gap_size);

    end

    clear images;
    
elseif strcmp(dot_piv.type,'multiff')
    
    % Find correct file names
%     image_info = imfinfo([dot_piv.directory dot_piv.imname '.tif']);
    channel = dot_piv.channel;
    num_channel = dot_piv.num_channel;
    
    % Sort through images to only keep desired set
    in_channel = dot_piv.firstframe*channel:num_channel:dot_piv.lastframe*num_channel;
%     in_channel = channel:num_channel:length(image_info);
%     in_channel = in_channel(firstframe:dot_piv.lastframe);
    
    im_mask = cell(1,numel(in_channel));
    direc = dot_piv.directory;
    imname = dot_piv.imname;
    
    parfor i = 1:numel(in_channel)
        im = imread([direc imname '.tif'],in_channel(i));
        im = im(top_cut+1:end,:);
        im_mask{i} = cell_area_filter(im,size_exclude,erode_size,gap_size);        
    end
    
else
    
    list = dir([dot_piv.directory dot_piv.imname '*' dot_piv.imname2 '.' dot_piv.type]);
    
    in_channel = dot_piv.channel:dot_piv.num_channel:numel(list);
    in_channel = in_channel(dot_piv.firstframe:dot_piv.lastframe);
    list = list(in_channel);
    
    im_mask = cell(1,numel(in_channel));
    direc = dot_piv.directory;
    
    parfor i = 1:numel(list)
        im = imread([direc list(i).name]);
        im = im(top_cut+1:end,:);
        im_mask{i} = cell_area_filter(im,size_exclude,erode_size,gap_size);        
    end
    
end

% Filter dot_piv using masks
dot_piv_filtered = filter_piv_cell_area(dot_piv, im_mask, save_dir);
