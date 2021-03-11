% This function, with its nested functions, finds the leading edges of an
% image sequence that depicts an undulating cell sheet.
% Usage: edgedat = sheet_edger_combo(imname,directory,im_format,intframe,finframe,diskrad,threshm,thresh_in,threshs,savedir,frameskip,saveimpose)

% This functions nested functions are:
% -- open_edger
% -- matdijkstra
%   -- dijkstra_mat
% -- mexdijkstra
%   -- Kadjacency
%   -- dijkstra
% -- cutter2

function  edgedat = sheet_edger_combo(imname,directory,im_format,intframe,finframe,diskrad,threshm,thresh_in,threshs,savedir,frameskip,saveimpose)

% INTPUT VARIABLES
% ~ 'imname' is the base name of the input image sequence. For example, a
%       sequence of enumerated images with the names 'examplename###.tif'
%       has a base name of 'examplename'
% ~ 'directory' contains the directory location of the image sequence. It
%       must end with a / or \ character depending on the storage format.
% ~ 'intframe' is the first desired frame of the enumerated image sequence
% ~ 'finframe' is the last desired frame of the enumerated image sequence
% ~ 'im_format' specifies the format of the image sequence. jpg, bmp, tif,
%       etc.
% ~ 'diskrad' determines the radius of the disk used to morphologically
%       open the image mask.
% ~ 'threshm' specifies the thresholding technique used the generate
%       a mask from the edge enhanced image. There are two acceptable inputs,
%       'auto', which specifies Otsu Method thresholding and 'manual' which
%       uses a user defined threshold.
% ~ 'thresh_in' is the constant, user defined threshold used when
%       'threshm' is set to manual.
% ~ 'threshs' specifies a constant bias in the Otsu threshold if
%       'threshm' is set to 'auto'. Negative values bias towards higher
%       degradation and positive values bias towards lower degradation
% ~ 'savedir' saves edgedat to a different directory
% ~ 'frameskip' specifies how often to plot an overlay
%       (see SuperImpose_Edge_crop)
% ~ 'saveimpose' specifies where to save the overlays

% OUTPUT VARIABLES
% 'edgedat' is a data structure containing the following fields
% ~~ '.points' is an N by 2 by C uint16 array that stores the pixel coordinates
%       of the edges of all images in the sequence. N is the length of the
%       longest edge in the sequence, and C is the number of frames analyzed
% ~~ '.diskrad' stores the input variable diskrad for later reference.
% ~~ '.threshm' stores the input variable 'threshm' for later reference
% ~~ '.threshv' only exists in the output if the user chooses a constant
%       threshold. It stores the input threshold.
% ~~ '.threshs' stores the input variable 'threshs' for later reference.
% ~~ '.side' stores a string saying wether the cells start on the right or
%       left side of the image.
% ~~ '.imname' stores the input variable 'imname' for later reference
% ~~ '.imsize' contains the m by n resolution of the source images
% ~~ '.info' contains the data structure containing all available
%       information on the source images. For example, resolution, bit
%       depth, publisher, etc.
% ~~ '.date' contains the date and time that the output data structure was
%       created.

% LOG
% File created by Peter Kordell
% Edited by Peter Kordell, 8/13/2011
%   increased the preallocated memory to edgedat.points (line 206) at the
%   time...
% Edited by Peter Kordell, 9/21/2011
%   Added support for zvi files.
%   Default methods for automatically identifying image sequences added
%   Input paramaters rearranged to allow for the method above

% Edited by Rachel Lee, May 2011
% Edited by Rachel Lee, 5/21/2012
%   Changed the way the output was saved (now goes to same directory as
%   images)
% Same as sheet_edger, but crops off faulty microscope pixels! RML 2012/08/09 
% Edited by Rachel Lee, 7/26/2013
%   Fixed frame specification for zvi -- previously kept extra zeros when
%   starting from a frame different than 1.
% Edited by Rachel Lee, 9/30/2013
%   Allowed saving in a different directory.  Changed to read an image list
%   automatically (instead of loading in images by strings with known
%   size). Changed the way default inputs were controlled (originally based
%   on if nargin == x, now based on each variable).
% Edited by Rachel Lee, 3/17/2015
%   Cleaned up formatting, commenting
%   Added additional file type (czi) and updated bfopen.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% PARAMETERS:
diskrad_def = 12;
threshm_def = 'auto';
thresh_in_def = 'NA';
threshs_def = 0;
savedir_def = directory;
top_cut = 0; %This crops off the bad pixels on the Colibri, set to zero if unnecessary

% set a limited number of default paramaters if the user does not specify
if(nargin<3)
    error('ERROR!! Enter at least three input arguments: sheet_edger_combo(imname,directory,im_format,intframe,finframe,....');
end
if ~exist('diskrad','var') || isempty(diskrad)
    diskrad = diskrad_def;
end
if ~exist('threshm','var') || isempty(threshm)
    threshm = threshm_def;
end
if ~exist('thresh_in','var') || isempty(thresh_in)
    thresh_in = thresh_in_def;
end
if ~exist('threshs','var') || isempty(threshs)
    threshs = threshs_def;
end
if ~exist('savedir','var') || isempty(savedir)
    savedir = savedir_def;
end

% Process edge detection based on image type
% First Zeiss formats:
if strcmp(im_format,'zvi') || strcmp(im_format,'czi')
    
    % if the zvi file doesn't exist, throw an error
    if(~exist([directory,imname,'.',im_format],'file'));
        error(['ERROR!!! Was not able to locate the file:',...
            directory, imname, '.', im_format]);
    end
    
    % start the wait bar
    h = waitbar(0,'Opening Bio-Formats Reader...');
    % read in the zvi file into matlab
    bfIms = bfopen4PIV([directory,imname,'.',im_format],h);
    % extract the images and the metadata from the file
    inf = bfIms{2};
    % get rid of all other non-image data
    bfIms = bfIms{1};
    
    [Height,Width] = size(bfIms{1});
    % either get the number of frames in the file or check that the user
    % number of frames is in agreement with the number of frames
    if ~exist('intframe','var') || isempty(intframe)
        intframe = 1;
    end
    if ~exist('finframe','var') || isempty(finframe)
        finframe = size(bfIms,1);
    end    


else
    
    %start the waitbar
    h = waitbar(0,'Checking Image Sequence');
    % Set useable frames   
    if ~exist('intframe','var') || isempty(intframe)
        intframe = 1;
    end
    if ~exist('finframe','var') || isempty(finframe)
        finframe = size(dir([directory imname '*.' im_format]),1);
    end   
    % Find the list of images   
    img_list = dir([directory imname '*.' im_format]);
    % acuire image info
    inf = imfinfo([directory img_list(intframe).name]);
    %inf = imfinfo([directory sprintf('%s%03u.%s',imname,intframe,im_format)]);
    Height = inf.Height;
    Width = inf.Width;
    
end

delete(h)

% declare the output data structure
edgedat = struct('points',[], 'diskrad', [], 'threshm',[], 'threshv',[],...
    'threshs', [], 'side',[],'imname',[],'imsize',[], 'info',[],'date',[]);

% store available data in edgedat
edgedat.imsize = [Height-top_cut, Width]; %Subtract top_cut because cropping later!
edgedat.imname = imname;
edgedat.date = date;
edgedat.info = inf;
edgedat.diskrad = diskrad;
if strcmp(threshm,'auto')
    edgedat.threshm = 'auto';
    edgedat.threshs = threshs;
    edgedat = rmfield(edgedat,'threshv');
else
    edgedat.threshm = threshm;
    edgedat.threshv = thresh_in;
    edgedat.threshs = threshs;
end
edgedat.firstframe = intframe;

% Preallocate Memory for output field '.points'
m = 11000; %2014/02/17 increased from 9000 %2014/05/15 increased from 10000
edgedat.points = zeros(m,2,finframe-intframe+1,'uint16'); %RML 2013/07/26

h = waitbar(0,'Begin Edge Detection');

% Implement auto thresholding edging
if strcmp(threshm,'auto');
    
    % Start up the waitbar with an initial guess of 2 seconds per fame
    % waitbar(0,h,'Processing Frame: 1');
    t_ave = 2;
    
    % Execute code only for images contained by 'intframe' and 'finframe'
    frame_count = 1;
    for frame = intframe:finframe
        
        % begin recording the average time of computation for this frame
        tic
        
        % calculate the estimated time of finishing the sequence in minues
        % and seconds, and display the estimate in the waitbar
        mins = floor(t_ave*(finframe-frame+1)/60);
        secs = round(floor((t_ave*(finframe-frame+1))) - mins*60);
        waitbar((frame-intframe)/(finframe-intframe+1),h,sprintf('%s%u%s%u%s%u',...
            'Processing frame: ',...
            frame,...
            '      Est time remaining: ',...
            mins,...
            'mins ',...
            secs,...
            'secs'));
        
        % load the image of number 'frame' into the workspace
        if strcmp(im_format,'zvi') || strcmp(im_format,'czi')
            im = bfIms{frame,1};
        else
            im = imread([directory img_list(frame).name],im_format);
        end        
        im = im(top_cut+1:end,:); %crops top pixels if necessary from image
        
        % if 'auto' is specified, find edge image with function 'open_edger'  using Otsu thresholding
        [im_edge, thresh_out, lr] = open_edger(im,diskrad,'auto','NA',threshs);
        
        % refuse to accept a automatic threshold of 0. use previous
        % threshold if edger assigned a threshold of less than one
        if (thresh_out) <= .001
            thresh_old = 0.001;
            [im_edge, thresh_old, lr] = open_edger(im,diskrad,'manual',thresh_old,threshs);
                  
            % if the output threshold is greater than 1, save in the variable
            % thresh_old
        else
            thresh_old = thresh_out;
        end
        
        % This piece of code added to addess the "numel(s) = 0" error
        % It used to be inside the (thresh_out) <= 0.001 loop
        %(RML 2013/07/24)
            test_thresh = find(im_edge == 1);
            if isempty(test_thresh)
                disp('Bad Thresholding!')
                [im_edge, thresh_old, lr] = open_edger(im,diskrad,'manual',0.002,threshs);
            end  
        
        % find the shortest path along the true pixels of the image with
        % function 'shortest_edge_fast' and assign sequence to 'coords'
        coords = mexdijkstra(im_edge);
        
        % if the mex function dijkstra function returns an empty set
        % (there's a bug in there somewhere...) use the old .m dijkstra function instead
        a = size(coords,1);
        if(~a)
            display(sprintf('%s%u','WARNING!: .mex dijkstra algorithm failed on frame ', frame))
            display('   Trying .m version of dijkstra algorithm...')
            clear coords
            coords = matdijkstra(im_edge);
            a = size(coords,1);
            if(a)
                display(sprintf('%s%u','   .m version succeeded in processing frame ', frame))
            elseif(~a)
                display(sprintf('%s%u','   .m version failed in processing frame ', frame))
                display(sprintf('%s%u%s','   WARNING!: edgedat.points(:,:,',frame,') is empty'))
            end
        end
        
        % add zeros to the output of of coords so that it fits in the
        % preallocated memory points_seq
        coords_resize = [coords; zeros(m-a, 2)];
        edgedat.points(:,1,frame_count) = coords_resize(:,1);  %IF subscripted assignment mismatch on this line:
        % Increase the size of m (too complicated of an edge).
        edgedat.points(:,2,frame_count) = coords_resize(:,2);
        
        % halt edging process when any part of the edge touches the left
        % or right hand side of the image. Prevents an incomplete edge
        % RML 2014/06/12 This is buggy; I'm not sure why!
        a = sum(im_edge(:,2)==1);
        b = sum(im_edge(:,Width-1)==1);
%         if frame > 200 %Use this to skip a frame or a few manually - ONLY IF NECESSARY
        if (a>0) || (b>0) 
            display(sprintf('%s%u','The cells touched the edge on image ',frame))
            break
        end
%         end
        
        if ~mod(frame-1,frameskip)
            %plot it
            points_seq_slice = cutter2([coords_resize(:,1) coords_resize(:,2)]);
            SuperImpose_Edge_combo(imname, points_seq_slice, im ,saveimpose,frame)
        end


        % record running average computation time per frame
        t = toc;
        t_ave = (t_ave+t)/2;
        
        frame_count = frame_count+1;
        
    end
    
% Implement user chosen paramater edging
elseif strcmp(threshm,'manual')
    
    h = waitbar(0,'Processing Frame: 1');
    
    t_ave = 2;
    
    % Exacute code only for images contained by 'intframe' and 'finframe'
    frame_count = 1;
    for frame = intframe:finframe
        
        tic
        
        mins = floor(t_ave*(finframe-frame+1)/60);
        secs = round(floor((t_ave*(finframe-frame+1))) - mins*60);
        waitbar((frame-intframe)/(finframe-intframe+1),h,sprintf('%s%u%s%u%s%u',...
            'Processed frame: ',...
            frame,...
            '      Est time remaining: ',...
            mins,...
            'mins ',...
            secs,...
            'secs'));
        
        % load the image of number 'frame' into the workspace
        if strcmp(im_format,'zvi') || strcmp(im_format,'czi')
            im = bfIms{frame,1};
        else
            im = imread([directory img_list(frame).name],im_format);
        end
        im = im(top_cut+1:end,:); %crops top pixels if necessary from image
        
        [im_edge, thresh_out, lr] = open_edger(im,diskrad,'manual',thresh_in,threshs);
        
        % find the shortest path along the true pixels of the image with
        % function 'shortest_edge_fast' and assign sequence to 'coords'
        coords = mexdijkstra(im_edge);
        
        % add zeros to the output of of coords so that it fits in the
        % preallocated memory points_seq
        a = size(coords,1);
        coords_resize = [coords; zeros(m-a, 2)];
        edgedat.points(:,1,frame_count) = coords_resize(:,1);
        edgedat.points(:,2,frame_count) = coords_resize(:,2);
        
        % halt edging process when any part of the edge touches the left
        % or right hand side of the image. Prevents an incomplete edge
        a = sum(im_edge(:,2)==1);
        b = sum(im_edge(:,inf.Width-1)==1);
        if (a>0) || (b>0)
            display(sprintf('%s%u','The cells touched the edge on image ',frame))
            break
        end
        
        if ~mod(frame-1,frameskip)
            %plot it
            points_seq_slice = cutter2([coords_resize(:,1) coords_resize(:,2)]);
            SuperImpose_Edge_combo(imname, points_seq_slice, im ,saveimpose,frame)
        end
        
        % record running average computation time per frame
        t = toc;
        t_ave = (t_ave+t)/2;
        
        frame_count = frame_count+1;
        
    end
    
end

waitbar(1,h,'Saving...')

% store wether or not the cells were on the right or left
edgedat.side = lr;

% eliminate empty space from preallocated variable points_seq
% make index unit8
edgedat.points = edgedat.points(:,:,1:frame_count-1);
[a, ~] = find(edgedat.points);
edgedat.points = edgedat.points(1:max(a),:,:);

% Save file
save(sprintf('%s%s%s',savedir,imname,'_edgedat.mat'),'edgedat')

delete(h);