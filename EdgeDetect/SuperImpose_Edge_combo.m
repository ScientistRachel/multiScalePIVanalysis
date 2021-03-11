% this function superimposes the edge contained in the edgedat structure on
% a (usually the source sequence) image sequence.
% Only superimpose a single image which is now an input - RML, 2014/02/08

function SuperImpose_Edge_combo(imname, points_seq_slice, im_cells,savedir,frame)

% INPUT VARIABLES
% ~ 'imname' contains base name of the image sequence (see imname in
%       sheet_edge
% ~ 'directory' defines the location of the image sequece, must end in a /
%       or \ character depending on the opperating system
% ~ 'format' specifies the format of the image sequence. jpg, bmp, tif,
%       etc.
% ~ 'edgedat' the output data structure from the function sheet_edge
% ~ 'frameskip' how often to plot frames (i.e. every 10th frame).
% ~ 'firstframe' the first frame of the movie to plot (not necessarily in
%       line with edgedat.points(:,:,1))
% ~ 'savedir' allows for saving to a different directory

% OUTPUT VARIABLES
% ~ '' null. This function writes a series of image files in the current
%       MATLAB directory

% LOG
% File created by Peter Kordell
% Edited by Peter Korde, 10/13/2011
%   added support for zvi files
%   made the save directories multiplatform
% Edited by Rachel Lee, May 2011
% 2013/10/02 RML read images using dir instead of three digit file names
% 2014/02/05 RML add support for starting not on the first frame
% 2014/02/08 RML clean up skipping frames, add different save directory
% 2015/03/17 RML cleaned up comments/formatting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% PARAMETERS
% Define the color to be written on the edge (free to be changed by the
% user. Green seems to work best, so these values were not included
% in the input paramaters).
red_intensity = 0/255;
gre_intensity = 255/255;
blu_intensity = 0/255;

% see what filesystem is being used, so that directories are handled
% properly
if ispc
    dirchar = '\';
else
    dirchar = '/';
end

% define a default colormap
cmap = repmat(linspace(0,1,255)',1,3);

% Format image correctly
im_cells = rescale_image(im_cells,8);
im_cells_RGB = ind2rgb(im_cells, cmap);
        
% place a green dot in the image im_cells_RGB at the coordinates
% specified in points_seq_slice
q = size(points_seq_slice,1);
for i = 1:q
    im_cells_RGB(points_seq_slice(i,1),points_seq_slice(i,2),1) = red_intensity;
    im_cells_RGB(points_seq_slice(i,1),points_seq_slice(i,2),2) = gre_intensity;
    im_cells_RGB(points_seq_slice(i,1),points_seq_slice(i,2),3) = blu_intensity;
end
        
%create a folder for the images, if necessary, and write the images
A = exist(sprintf('%s%s%s',savedir,'EdgeImposed_',imname),'file');

if ~A
    mkdir(sprintf('%s%s%s',savedir,'EdgeImposed_',imname))
end
        
% write the image with the edge written on top
imwrite(im_cells_RGB,[savedir 'EdgeImposed_' imname dirchar 'SuperImposed_' imname num2str(frame,'%03u') '.tif'],'tif');
