% Get stage location from zvi file metadata
% Usage: [x, y] = get_zvi_location(filename)

function [x, y] = get_zvi_location(filename)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS
% filename = full path to zvi file
% OUTPUTS
% x,y      = x and y coordinates of stage location (in microns)
%
% Created by Rachel Lee, 2014/02/04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load the images
h = waitbar(0,'Opening Bio-Formats Reader...');
images = bfopen4PIV(filename,h,1,1); %Only open the first frame for speed
delete(h)

% Get metadata
metadata = images{2};
clear images

% Grab stage locations
x = str2double(metadata.get('X position for position #1'));
y = str2double(metadata.get('Y position for position #1'));