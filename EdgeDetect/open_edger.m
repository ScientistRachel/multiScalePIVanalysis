% 'open_edger' creates a logical image showing the coordinates of the
%  leading edge from an image of a cell group using a combination of the
%  sobel gradient filter, an intensity reciprocal, and either a
%  manual or automated otsu thresholding.

function [im_edge,thresh_out,lr] = open_edger(im,smooth_size,threshm,thresh_in,threshs)

% INPUT VARIABLES
% ~ 'im' is the matrix form of the image to be procesed
% ~ 'smooth_size' sets the radius of the disk used in the morphological
%    opening
% ~ 'threshm' specifies the thresholding technique used the generate
%    a mask from the edge enhanced image. There are two acceptable inputs,
%    'auto', which specifies otsu thresholding and 'manual', which uses a
%    user defined threshold
% ~ 'thresh_in' is the constant threshold used when 'threshm' is set to
%    maunal
% ~ 'threshs' specifies a constant bias in the otsu threshold if
%   'threshm' is set to 'auto'. Negative values bias towards higher
%    degradation and positive values bias towards lower degradation

% OUTPUT VARIABLES
% ~ 'im_edge' is a logical image with 1's at the location of the leading
%       edge and 0's everywhere else.
% ~ 'thresh_out' returns the threshold used in the image
% ~ 'lr' returns whether the cells were on the left or the right side of
%       the detected edge

% LOG
% File created by Peter Kordell
% Edited by Peter Kordell, 9/21/2011
%   turned off warnings that prevented division by zero
% Edited by Rachel Lee, May 2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% disable warnings that result form division by zero, and imopen
warning off

% set a limited number of default paramaters if the user does not specify
if(nargin==1)
    smooth_size = 12;
    threshm = 'auto';
    thresh_in = 'NA';
    threshs = 0;
elseif(nargin==2)
    threshm = 'auto';
    thresh_in = 'NA';
    threshs = 0;
elseif(nargin==3)
    thresh_in = 'NA';
    threshs = 0;
elseif(nargin==4)
    threshs = 0;
end


% Define the sobel matrix and perform a two dimensional convolution with the
% input image. The results of the concolutions are added in quatriture to 
% find the magnitude of the first order gradient of the image
sob = [-1 -2 -1; 0 0 0; 1 2 1];
ed1 = conv2(double(im),double(sob),'same');  %y direction gradient
ed2 = conv2(double(im),double(sob'),'same');  %x direction gradient
im_edge = sqrt((ed1.^2)+(ed2.^2));

% Take the reciprocal of the gradient image and eliminate infinities
im_edge = 1./im_edge;
im_edge(im_edge==Inf) = 0;

% Pass the reciprocated image through a median filter (this reduces the
% "salt and pepper" noise in the image).  [3 3] sets the size of the
% neighborhood over which pixels are averaged.
im_edge = medfilt2(im_edge, [3 3]);
im_edge = im_edge-min(im_edge(:));
im_edge = im_edge./max(im_edge(:));

% Identify the threshold to apply to the image as either a user input or
% or the result of the Otsu Method (chooses the threshold to minimize the
% intraclass variance of the black and white pixels).
if strcmp(threshm,'auto')

    thresh = graythresh(im_edge);   %This is the function for the Otsu Method
    thresh_out = (thresh-(threshs));

elseif strcmp(threshm,'manual')
    
    thresh_out = thresh_in - threshs;

end

% Apply the threshold (convert the grayscale image to a binary image)
im_edge = ~im2bw(im_edge,thresh_out);  %'~' is used to invert black and white

% Eliminate border edges by replacing with immediate neighbors 
[a, b] = size(im_edge);
im_edge(1,1:b) = im_edge(2,1:b);
im_edge(1:a,b) = im_edge(1:a,(b-1));
im_edge(a,1:b) = im_edge((a-1),1:b);
im_edge(1:a,1) = im_edge(1:a,2);

% Fill iosolated regions of both image objects
im_edge = bwfill(im_edge,'holes');
im_edge = ~im_edge;
im_edge = bwfill(im_edge,'holes');
im_edge = ~im_edge;

% Perform Morphologicl Opening
se = strel('disk',smooth_size);  %smooth size sets the radius of the morphological structuring disk
im_edge = imopen(im_edge,se);


% % Identify if the edge has cells on the right or left hand side
% if(im_edge(1,b))
%     lr = 'right';
% else
%     lr = 'left';
% end


% Use the first frame to decide which side the cells will be on.
sides(:,1) = sum(im_edge(:,1)); % Left
sides(:,2) = sum(im_edge(:,end)); % Right
sides(:,3) = sum(im_edge(1,:)); % Top
sides(:,4) = sum(im_edge(end,:)); % Bottom

% disp(sides)
sideEmpty = find(sides == min(sides));
if sideEmpty == 1
    lr = 'left';
elseif sideEmpty == 2
    lr = 'right';
elseif sideEmpty == 3
    lr = 'top';
elseif sideEmpty == 4
    lr = 'bottom';
else
    warning('Monolayer Side Not Detected')
    lr = 'unknown';
end

    
%find the outline of the mask
im_edge = bwperim(im_edge,8);

%deletes the borderlines of the outline
im_edge(1,:)=0;
im_edge(a,:)=0;
im_edge(:,1)=0;
im_edge(:,b)=0;

%reassigns the upper and lower borders to be their closest row nieghbors
im_edge(1,(im_edge(2,:)==1)) = 1;
im_edge(a,(im_edge(a-1,:)==1)) = 1;

%eliminates all edge bodies of fewer than 1000 connections
im_edge = bwareaopen(im_edge, size(im_edge,1));

% re enable warnings
warning on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%