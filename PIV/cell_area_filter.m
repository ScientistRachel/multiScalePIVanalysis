% 'cell_area_filter' thresholds an image of cells to distinguish between
% areas covered by cells and empty space.  This information can be use to 
% filter the output of dot_matpiv.  Thresholding is based on 'open_edger'
% by Peter Kordell.

function im_out = cell_area_filter(im,size_exclude,erode_size,gap_size)

% INPUT VARIABLES
% ~ 'im' is the matrix form of the image to be procesed
% ~ (Optional) 'size_exclude' is the smallest sized object to keep, in 
%               pixels squared.  Defaults to 1200.
% ~ (Optional) 'erode_size' defaults to zero.  If 1, will erode the image
%               to get rid of objects such as out of focus dirt.
%
% OUTPUT VARIABLES
% ~ 'im_out' is a logical image with 1's at the location of cells and 0's
%   everywhere else.
%
% 2013/04/01 Rachel Lee, file created
% 2013/04/17 RML added size_exclude and erode toggle as an input
% 2019/08/28 RML added gap_size to get rid of small holes in monolayers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% PARAMETERS %%%%%
size_exclude_def = 1200;
erode_size_def = 4;
gap_size_def = 0;

if ~exist('size_exclude','var') || isempty(size_exclude)
    size_exclude = size_exclude_def;
end
if ~exist('erode_size','var') || isempty(erode_size)
    erode_size = erode_size_def;
end
if ~exist('gap_size','var') || isempty(gap_size)
    disp('using def gap size')
    gap_size = gap_size_def;
end


im = rescale_image(im,8);
% Define the sobel matrix and perform a two dimensional convolution with the
% input image. The results of the concolutions are added in quatriture to 
% find the magnitude of the first order gradient of the image
sob = [-1 -2 -1; 0 0 0; 1 2 1];
ed1 = conv2(double(im),double(sob),'same');  %y direction gradient
ed2 = conv2(double(im),double(sob'),'same');  %x direction gradient
im_out = sqrt((ed1.^2)+(ed2.^2));

%figure;imshow(im_out), title('Sobel')

% Take the reciprocal of the gradient image and eliminate infinities
im_out = 1./im_out;
im_out(im_out==Inf) = 0;
im_out = im_out/max(im_out(:));

%figure;imshow(im_out), title('Inverse')

% Pass the reciprocated image through a median filter (this reduces the
% "salt and pepper" noise in the image).  [3 3] sets the size of the
% neighborhood over which pixels are averaged.
im_out = medfilt2(im_out, [3 3]);

% Identify the threshold to apply to the image as either a user input or
% or the result of the Otsu Method (chooses the threshold to minimize the
% intraclass variance of the black and white pixels).
thresh_out = graythresh(im_out);   %This is the function for the Otsu Method
% Apply the threshold (convert the grayscale image to a binary image)
im_out = ~im2bw(im_out,thresh_out);  %'~' is used to invert black and white

% figure;imshow(im_out)
% title('Black and White')

% Eliminate border edges by replacing with immediate neighbors 
[a, b] = size(im_out);
im_out(1,1:b) = im_out(2,1:b);
im_out(1:a,b) = im_out(1:a,(b-1));
im_out(a,1:b) = im_out((a-1),1:b);
im_out(1:a,1) = im_out(1:a,2);

if erode_size
    se = strel('disk',erode_size);
    im_out = imerode(im_out,se);
    im_out = bwareaopen(im_out,size_exclude);
    im_out = imdilate(im_out,se);
    
    im_out = imdilate(im_out,se);
    im_out = imfill(im_out,'holes');
%     im_out = imfill(im_out,8,'holes'); %Switched to 8 connectivity on 2019/08/28
    im_out = imerode(im_out,se);
else
    im_out = imfill(im_out,'holes');
end

%2019/08/28 Add requirement that monolayer not have small gaps
if gap_size > 0
    
    im_out = ~im_out;
    im_out = bwareaopen(im_out,gap_size);
    im_out = ~im_out;
    
end

% figure; imshow(im_out)
% title('Final')
