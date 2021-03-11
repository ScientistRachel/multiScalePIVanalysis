function ims = sobel(im)

% enhance edges in the image using sobel in both directions and adding them
ha = fspecial('sobel');
hb=ha';
im1a = imfilter(im,ha,'replicate');
im1b = imfilter(im,hb,'replicate');
ims=sqrt(double (im1a).^2+double (im1b).^2);