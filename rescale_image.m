function new_image = rescale_image(im,type)

if nargin == 1
    type = 8;
end


if type == 8
    
    image = double(im);
    image = image - min(image(:));
    image = image./max(image(:));
    new_image = uint8(255.*image);
    
elseif type == 16
    
    image = double(im);
    image = image - min(image(:));
    image = image./max(image(:));
    new_image = uint16(65535.*image);
    
end