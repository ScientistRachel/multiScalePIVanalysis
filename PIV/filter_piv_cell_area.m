% Uses a mask to filter the velocity information in dot_piv.mat
%
% Usage: dot_piv_filtered = filter_piv_cell_area(dot_piv,im_mask,save_dir)

function dot_piv_filtered = filter_piv_cell_area(dot_piv,im_mask,save_dir)

%%%%% INPUTS
%  dot_piv = output from dot_matpiv.m
%       (see that function for details)
%  im_mask = binary matrix used to filter the piv data
%       See cell_area_filter for a simple threshold version
%  save_dir (optional) = directory to save filtered piv .mat file
%%%%% OUTPUTS
%  dot_piv_filtered = structure with same fields as dot_piv, but with NaN
%       in regions outside of mask
%
% Based on 'filter_piv', created by Rachel Lee, 2012/12/07
% Edited to work based on an area mask, 2013/04/01, RML
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 2
    save_dir = dot_piv.directory;
end


%Preallocate storage
mask = ones(size(dot_piv.u_fi));

max_frame = min(size(dot_piv.u_fi,3),numel(im_mask));

[X, Y] = meshgrid(1:size(im_mask{1},2),1:size(im_mask{1},1));

% For each frame find interpolate mask to PIV grid
for i = 1:max_frame
    
    IN = interp2(X,Y,double(im_mask{i}),dot_piv.x,dot_piv.y);
%     IN = double(imfill(IN>0,8,'holes')); %2019/08/28
    IN(IN == 0) = NaN;  %Make locations outside cell area not a number
    
    mask(:,:,i) = IN;
    
end

% Filter raw piv with mask
dot_piv.u = dot_piv.u.*mask;
dot_piv.v = dot_piv.v.*mask;
dot_piv.u_f = dot_piv.u_f.*mask;
dot_piv.v_f = dot_piv.v_f.*mask;
dot_piv.u_fi = dot_piv.u_fi.*mask;
dot_piv.v_fi = dot_piv.v_fi.*mask;
dot_piv.snr = dot_piv.snr.*mask;
dot_piv.pkh = dot_piv.pkh.*mask;

dot_piv.is_filtered = 1;

% Save data
dot_piv.date = datestr(now);
save(sprintf('%s%s%s',save_dir,dot_piv.imname,'dot_piv_filtered.mat'),'dot_piv')
% Output
dot_piv_filtered = dot_piv;