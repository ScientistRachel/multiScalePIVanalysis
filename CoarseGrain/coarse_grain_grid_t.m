% Coarse grain a grid, over time
% Usage: grid_c = coarse_grain_grid_r(grid_values,coarse_size)
%
% INPUTS
% grid_values = grid to be coarse grained
% coarse_size = number of frames to average over
%               each frame takes into account coarse_size/2 frames before
%               and after the frame in consideration
% OUTPUS
% grid_c      = Coarse grained grid
%
% note: NaN values result in the entire grid point becoming NaN
%
% RML 2014/08/06

function grid_c = coarse_grain_grid_t(grid_values,coarse_size)

% Center the averaged frame in a bit of size coarse_size
coarse_dist = floor(coarse_size/2);

% Preallocate
[a,b,c] = size(grid_values);
grid_new = NaN*ones(a,b,length((1+coarse_dist):(c-coarse_dist)));

count = 1;
% Loop through all frames with enough frames on either side for complete
% averaging
for kk = (1+coarse_dist):(c-coarse_dist)
    
    slice = grid_values(:,:,(kk-coarse_dist:kk+coarse_dist));
    grid_new(:,:,count) = mean(slice,3);
    count = count+1;
    
end

grid_c = grid_new;