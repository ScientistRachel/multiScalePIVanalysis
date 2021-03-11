% Coarse grain a grid, spatially
% Usage: grid_c = coarse_grain_grid_r(grid_values,coarse_size)
%
% INPUTS
% grid_values = grid to be coarse grained
% coarse_size = number of vectors to average over
%               coarse_size(1) averages over the rows, while coarse_size(2)
%               averages over the columns.  If a single value, the grid is
%               averaged isotropically.
% OUTPUS
% grid_c      = Coarse grained grid
%
% note: NaN values result in the entire grid point becoming NaN
%
% RML 2014/08/05
% 2014/08/07 RML -- changed loop over time to a parfor loop
% 2016/04/25 RML -- changed to keep from cutting unnecessary cols/rows off of the grid 

function grid_c = coarse_grain_grid_r(grid_values,coarse_size)

if length(coarse_size) == 1
    coarse_size(2) = coarse_size(1);
end
A = coarse_size(1);
B = coarse_size(2);

[a,b,c] = size(grid_values);

% Set up bins to average; don't keep bins too far out
% This will result in lost vectors near the far edges
row_bits = 1:A:a;
if row_bits(end)+A > a+1 %%%% 2016/04/25 a -> a + 1;
    row_bits(end) = [];
end

col_bits = 1:B:b;
if col_bits(end)+B > b+1 %%%% 2016/04/25 b -> b + 1;
    col_bits(end) = [];
end

% grid_new = NaN*ones(length(row_bits),length(col_bits),c);
grid_new = cell(c,1);

parfor kk = 1:c
    
    slice = grid_values(:,:,kk);
    col_count = 1;
    
    temp_grid = NaN*ones(length(row_bits),length(col_bits));
  
    for jj = col_bits
        row_count = 1;
        
        for ii = row_bits
            
            my_bit = slice(ii:(ii+A-1),jj:(jj+B-1));
    
            temp_grid(row_count,col_count) = mean(my_bit(:));
            
            row_count = row_count+1;
            
        end
        
        col_count = col_count+1;
        
    end
    
    grid_new{kk} = temp_grid;
    
end

grid_c = NaN*ones(length(row_bits),length(col_bits),c);

for kk = 1:c
    grid_c(:,:,kk) = grid_new{kk};
end
