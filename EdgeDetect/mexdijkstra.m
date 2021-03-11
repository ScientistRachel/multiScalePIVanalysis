% This function plugs the appropriate values from im_edge into the .mex
% functions 'Kadjacency' and 'dijkstra' and formats the outputs correctly
% for its parent function sheet_edger

function coords = mexdijkstra(im_edge,lr)

% INPUT VARIABLES
% ~ 'im_edge' is a logical image with 1's defining the edge of cells on the
%       sheet. The edge runs from top to bottom
% ~ 'lr' is a string indicated what side of the image is not covered with
%    cells (valid values: 'left','right','top','bottom')
%
% OUTPUT VARIABLES
% ~ 'coords' a coordinate sequence that spans the length of the edge
% File created by Peter Kordell
% 2019/08/27 RML allow for vertical or horizontal edges by including lr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% find the begining and ending nodes for the edge
[lx, ly] = find(im_edge);

if strcmp(lr,'left') || strcmp(lr,'right')
    s = find(max(lx)==lx);
    s = s(1);
    e = find(min(lx)==lx);
    e = e(1);
    
    % plug values into .mex functions
    A = Kadjacency([lx, ly]' , 2);
    path = dijkstra(A , s , e);

    % format and return values
    coords = [lx(path),ly(path)];    
    
elseif strcmp(lr,'top') || strcmp(lr,'bottom')
    
%     disp(lr)
    
    s = find(max(ly)==ly);
    s = s(1);
    e = find(min(ly)==ly);
    e = e(1);
    
    % plug values into .mex functions
    A = Kadjacency([lx, ly]' , 2);
    path = dijkstra(A , s , e);

    % format and return values
    coords = [lx(path),ly(path)];    
    
else
    error('lr not recognized: which side is the monolayer?')
end

