%This function reformats a logical input image so that the shortest path
%from the bottom to the top of the image can be computed by the djikstra
%algorithm containted in djikstra_mat

function coords = matdijkstra(im,lr)

% INPUT VARIABLES
% ~ 'im' is a logical matrix that has at least one possible path of true
%    values from the bottom of the image to the top
% ~ 'lr' is a string indicated what side of the image is not covered with
%    cells (valid values: 'left','right','top','bottom')

% OUTPUT VARIABLES
% ~ 'coords' is an N by 2 matrix that contains the pixel coordinates of the
%    shortest edge along the possible paths from input 'im'

% File created by Peter Kordell
% Edited by Rachel Lee, May 2011
% 2013/07/25 RML - added if statement to skip especially difficult frames
% 2019/08/27 RML - added ability to have horizontal or vertical edges by
% adding lr as an input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%convert logical trues in the input image into a point sequance with 
%trailing unique ID (found in the first colomn)---
[row,col]=find(im==1);
numnodes=size(row,1);
nodes=zeros(numnodes,3);
for i=1:size(row,1)
    nodes(i,1)=i;
    nodes(i,2)=row(i,1);
    nodes(i,3)=col(i,1);
end

%calculate the matrix 'segements' that stores row-wize the unique ID's of
%points that connect

%initialize and prealocate memory
segments = zeros(numnodes,3);
count = 0;    
    
for i=1:numnodes
    ix=nodes(i,2);
    iy=nodes(i,3);
    for j=i+1:numnodes
        jx=nodes(j,2);
        jy=nodes(j,3);
        if (abs(jx-ix)<=1 && abs(jy-iy)<=1)
            count = count + 1;
            segments(count,:) = [count, nodes(i,1), nodes(j,1)];
        end
    end
end

segments = segments(1:count,:);

%locate the unique ID's of the start and end nodes of the sequence
n1 = nodes(:,1);
n2 = nodes(:,2);
n3 = nodes(:,3);

if strcmp(lr,'left') || strcmp(lr,'right')
    start_id = n1.*((n3.*(n2 == max(n2))) == max(n3.*(n2 == max(n2))));
    start_id(start_id == 0) = [];
    finish_id = n1.*((n3.*(n2 == min(n2))) == max(n3.*(n2 == min(n2))));
    finish_id(finish_id == 0) = [];
elseif strcmp(lr,'top') || strcmp(lr,'bottom')
    start_id = n1.*((n2.*(n3 == max(n3))) == max(n2.*(n3 == max(n3))));
    start_id(start_id == 0) = [];
    finish_id = n1.*((n2.*(n3 == min(n3))) == max(n2.*(n3 == min(n3))));
    finish_id(finish_id == 0) = [];    
end

% feed matrix nodes, segements and the ID's of the start and
% end nodes into dijkstra_mat
[~, path] = dijkstra_mat(nodes,segments,start_id,finish_id);
if ~sum(isnan(path))
    % get x y sequance from point ID's given in path from dijkstra_mat (see function dijkstra_mat)
    coords(:,1) = nodes(path,2);
    coords(:,2) = nodes(path,3);
else
    coords = [];
end