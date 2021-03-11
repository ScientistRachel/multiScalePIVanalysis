%new method used to eliminate storage zeros from arrays and matrices 

function points_seq_slice = cutter2(points_seq_slice)

% INPUT VARIABLES
% ~ 'points_seq_slice' is the array/matrix which will have storage zeros
%    eliminated

% OUTPUT VARIABLES
% ~ 'points_seq_slice' is the array/matrix with the zeros removed

% File created by Peter Kordell
% Edited by Rachel Lee, May 2011

[row,~] = find(points_seq_slice,1,'last');
points_seq_slice = points_seq_slice(1:row,:);