% This funcation calculates the angular deviation of a set of angles.  The
% angular deviation varies from zero to sqrt(2).
%
% Usage: ang_dev = calc_ang_dev(angles)

function ang_dev = calc_ang_dev(angles)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS
% angles  = set of angles to analzye
% OUTPUTS
% ang_dev = angular deviation
%
% Created by Rachel Lee, 2014/02/05
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

term1 = sum(cos(angles(~isnan(angles))))/numel(find(~isnan(angles)));
term2 = sum(sin(angles(~isnan(angles))))/numel(find(~isnan(angles)));
R = sqrt(term1.^2 + term2.^2);
ang_dev = sqrt(2*(1-R));
