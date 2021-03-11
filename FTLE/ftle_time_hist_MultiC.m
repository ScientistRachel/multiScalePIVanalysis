% Function to provide some automated analysis of FTLE data.
% Please do not use this function without thinking!
% - Watch the original movie
% - Check the hard coded parameters
% - Watch that the results seem consistent!
%
% Usage: ftle_time_hist_MultiC(directory, imname, r_scale, t_scale, firstframe, lastframe)

function ftle_time_hist_MultiC(directory, imname, ftleLRC, r_scale, t_scale, firstframe, lastframe)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% INPUTS %%%
% directory     = folder of images
% imname        = file names (missing L, R, or C at the end)
% r_scale       = magnification of images in microns per pixel
% t_scale       = time scale of movie, in minutes per frame
% firstframe    = Optional. First frame to analyze. Defaults to 1.
% lastframe     = Optional. Last frame to analyze. Defaults to full image series.
%
% Note: frame range will be limited by edge detection results.
%
%%% OUTPUTS %%%
%
% Note: several figures and .mat files are also saved to the directory
% Figures:
% .mat Files: 
%
%%% CHANGE LOG %%%
% Created by Rachel Lee, 2014/02/19, based on piv_analytics
% 2014/03/27 RML Added FTLE heatmap, Added FTLE vs time
% 2015/06/15 used ftle_analytics as a base for time histograms
% 2015/12/14 edgedat file was loaded incorrectly!!! Fixed
% 2016/06/21 Changed the parameter time_bin to always bin by the hour
% regarldess of t_scale (previously was hard coded to be 20)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%% ****PARAMETERS**** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These parameters are often constant between experiments, but should be
% confirmed for all data sets.
param.edge_cut = 0.75; % Only include vectors with r/R greater than this in edge metrics (set to 0 to analyze all)
param.center_cut = 0.5; % Only include vectors with r/R less than this in center metrics (set to 1 to analzye all)
param.time_bin = floorR(60/t_scale); % Bin each hour

param.ftle_bins = -3.5:.05:1.5; % Bins (in hr^-1) for FTLE histograms
param.pos_list = 1; % 1 = excel file exists with locations, 0 = find from zvi files (slower)

if ~exist('firstframe','var') || isempty(firstframe)
    firstframe = 1;
end

if ~exist('lastframe','var') || isempty(lastframe)
    lastframe = inf;
end

%%%%%%%%%%%%%%%%%%%%%%% Make Folder for Saving %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist([directory 'FTLE Analytics\'],'file')
    mkdir(directory,'FTLE Analytics')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% Circle Fitting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %If piv_analytics was already run, save time by loading that radius data
if exist([directory 'Analytics Data\' imname '_dot_RvsT.mat'],'file')
    load([directory 'Analytics Data\' imname '_dot_RvsT.mat'])
else

    edgelist = dir([directory '*' imname '*edgedat.mat']);
    DotE = cell(2,1);
    for kk = 1:2
        load([directory edgelist(kk).name])
        DotE{kk} = edgedat;
        clear edgedat
    end

    % All outputs in um
    [X,Y,R,x_pos,y_pos] = do_edgedat_to_circleMultiC(DotE{1}, DotE{2}, pos_list, directory, imname,r_scale);
    clear DotE edgelist
end

%%%%%%%%%%%%%%%%%%%%%%% Convert FTLE to Radial %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Shift PIV locations and set to appropriate scale
% Function assumes that all three PIV files (L, R, C), have the same
% positions, so only loads one file to get x and y values
load([directory imname 'Cdot_piv.mat'])
% Other options for files:
% load([directory imname 'Ldot_piv_filtered.mat'])
% load([directory imname 'Rdot_piv_filtered.mat'])

% ftles = hr^-1; r_scaled = unitless; beginning,ending = frames
[ftles,r_scaled,beginning, ending] = do_scale_shift_ftleMultiC(ftleLRC,dot_piv,x_pos,y_pos,X,Y,R,r_scale,t_scale);

% Update frame range to not exceed rational ranges
firstframe = max(firstframe,beginning);
lastframe = min(lastframe,ending);
ftles = ftles(:,:,firstframe:lastframe);
r_scaled = r_scaled(:,:,firstframe:lastframe);

%%%%%%%%%%%%%%%%%%%%%% Set up different regions %%%%%%%%%%%%%%%%%%%%%%%%%%%
ftle_edge = ftles;
ftle_edge(r_scaled < param.edge_cut) = inf;

ftle_center = ftles;
ftle_center(r_scaled > param.center_cut) = inf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Calc Means %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mean_ftle(1) = mean(ftles(~isnan(ftles)));
% SEM_ftle(1) = std(ftles(~isnan(ftles)))/numel(find(~isnan(ftles)));
per_pos(1) = numel(find(ftles > 0))/numel(find(~isnan(ftles)))*100;

ftle_slice = ftles(r_scaled > param.edge_cut);
% mean_ftle(2) = mean(ftle_slice(~isnan(ftle_slice)));
% SEM_ftle(2) = std(ftle_slice(~isnan(ftle_slice)))/numel(find(~isnan(ftle_slice)));
per_pos(2) = numel(find(ftle_slice > 0))/numel(find(~isnan(ftle_slice)))*100;

ftle_slice = ftles(r_scaled < param.center_cut);
% mean_ftle(3) = mean(ftle_slice(~isnan(ftle_slice)));
% SEM_ftle(3) = std(ftle_slice(~isnan(ftle_slice)))/numel(find(~isnan(ftle_slice)));
per_pos(3) = numel(find(ftle_slice > 0))/numel(find(~isnan(ftle_slice)))*100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Histograms %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All vectors
figure(1)
savename = [directory 'FTLE Analytics\' imname 'hist_time_total'];
title_text = ['Total Positive: ' num2str(per_pos(1),'%0.0f') ' %'];
plot_ftle_pdf_time(ftles,param.ftle_bins,param.time_bin,t_scale,savename,title_text)


% Edge vectors only
figure(2)
savename = [directory 'FTLE Analytics\' imname 'hist_time_edge'];
title_text = ['Edge Positive: ' num2str(per_pos(2),'%0.0f') ' %'];
plot_ftle_pdf_time(ftle_edge,param.ftle_bins,param.time_bin,t_scale,savename,title_text)

% Center vectors only
figure(3)
savename = [directory 'FTLE Analytics\' imname 'hist_time_center'];
title_text = ['Center Positive: ' num2str(per_pos(3),'%0.0f') ' %'];
plot_ftle_pdf_time(ftle_center,param.ftle_bins,param.time_bin,t_scale,savename,title_text)


%%%%%%%%%%%%%%%%%%%%%%%%%%%% Save Outputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save([directory 'FTLE Analytics\' imname 'FTLE_edge.mat'],'ftles'...
    ,'ftle_edge', 'ftle_center','firstframe','lastframe','param')
