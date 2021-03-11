% Function to provide some automated analysis of PIV data.
% Please do not use this function without thinking!
% - Watch the original movie
% - Check the hard coded parameters
% - Watch that the results seem consistent!
%
% Usage: [ang_dev, mean_speed, SEM_speed, deltaR, firstframe, lastframe] = piv_analytics(directory, imname, r_scale, t_scale, firstframe, lastframe,pos_list)

function piv_corr_r(directory, imname, r_scale, t_scale, firstframe, lastframe,pos_list)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% INPUTS %%%
% directory     = folder of images
% imname        = file names (missing L, R, or C at the end)
% r_scale       = magnification of images in microns per pixel
% t_scale       = time scale of movie, in minutes per frame
% firstframe    = Optional. First frame to analyze. Defaults to 1.
% lastframe     = Optional. Last frame to analyze. Defaults to full image series.
% pos_list      = Specify how to calculation image locations.  Defaults to 1.
%                   0 = Use zvi files to extract position (robust but slow)
%                   1 = Use 'PosList.csv' to find postions (must be in same directory)
%                   2 = 'dot_RvsT.mat' file already exists, use previous analysis
%
% Note: frame range will be limited by edge detection results.
%
%%% OUTPUTS %%%
% Several figures and .mat files are also saved to the directory
% THIS LIST IS WRONG (FROM OLD FILE)
% Figures: 2 heatmaps, 3 rose plots, R vs time, speed and vel vs r/R
% .mat Files: dotRvsT, velcomp, Analytics
%
%%% CHANGE LOG %%%
% 2015/07/09 RML - created piv_corr_t based on piv_analyticsMultiC
% 2015/12/14 edgedat file was loaded incorrectly!!! Fixed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%% ****PARAMETERS**** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These parameters are often constant between experiments, but should be
% confirmed for all data sets.
param.edge_cut = 0.75; % Only include vectors with r/R greater than this in edge metrics (set to 0 to analyze all)
param.center_cut = 0.5; % Only include vectors with r/R less than this in center metrics (set to 1 to analzye all)
max_dist = 400; % in um, put inf to use as much as possible based on image size

if exist('pos_list','var') && ~isempty(pos_list)
    param.pos_list = pos_list;
elseif exist([directory 'Analytics Data\' imname '_dot_RvsT.mat'],'file')
    param.pos_list = 2;
else
    param.pos_list = 1; % 1 = excel file exists with locations, 0 = find from zvi files (slower)
end

if ~exist('firstframe','var') || isempty(firstframe)
    firstframe = 1;
end

if ~exist('lastframe','var') || isempty(lastframe)
    lastframe = inf;
end

%%%%%%%%%%%%%%%%%%%%%%% Make Folder for Saving %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist([directory 'Correlation Data\'],'file')
    mkdir(directory,'Correlation Data')
end

if ~exist([directory 'Correlation Data\'],'file')
    mkdir(directory,'Correlation Data')
end
    

%%%%%%%%%%%%%%%%%%%%%%%%%%% Circle Fitting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if param.pos_list == 2
    
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

%%%%%%%%%%%%%%%%%%%%%%% Convert PIV to Radial %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Shift PIV locations and set to appropriate scale
list = [dir([directory imname 'C*dot_piv.mat']) ; dir([directory imname '*dot_piv_filtered.mat'])];
u_raw = []; v_raw = u_raw;
x = cell(numel(list),1); y = x;
for kk = 1:numel(list)
    load([directory list(kk).name])
    u_raw = [u_raw dot_piv.u_fi];
    v_raw = [v_raw dot_piv.v_fi];
    x{kk} = dot_piv.x;
    y{kk} = dot_piv.y;
end
clear dot_piv

% x,y = um; u,v = um/min; r_scaled = unitless; beginning,ending = frames
[x,y,u,v,r_scaled, beginning, ending] = do_scale_shift_MultiC(u_raw,v_raw,x,y,x_pos,y_pos,X,Y,R,r_scale,t_scale);

% Update frame range to not exceed rational ranges
firstframe = max(firstframe,beginning);
lastframe = min(lastframe,ending);
r_scaled = r_scaled(:,:,firstframe:lastframe);

%%%%%%%%%%%%%%%%%%%%%%%%% Useful quantities %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u = u(:,:,firstframe:lastframe);
v = v(:,:,firstframe:lastframe);
x = x(:,:,firstframe:lastframe);
y = y(:,:,firstframe:lastframe);

% PIVtopixel = x(1,2) - x(1,1);

speed = sqrt(u.^2+v.^2);
r = sqrt(x.^2 + y.^2);

v_r = (u.*x + v.*y)./r; % Dot product of velocity and r_hat (polar coordinates)
v_th = (-u.*y + v.*x)./r; % Dot product of velocity and theta_hat (polar coordinates)

% angles = atan2(v_th,v_r); % Angle between v_r and v_th

% Keep only the edge data
speed(r_scaled < param.edge_cut) = NaN;
v_r(r_scaled < param.edge_cut) = NaN;
v_th(r_scaled < param.edge_cut) = NaN;
x(r_scaled < param.edge_cut) = NaN;
y(r_scaled < param.edge_cut) = NaN;

%%%%%%%%%%%%%%%%%%%%%%%%% Time correlations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(1)
savedir = [directory 'Correlation Data\' imname '_speed_radial_corr'];
correlatePIV_r(speed,x,y,t_scale,r_scale,0,lastframe,max_dist,savedir,'Autocorrelation of speed')

% close all

% figure(2)
savedir = [directory 'Correlation Data\' imname '_vr_radial_corr'];
correlatePIV_r(v_r,x,y,t_scale,r_scale,0,lastframe,max_dist,savedir,'Autocorrelation of v_r')

% close all

% figure(3)
savedir = [directory 'Correlation Data\' imname '_vth_radial_corr'];
correlatePIV_r(v_th,x,y,t_scale,r_scale,0,lastframe,max_dist,savedir,'Autocorrelation of v_\theta')

% close all

