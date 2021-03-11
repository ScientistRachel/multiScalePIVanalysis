% Function to provide some automated analysis of PIV data.
% Please do not use this function without thinking!
% - Watch the original movie
% - Check the hard coded parameters
% - Watch that the results seem consistent!
%
% Usage: [ang_dev, mean_speed, SEM_speed, deltaR, firstframe, lastframe] = piv_analytics(directory, imname, r_scale, t_scale, firstframe, lastframe,pos_list)

function plot_piv_time_hist(directory, imname, r_scale, t_scale, firstframe, lastframe,pos_list)
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
%%% SUBFUNCTIONS %%%
% bfopenForPIV.m, calc_ang_dev.m, circfit.m, cutter2.m,
% do_edgedat_to_circle.m, do_scale_shift_piv.m, get_zvi_location.m,
% plot_kymo.m, rescale_image.m, rose_normal.m
%
%%% CHANGE LOG %%%
% piv_analytics.m created by Rachel Lee, 2014/02/04
% 2014/02/12 RML    Add edge and center values to ang_dev and mean_speed
%                   Plot two rose plots -- total and edge
% 2014/03/30 RML    Fixed frame range error on kymograph plots
%                   Made type of position list an input
% 2014/04/02 RML    Used piv_analytics as staring point for code to plot
%                   directionality and speed versus time and r/R
% 2015/06/29 RML    Updated to be compatible with MultiC loading
%                   Changed colormap to parula
% 2015/08/15 RML    Used plot_piv_time_r as the base for plot_piv_time_hist
% 2015/12/15 edgedat file was loaded incorrectly!!! Fixed
% 2015/12/15 changed an isinf to an isnan
% 2016/06/21 Changed the parameter time_bin to always bin by the hour
% regarldess of t_scale (previously was hard coded to be 20)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%% ****PARAMETERS**** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These parameters are often constant between experiments, but should be
% confirmed for all data sets.
param.edge_cut = 0.75; % Only include vectors with r/R greater than this in edge metrics (set to 0 to analyze all)
param.center_cut = 0.5; % Only include vectors with r/R less than this in center metrics (set to 1 to analzye all)

param.speed_bins = 0:.1:3;
angle_bins = 20;
% param.ang_bins = (0:angle_bins-1)*2*pi/angle_bins - pi;
param.ang_bins = linspace(-pi,pi,angle_bins);
param.time_bin = floorR(60/t_scale); % Bin each hour

if exist('pos_list','var') && ~isempty(pos_list)
    param.pos_list = pos_list;
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
if ~exist([directory 'PIV Hist\'],'file')
    mkdir(directory,'PIV Hist')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% Circle Fitting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%%%%%%%%%%%%%%%%%%%%%%% Convert PIV to Radial %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Shift PIV locations and set to appropriate scale
list = [dir([directory imname 'C*dot_piv.mat']) ; dir([directory imname '*dot_piv_filtered.mat'])];
u = []; v = u;
x = cell(numel(list),1); y = x;
for kk = 1:numel(list)
    load([directory list(kk).name])
    u = [u dot_piv.u_fi];
    v = [v dot_piv.v_fi];
    x{kk} = dot_piv.x;
    y{kk} = dot_piv.y;
end
clear dot_piv

% x,y = um; u,v = um/min; r_scaled = unitless; beginning,ending = frames
[x,y,u,v,r_scaled, beginning, ending] = do_scale_shift_MultiC(u,v,x,y,x_pos,y_pos,X,Y,R,r_scale,t_scale);

% Update frame range to not exceed rational ranges
firstframe = max(firstframe,beginning);
lastframe = min(lastframe,ending);
r_scaled = r_scaled(:,:,firstframe:lastframe);

%%%%%%%%%%%%%%%%%%%%%%% Relevant Quantities %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Heatmap Kymograph: Speed
speed = sqrt(u.^2+v.^2);
r = sqrt(x.^2 + y.^2);

% Restrict to the correct time frame
speed = speed(:,:,firstframe:lastframe);

% Heatmap Kymograph: cos(theta)
v_r = (u.*x + v.*y)./r; % Dot product of velocity and r_hat (polar coordinates)
v_th = (-u.*y + v.*x)./r; % Dot product of velocity and theta_hat (polar coordinates)

angles = atan2(v_th,v_r); % Angle between v_r and v_th
angles = angles(:,:,firstframe:lastframe);

speed_edge = speed;
speed_edge(r_scaled < param.edge_cut) = NaN;

angles_edge = angles;
angles_edge(r_scaled < param.edge_cut) = NaN;


%%%%%%%%%%%%%%%%%%%%%%%%%%%% Histograms %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time_pieces = floor(size(speed_edge,3)/param.time_bin);
time_step = param.time_bin*t_scale/60;
time_vec = 1:time_step:time_step*time_pieces;
color_label = cell(1,length(time_vec));
for kk = 1:length(time_vec)
    color_label{kk} = num2str(time_vec(kk));
end
x_tick_places = 0:(1/time_pieces):((time_pieces-1)/time_pieces);
x_tick_places = x_tick_places + mean(diff(x_tick_places))/2;

% Speed
savename = [directory 'PIV Hist\' imname 'speed_edge_hist'];
figure(1)
% Take histogram and normalize
n = NaN*ones(length(param.speed_bins),time_pieces);
for kk = 1:time_pieces
    
    frames = (param.time_bin*(kk-1)+1):(param.time_bin*kk);
    
    slice = speed_edge(:,:,frames);
    slice = slice(~isnan(slice)); %2015/12/15 changed from isinf to isnan
    
    n(:,kk) = hist(slice,param.speed_bins);
    n(:,kk) = n(:,kk)/trapz(param.speed_bins,n(:,kk));
    
end
cmap = parula(size(n,2));
hold on
for kk = 1:size(n,2)
    plot(param.speed_bins,n(:,kk),'LineWidth',2,'Color',cmap(kk,:))
end
xlabel('Speed (\mum/min)','FontSize',16)
ylabel('Probability Density','FontSize',16)
set(gca,'FontSize',16)
xlim([min(param.speed_bins) max(param.speed_bins)])
colormap(cmap)
h = colorbar;
set(h,'Ticks',x_tick_places,'TickLabels',color_label,'TickLength',0)
set(get(h,'Label'),'String','Hour','FontSize',16)
saveas(gcf,[savename '.tif'],'tif')

% Angles
savename = [directory 'PIV Hist\' imname 'angle_edge_hist'];
figure(2)
% Take histogram and normalize
n = NaN*ones(length(param.ang_bins),time_pieces);
for kk = 1:time_pieces
    
    frames = (param.time_bin*(kk-1)+1):(param.time_bin*kk);
    
    slice = angles_edge(:,:,frames);
    slice = slice(~isinf(slice));
    
    n(:,kk) = hist(slice,param.ang_bins);
    n(:,kk) = n(:,kk)/trapz(param.ang_bins,n(:,kk));
    
end
cmap = parula(size(n,2));
hold on
for kk = 1:size(n,2)
    plot(param.ang_bins,n(:,kk),'LineWidth',2,'Color',cmap(kk,:))
end
xlabel('Velocity Direction','FontSize',16)
ylabel('Probability Density','FontSize',16)
set(gca,'FontSize',16)
xlim([min(param.ang_bins) max(param.ang_bins)])
colormap(cmap)
h = colorbar;
set(h,'Ticks',x_tick_places,'TickLabels',color_label,'TickLength',0)
set(get(h,'Label'),'String','Hour','FontSize',16)
saveas(gcf,[savename '.tif'],'tif')

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Save Outputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save([directory 'PIV Hist\' imname 'PIV_edge.mat'],'speed','angles'...
    ,'speed_edge', 'angles_edge','firstframe','lastframe','param')