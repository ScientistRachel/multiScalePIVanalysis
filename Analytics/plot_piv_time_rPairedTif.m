% Function to provide some automated analysis of PIV data.
% Please do not use this function without thinking!
% - Watch the original movie
% - Check the hard coded parameters
% - Watch that the results seem consistent!
%
% Usage: [ang_dev, mean_speed, SEM_speed, deltaR, firstframe, lastframe] = piv_analytics(directory, imname, r_scale, t_scale, firstframe, lastframe,pos_list)

function plot_piv_time_rPairedTif(directory, imnames, r_scale, t_scale, firstframe, lastframe,pos_list)
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
% 2015/12/14 edgedat file was loaded incorrectly!!! Fixed
% 2016/01/06 RML    Added v_r plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%% ****PARAMETERS**** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These parameters are often constant between experiments, but should be
% confirmed for all data sets.
param.edge_cut = 0.75; % Only include vectors with r/R greater than this in edge metrics (set to 0 to analyze all)
param.center_cut = 0.5; % Only include vectors with r/R less than this in center metrics (set to 1 to analzye all)

param.R_bin_size = 0.05; % bin metrics in units of r/R
param.R_bin_max = 1; % Maximum value of r/R to plot
param.size_limit = 5000; % Number of vectors needed to keep a bin for r/R plots
% param.time_bin = 20; %Number of frames per time bin
param.time_bin = 30; %Number of frames per time bin

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
if ~exist([directory 'PIV Time and R\'],'file')
    mkdir(directory,'PIV Time and R')
end

imname = [imnames{1} '_' imnames{2}];

%%%%%%%%%%%%%%%%%%%%%%%%%%% Circle Fitting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if exist([directory 'Analytics Data\' imname '_dot_RvsT.mat'],'file')
    
    load([directory 'Analytics Data\' imname '_dot_RvsT.mat'],'X','Y','R','x_pos','y_pos')
    
else
    load([directory imnames{1} '_edgedat.mat'],'edgedat')
    Dot1 = edgedat;
    clear edgedat

    load([directory imnames{2} '_edgedat.mat'],'edgedat')
    Dot2 = edgedat;
    clear edgedat

    % All outputs in um
    [X,Y,R,x_pos,y_pos] = do_edgedat_to_circlePairedtif(Dot1, Dot2, pos_list, directory, imnames,r_scale);
    clear Dot1 Dot2 edgelist
end

%%%%%%%%%%%%%%%%%%%%%%% Convert PIV to Radial %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Shift PIV locations and set to appropriate scale
% Shift PIV locations and set to appropriate scale
list = [dir([directory imnames{1} '*dot_piv_filtered.mat']) ; dir([directory imnames{2} '*dot_piv_filtered.mat'])];
u = []; v = u;
x = cell(numel(list),1); y = x;
for kk = 1:numel(list)
%     load([directory list(kk).name])
%     u = [u dot_piv.u_fi];
%     v = [v dot_piv.v_fi];
    load([directory list(kk).name],'dot_piv')
%     disp(list(kk).name)
    c = size(dot_piv.u_fi,3); %2016/03/02 In case PIV is differently sized...
%     disp(c)
    if size(u,3) > c
        u = u(:,:,1:c);
        v = v(:,:,1:c);
    end
    if kk == 1
        b = c;
    else
        b = size(u,3);
    end
    u = [u dot_piv.u_fi(:,:,1:b)];
    v = [v dot_piv.v_fi(:,:,1:b)];
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


%%%%%%%%%%%%%%%%%%%%%%% Velocity Components %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bins = 0:param.R_bin_size:param.R_bin_max;
time_pieces = floor(size(speed,3)/param.time_bin);

% Preallocate storage
bin_mid = NaN*ones(time_pieces,length(bins) - 1);
N = bin_mid;
r_mean = bin_mid; r_std = bin_mid;
speed_mean = bin_mid; speed_std = bin_mid;
cos_mean = bin_mid; cos_std = bin_mid;
v_r_mean = speed_mean; v_r_std = speed_mean;

% Perform binning
for kk = 1:time_pieces
    
    time_chunk = (param.time_bin*(kk-1)+1):(param.time_bin*kk);
    speed_piece = speed(:,:,time_chunk);
    ang_piece = angles(:,:,time_chunk);
    r_piece = r_scaled(:,:,time_chunk);
    v_r_piece = v_r(:,:,time_chunk);
    
    for j = 1:(length(bins) - 1)

        now1 = r_piece > bins(j);
        now2 = r_piece < bins(j+1);
        now3 = now1.*now2;
        now = find(now3);

        bin_mid(kk,j) = (bins(j)+bins(j+1))/2;

        speed_now = speed_piece(now);
        speed_now = speed_now(~isnan(speed_now));
        
        ang_now = ang_piece(now);
        ang_now = ang_now(~isnan(ang_now));
        
        v_r_now = v_r_piece(now);
        v_r_now = v_r_now(~isnan(v_r_now));
        
        N(kk,j) = numel(speed_now);

        if N(kk,j) > param.size_limit

            r_mean(kk,j) = mean(r_piece(now));
            r_std(kk,j) = std(r_piece(now));

            speed_mean(kk,j) = mean(speed_now);
            speed_std(kk,j) = std(speed_now);
            
            cos_mean(kk,j) = mean(cos(ang_now));
            cos_std(kk,j) = std(cos(ang_now));
            
            v_r_mean(kk,j) = mean(v_r_now);
            v_r_std(kk,j) = std(v_r_now);

        end

    end
end

%%%%%%%%%% Plotting
cmap = parula(time_pieces);
% Colormap labeling:
time_step = param.time_bin*t_scale/60;
time_vec = 1:time_step:time_step*time_pieces;
color_label = cell(1,length(time_vec));
for kk = 1:length(time_vec)
    color_label{kk} = num2str(time_vec(kk),'%0.1f');
end
x_tick_places = 0:(1/time_pieces):((time_pieces-1)/time_pieces);
x_tick_places = x_tick_places + mean(diff(x_tick_places))/2;

figure(1)
hold on
for kk = 1:time_pieces
    plot(bin_mid(kk,:),speed_mean(kk,:),'-','LineWidth',2,'Color',cmap(kk,:))
end
hold off
xlabel('r/R','FontSize',16)
ylabel('Speed (\mum/min)','FontSize',16)
title([imname 10 'Speed'],'FontSize',16,'Interpreter','none')
set(gca,'FontSize',16)
xlim([0 max(bins)+param.R_bin_size])
% a = get(gca,'YLim');
% disp(a)
% ylim([0 1])
colormap(cmap)
h = colorbar;
set(h,'Ticks',x_tick_places,'TickLabels',color_label,'TickLength',0)
set(get(h,'Label'),'String','Hour','FontSize',16)
% Save figure
savename = [directory 'PIV Time and R\' imname 'Speed'];
saveas(gcf,[savename '.fig'],'fig')
saveas(gcf,[savename '.tif'],'tif')

figure(2)
hold on
for kk = 1:time_pieces
    plot(bin_mid(kk,:),cos_mean(kk,:),'-','LineWidth',2,'Color',cmap(kk,:))
end
hold off
xlabel('r/R','FontSize',16)
ylabel('<cos \theta>','FontSize',16)
title([imname 10 '<cos \theta>'],'FontSize',16,'Interpreter','none')
set(gca,'FontSize',16)
xlim([0 max(bins)+param.R_bin_size])
ylim([-1 1])
colormap(cmap)
h = colorbar;
set(h,'Ticks',x_tick_places,'TickLabels',color_label,'TickLength',0)
set(get(h,'Label'),'String','Hour','FontSize',16)
% Save figure
savename = [directory 'PIV Time and R\' imname 'Directionality'];
saveas(gcf,[savename '.fig'],'fig')
saveas(gcf,[savename '.tif'],'tif')

figure(3)
hold on
for kk = 1:time_pieces
    plot(bin_mid(kk,:),v_r_mean(kk,:),'-','LineWidth',2,'Color',cmap(kk,:))
end
hold off
xlabel('r/R','FontSize',16)
ylabel('V_r (\mum/min)','FontSize',16)
title([imname 10 'Radial Velocity'],'FontSize',16,'Interpreter','none')
set(gca,'FontSize',16)
xlim([0 max(bins)+param.R_bin_size])
colormap(cmap)
h = colorbar;
set(h,'Ticks',x_tick_places,'TickLabels',color_label,'TickLength',0)
set(get(h,'Label'),'String','Hour','FontSize',16)
% Save figure
savename = [directory 'PIV Time and R\' imname 'Radial Velocity'];
saveas(gcf,[savename '.fig'],'fig')
saveas(gcf,[savename '.tif'],'tif')