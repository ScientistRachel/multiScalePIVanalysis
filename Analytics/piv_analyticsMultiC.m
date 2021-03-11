% Function to provide some automated analysis of PIV data.
% Please do not use this function without thinking!
% - Watch the original movie
% - Check the hard coded parameters
% - Watch that the results seem consistent!
%
% Usage: [ang_dev, mean_speed, SEM_speed, deltaR, firstframe, lastframe] = piv_analytics(directory, imname, r_scale, t_scale, firstframe, lastframe,pos_list)

function [ang_dev, mean_speed, SEM_speed, deltaR, firstframe, lastframe] = piv_analyticsMultiC(directory, imname, r_scale, t_scale, firstframe, lastframe,pos_list)
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
% ang_dev       = angular deviation ( 0 <= ang_dev <= sqrt(2) )
%               3x1 vector for all vectors, edge vectors, center vectors
% mean_speed    = mean speed over the entire movie (um/min)
%               3x1 vector for all vectors, edge vectors, center vectors
% SEM_speed     = standard error of the mean of the speed (um/min)
%               3x1 vector for all vectors, edge vectors, center vectors
% deltaR        = change in dot radius from first to last frame (um)
% firstframe    = first frame actually analyzed
% lastframe     = last frame actually analyzed
%
% Note: several figures and .mat files are also saved to the directory
% Figures: 2 heatmaps, 3 rose plots, R vs time, speed and vel vs r/R
% .mat Files: dotRvsT, velcomp, Analytics
%
%%% SUBFUNCTIONS %%%
% bfopenForPIV.m, calc_ang_dev.m, circfit.m, cutter2.m,
% do_edgedat_to_circle.m, do_scale_shift_piv.m, get_zvi_location.m,
% plot_kymo.m, rescale_image.m, rose_normal2.m
%
%%% CHANGE LOG %%%
% Created by Rachel Lee, 2014/02/04
% 2014/02/12 RML    Add edge and center values to ang_dev and mean_speed
%                   Plot two rose plots -- total and edge
% 2014/03/30 RML    Fixed frame range error on kymograph plots
%                   Made type of position list an input
% 2015/12/14 edgedat file was loaded incorrectly!!! Fixed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%% ****PARAMETERS**** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These parameters are often constant between experiments, but should be
% confirmed for all data sets.
param.edge_cut = 0.75; % Only include vectors with r/R greater than this in edge metrics (set to 0 to analyze all)
param.center_cut = 0.5; % Only include vectors with r/R less than this in center metrics (set to 1 to analzye all)

param.r_step = 50; % microns for non-scaled bins
param.numel_limit = 42; % Number of vectors needed to keep a bin for heatmap
param.speed_limits = [0 .7]; % Range for speed kymograph colormap, in um/min
param.cos_limits = [-1 1]; % Range for cos(theta) kymograph colormap -- recommend [0 1] or [-1 1]

param.angle_bins = 20; % Number of bins for rose plot
param.rose_limit = .25; % Set axis limit for rose plot

param.R_bin_size = 0.05; % bin metrics in units of r/R
param.R_bin_max = 1; % Maximum value of r/R to plot
param.size_limit = 1100; % Number of vectors needed to keep a bin for r/R plots

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
if ~exist([directory 'Analytics Data\'],'file')
    mkdir(directory,'Analytics Data')
end

if ~exist([directory 'Analytics Figures\'],'file')
    mkdir(directory,'Analytics Figures')
end

%%%%%%%%%%%%%%%%%%%%%%%% Validate User Input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function will assume that edgedat.mat and dot_piv.mat already exist
% File names need to have the last characters: L, R, C
% Check for the existence of relevant files -- some might not exist due to
% experimental errors, etc --> output NaN.
% bad_file = 0; %assume it is a good file for starters...
% if ~exist([directory imname 'L_edgedat.mat'],'file')
%     warning(['The file ' directory imname 'L_edgedat.mat does not exist.' 10 ...
%         'All analytics output for this well will equal NaN'])
%     bad_file = 1;
% end
% 
% if ~exist([directory imname 'R_edgedat.mat'],'file')
%     warning(['The file ' directory imname 'R_edgedat.mat does not exist.' 10 ...
%         'All analytics output for this well will equal NaN'])
%     bad_file = 1;
% end
% Be flexible on PIV parts for now, otherwise uncomment below:
% if ~exist([directory imname 'Cdot_piv.mat'],'file')
%     warning(['The file ' directory imname 'Cdot_piv.mat does not exist.' 10 ...
%         'All analytics output for this well will equal NaN'])
%     bad_file = 1;
% end
% 
% if ~exist([directory imname 'Ldot_piv_filtered.mat'],'file')
%     warning(['The file ' directory imname 'Ldot_piv_filtered.mat does not exist.' 10 ...
%         'All analytics output for this well will equal NaN'])
%     bad_file = 1;
% end
% 
% if ~exist([directory imname 'Rdot_piv_filtered.mat'],'file')
%     warning(['The file ' directory imname 'Rdot_piv_filtered.mat does not exist.' 10 ...
%         'All analytics output for this well will equal NaN'])
%     bad_file = 1;
% end
% 
% if bad_file
%     ang_dev = NaN*[1 1 1];
%     mean_speed = NaN*[1 1 1];
%     SEM_speed = NaN*[1 1 1];
%     deltaR = NaN;
%     firstframe = NaN;
%     lastframe = NaN;
%     
% else
    

%%%%%%%%%%%%%%%%%%%%%%%%%%% Circle Fitting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if param.pos_list == 2
    
    load([directory 'Analytics Data\' imname '_dot_RvsT.mat'])
    
else

%     load([directory imname 'L_edgedat.mat'])
%     DotL = edgedat;
%     clear edgedat
% 
%     load([directory imname 'R_edgedat.mat'])
%     DotR = edgedat;
%     clear edgedat
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
%     disp(list(kk).name)
    c = size(dot_piv.u_fi,3); %2016/03/02 In case PIV is differently sized...
%     disp(c)
    if size(u_raw,3) > c
        u_raw = u_raw(:,:,1:c);
        v_raw = v_raw(:,:,1:c);
    end
    if kk == 1
        b = c;
    else
        b = size(u_raw,3);
    end
    u_raw = [u_raw dot_piv.u_fi(:,:,1:b)];
    v_raw = [v_raw dot_piv.v_fi(:,:,1:b)];
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

%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot Heatmaps %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Heatmap Kymograph: Speed
speed = sqrt(u.^2+v.^2);
r = sqrt(x.^2 + y.^2);

speed_slice = speed(:,:,firstframe:lastframe);

title_text = ['Mean Speed: ' num2str(mean(speed_slice(~isnan(speed_slice))),'%0.03f') ' \mum/min'];
bar_text = '\mum/min';
savename = [directory 'Analytics Figures\' imname 'speed_kymo'];

plot_kymo(speed_slice,r(:,:,1),t_scale,param.r_step,param.numel_limit,1,lastframe-firstframe+1,param.speed_limits,title_text,savename,bar_text,imname)

% Restrict to the correct time frame
speed = speed(:,:,firstframe:lastframe);

% Calculate means
mean_speed(1) = mean(speed(~isnan(speed)));
SEM_speed(1) = std(speed(~isnan(speed)))/sqrt(numel(find(~isnan(speed))));

speed_temp = speed(r_scaled > param.edge_cut);
mean_speed(2) = mean(speed_temp(~isnan(speed_temp)));
SEM_speed(2) = std(speed_temp(~isnan(speed_temp)))/sqrt(numel(find(~isnan(speed_temp))));

speed_temp = speed(r_scaled < param.center_cut);
mean_speed(3) = mean(speed_temp(~isnan(speed_temp)));
SEM_speed(3) = std(speed_temp(~isnan(speed_temp)))/sqrt(numel(find(~isnan(speed_temp))));


% Heatmap Kymograph: cos(theta)
v_r = (u.*x + v.*y)./r; % Dot product of velocity and r_hat (polar coordinates)
v_th = (-u.*y + v.*x)./r; % Dot product of velocity and theta_hat (polar coordinates)

angles = atan2(v_th,v_r); % Angle between v_r and v_th

ang_dev(1) = calc_ang_dev(angles(:,:,firstframe:lastframe));
title_text = ['Total Angular Deviation: ' num2str(ang_dev(1),'%0.03f')];
bar_text = 'cos(\theta)';
savename = [directory 'Analytics Figures\' imname 'cos_kymo'];

plot_kymo(cos(angles),r(:,:,1),t_scale,param.r_step,param.numel_limit,1,lastframe-firstframe+1,param.cos_limits,title_text,savename,bar_text,imname)

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Rose Plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
angles = angles(:,:,firstframe:lastframe);
% Rose plot is normalized by the total count
% All vectors
savename = [directory 'Analytics Figures\' imname 'rose_total'];
title_text = ['Total Angular Deviation: ' num2str(ang_dev(1),'%0.03f')];
rose_normal2(angles,param.angle_bins,param.rose_limit,savename,[],title_text,imname)

% Edge vectors only
ang_dev(2) = calc_ang_dev(angles(r_scaled > param.edge_cut));
title_text = ['Edge Angular Deviation: ' num2str(ang_dev(2),'%0.03f')];
savename = [directory 'Analytics Figures\' imname 'rose_edge'];
rose_normal2(angles(r_scaled > param.edge_cut),param.angle_bins,param.rose_limit,savename,[],title_text,imname)

% Center vectors only
ang_dev(3) = calc_ang_dev(angles(r_scaled < param.center_cut));
title_text = ['Center Angular Deviation: ' num2str(ang_dev(3),'%0.03f')];
savename = [directory 'Analytics Figures\' imname 'rose_center'];
rose_normal2(angles(r_scaled < param.center_cut),param.angle_bins,param.rose_limit,savename,[],title_text,imname)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Dot Size %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R_smooth = smooth(R,'lowess');
deltaR = R_smooth(lastframe) - R_smooth(firstframe);
% Plot R vs t
time_vec = (firstframe:lastframe)*3/60; % hours
plot(time_vec,R(firstframe:lastframe)/1000,'k','LineWidth',2)
xlabel('Time (hours)','FontSize',16)
ylabel('Dot Radius (mm)','FontSize',16)
title(imname,'FontSize',16,'Interpreter','none')
set(gca,'FontSize',16)
xlim([time_vec(1) time_vec(end)])

savename = [directory 'Analytics Figures\' imname 'size'];
saveas(gcf,[savename '.fig'],'fig')
saveas(gcf,[savename '.tif'],'tif')

%%%%%%%%%%%%%%%%%%%%%%% Velocity Components %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bins = 0:param.R_bin_size:param.R_bin_max;

% Preallocate storage
bin_mid = NaN*ones(1,length(bins) - 1);
N = bin_mid; r_mean = bin_mid; r_std = bin_mid; speed_mean = bin_mid;
speed_std = bin_mid; v_r_mean = bin_mid; v_r_std = bin_mid; v_th_mean = bin_mid;
v_th_std = bin_mid; v_r_mean_abs = bin_mid; v_r_std_abs = bin_mid; 
v_th_mean_abs = bin_mid; v_th_std_abs = bin_mid;

% Perform binning
for j = 1:(length(bins) - 1)

    now1 = r_scaled > bins(j);
    now2 = r_scaled < bins(j+1);
    now3 = now1.*now2;
    now = find(now3);
    
    bin_mid(j) = (bins(j)+bins(j+1))/2;
    
    speed_now = speed(now);
    
    N(j) = numel(find(~isnan(speed_now)));
    
    if N(j) > param.size_limit
    
        r_mean(j) = mean(r_scaled(now));
        r_std(j) = std(r_scaled(now));
       
        v_r_now = v_r(now);
        v_th_now = v_th(now);

        speed_mean(j) = mean(speed_now(~isnan(speed_now)));
        speed_std(j) = std(speed_now(~isnan(speed_now)));

        % Averages without absolute value -- do they cancel?
        v_r_mean(j) = mean(v_r_now(~isnan(v_r_now)));
        v_r_std(j) = std(v_r_now(~isnan(v_r_now)));

        v_th_mean(j) = mean(v_th_now(~isnan(v_th_now)));
        v_th_std(j) = std(v_th_now(~isnan(v_th_now)));

        % Averages with absolute value -- what fraction of the motion?
        v_th_mean_abs(j) = mean(abs(v_th_now(~isnan(v_th_now))));
        v_th_std_abs(j) = std(abs(v_th_now(~isnan(v_th_now))));

        v_r_mean_abs(j) = mean(abs(v_r_now(~isnan(v_r_now))));
        v_r_std_abs(j) = std(abs(v_r_now(~isnan(v_r_now))));
    
    end

end

% Save components
save([directory 'Analytics Data\' imname '_vel_comp.mat'],'speed_mean','speed_std',...
    'v_r_mean','v_r_std','v_th_mean','v_th_std','v_th_mean_abs','v_th_std_abs',...
    'v_r_mean_abs','v_r_std_abs','bin_mid','N','r_mean','r_std','firstframe','lastframe')

errorbar(bin_mid,speed_mean,speed_std./sqrt(N),'k-s','LineWidth',2,'MarkerFaceColor','k')
hold on
errorbar(bin_mid,v_r_mean_abs,v_r_std_abs./sqrt(N),'r-s','LineWidth',2,'MarkerFaceColor','r')
errorbar(bin_mid,v_th_mean_abs,v_th_std_abs./sqrt(N),'b-s','LineWidth',2,'MarkerFaceColor','b')
hold off
xlabel('r/R','FontSize',16)
ylabel('Speed (\mum/min)','FontSize',16)
legend('Total','Radial','Rotational','Location','Northwest','Orientation','Horizontal')
title([imname 10 'Speed Components'],'FontSize',16,'Interpreter','none')
set(gca,'FontSize',16)
xlim([0 max(bins)+param.R_bin_size])
% Save figure
savename = [directory 'Analytics Figures\' imname 'Radial Speed'];
saveas(gcf,[savename '.fig'],'fig')
saveas(gcf,[savename '.tif'],'tif')


errorbar(bin_mid,v_r_mean,v_r_std./sqrt(N),'r-s','LineWidth',2,'MarkerFaceColor','r')
hold on
errorbar(bin_mid,v_th_mean,v_th_std./sqrt(N),'b-s','LineWidth',2,'MarkerFaceColor','b')
plot([0 max(bins)],[0 0],'--k')
hold off
xlabel('r/R','FontSize',16)
ylabel('Velocity Component (\mum/min)','FontSize',16)
legend('Radial','Rotational','Location','Northwest','Orientation','Horizontal')
title([imname 10 'Velocity Components'],'FontSize',16,'Interpreter','none')
set(gca,'FontSize',16)
xlim([0 max(bins)+param.R_bin_size])
% Save figure
savename = [directory 'Analytics Figures\' imname 'Radial Velocity'];
saveas(gcf,[savename '.fig'],'fig')
saveas(gcf,[savename '.tif'],'tif')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Save Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save([directory 'Analytics Data\' imname 'Analytics.mat'],'ang_dev',...
    'mean_speed','SEM_speed','deltaR','firstframe','lastframe','param')

end %This ends the loops that skips everything if files are missing