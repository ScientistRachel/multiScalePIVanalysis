% Function to provide some automated analysis of PIV data.
% Please do not use this function without thinking!
% - Watch the original movie
% - Check the hard coded parameters
% - Watch that the results seem consistent!
%
% Usage: [ang_dev, mean_speed, SEM_speed, deltaR, firstframe, lastframe] = piv_analytics(directory, imname, r_scale, t_scale, firstframe, lastframe,pos_list)

function piv_analytics_timeMultiC(directory, imname, r_scale, t_scale, firstframe, lastframe,pos_list)
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
% Note: a figure and .mat file are also saved to the directory
% Figures: speed vs time, ang_dev vs time
% .mat Files: time_analytics
%
%%% SUBFUNCTIONS %%%
% bfopenForPIV.m, calc_ang_dev.m, circfit.m, cutter2.m,
% do_edgedat_to_circle.m, do_scale_shift_piv.m, get_zvi_location.m,
% rescale_image.m
%
%%% CHANGE LOG %%%
% Created by Rachel Lee, 2014/07/03
% based on piv_analytics
% 2015/12/14 edgedat file was loaded incorrectly!!! Fixed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%% ****PARAMETERS**** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These parameters are often constant between experiments, but should be
% confirmed for all data sets.
param.edge_cut = 0.75; % Only include vectors with r/R greater than this in edge metrics (set to 0 to analyze all)
param.center_cut = 0.5; % Only include vectors with r/R less than this in center metrics (set to 1 to analzye all)

param.size_limit = 500; % Number of vectors needed to keep a bin for r/R plots

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
if exist([directory 'Analytics Data\' imname '_dot_RvsT.mat'],'file')
    
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
% Shift PIV locations and set to appropriate scale
list = [dir([directory imname 'C*dot_piv.mat']) ; dir([directory imname '*dot_piv_filtered.mat'])];
u = []; v = u;
x = cell(numel(list),1); y = x;
for kk = 1:numel(list)
%     load([directory list(kk).name])
%     u = [u dot_piv.u_fi];
%     v = [v dot_piv.v_fi];
    load([directory list(kk).name])
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

%%%%%%%%%%%%%%%%%%%%%%%%%%% Basic Speed Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
speed = sqrt(u.^2+v.^2);
r = sqrt(x.^2 + y.^2);

% Restrict to the correct time frame
speed = speed(:,:,firstframe:lastframe);

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Angles Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
v_r = (u.*x + v.*y)./r; % Dot product of velocity and r_hat (polar coordinates)
v_th = (-u.*y + v.*x)./r; % Dot product of velocity and theta_hat (polar coordinates)

angles = atan2(v_th,v_r); % Angle between v_r and v_th
angles = angles(:,:,firstframe:lastframe);

%%%%%%%%%%%%%%%%%%%%%%% Time Calculations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hour_vec= (firstframe:lastframe)*t_scale/60; %plot versus hours

% Preallocate storage
N = NaN*ones(3,length(hour_vec));
speed_mean = N; speed_std = N;
v_r_mean = N; v_r_std = N;
v_th_mean = N; v_th_std = N;
v_r_mean_abs = N; v_r_std_abs = N; 
v_th_mean_abs = N; v_th_std_abs = N;
ang_dev_t = N;

% Perform binning
for j = 1:length(hour_vec)
    
    speed_now = speed(:,:,j);
    ang_now = angles(:,:,j);
    r_now = r_scaled(:,:,j);
    v_r_now = v_r(:,:,j);
    v_th_now = v_th(:,:,j);
    
    %%%% All spatial locations
    N(1,j) = numel(find(~isnan(speed_now)));
    
    speed_mean(1,j) = mean(speed_now(~isnan(speed_now)));
    speed_std(1,j) = std(speed_now(~isnan(speed_now)));

    % Averages without absolute value -- do they cancel?
    v_r_mean(1,j) = mean(v_r_now(~isnan(v_r_now)));
    v_r_std(1,j) = std(v_r_now(~isnan(v_r_now)));

    v_th_mean(1,j) = mean(v_th_now(~isnan(v_th_now)));
    v_th_std(1,j) = std(v_th_now(~isnan(v_th_now)));

    % Averages with absolute value -- what fraction of the motion?
    v_th_mean_abs(1,j) = mean(abs(v_th_now(~isnan(v_th_now))));
    v_th_std_abs(1,j) = std(abs(v_th_now(~isnan(v_th_now))));

    v_r_mean_abs(1,j) = mean(abs(v_r_now(~isnan(v_r_now))));
    v_r_std_abs(1,j) = std(abs(v_r_now(~isnan(v_r_now))));
    
    ang_dev_t(1,j) = calc_ang_dev(ang_now(~isnan(ang_now)));
    
    %%%% Edge only
    speed_now2 = speed_now(r_now > param.edge_cut);
    ang_now2 = ang_now(r_now > param.edge_cut);
    v_r_now2 = v_r_now(r_now > param.edge_cut);
    v_th_now2 = v_th_now(r_now > param.edge_cut);
    
    N(2,j) = numel(find(~isnan(speed_now2)));
    
    speed_mean(2,j) = mean(speed_now2(~isnan(speed_now2)));
    speed_std(2,j) = std(speed_now2(~isnan(speed_now2)));

    % Averages without absolute value -- do they cancel?
    v_r_mean(2,j) = mean(v_r_now2(~isnan(v_r_now2)));
    v_r_std(2,j) = std(v_r_now2(~isnan(v_r_now2)));

    v_th_mean(2,j) = mean(v_th_now2(~isnan(v_th_now2)));
    v_th_std(2,j) = std(v_th_now2(~isnan(v_th_now2)));

    % Averages with absolute value -- what fraction of the motion?
    v_th_mean_abs(2,j) = mean(abs(v_th_now2(~isnan(v_th_now2))));
    v_th_std_abs(2,j) = std(abs(v_th_now2(~isnan(v_th_now2))));

    v_r_mean_abs(2,j) = mean(abs(v_r_now2(~isnan(v_r_now2))));
    v_r_std_abs(2,j) = std(abs(v_r_now2(~isnan(v_r_now2))));
    
    ang_dev_t(2,j) = calc_ang_dev(ang_now2(~isnan(ang_now2)));
    
    %%%% Center only
    speed_now2 = speed_now(r_now < param.center_cut);
    ang_now2 = ang_now(r_now < param.center_cut);
    v_r_now2 = v_r_now(r_now < param.center_cut);
    v_th_now2 = v_th_now(r_now < param.center_cut);
    
    N(3,j) = numel(find(~isnan(speed_now2)));
    
    speed_mean(3,j) = mean(speed_now2(~isnan(speed_now2)));
    speed_std(3,j) = std(speed_now2(~isnan(speed_now2)));

    % Averages without absolute value -- do they cancel?
    v_r_mean(3,j) = mean(v_r_now2(~isnan(v_r_now2)));
    v_r_std(3,j) = std(v_r_now2(~isnan(v_r_now2)));

    v_th_mean(3,j) = mean(v_th_now2(~isnan(v_th_now2)));
    v_th_std(3,j) = std(v_th_now2(~isnan(v_th_now2)));

    % Averages with absolute value -- what fraction of the motion?
    v_th_mean_abs(3,j) = mean(abs(v_th_now2(~isnan(v_th_now2))));
    v_th_std_abs(3,j) = std(abs(v_th_now2(~isnan(v_th_now2))));

    v_r_mean_abs(3,j) = mean(abs(v_r_now2(~isnan(v_r_now2))));
    v_r_std_abs(3,j) = std(abs(v_r_now2(~isnan(v_r_now2))));
    
    ang_dev_t(3,j) = calc_ang_dev(ang_now2(~isnan(ang_now2)));

end

% Save components
save([directory 'Analytics Data\' imname '_time_analytics.mat'],'speed_mean','speed_std',...
    'v_r_mean','v_r_std','v_th_mean','v_th_std','v_th_mean_abs','v_th_std_abs',...
    'v_r_mean_abs','v_r_std_abs','ang_dev_t','N','hour_vec','param')

%  Make Figures
figure(1)
plot(hour_vec,speed_mean(2,:),'ko','LineWidth',2)
hold on
plot(hour_vec,v_r_mean_abs(2,:),'ro','LineWidth',2)
plot(hour_vec,v_th_mean_abs(2,:),'bo','LineWidth',2)
hold off
xlabel('Hour','FontSize',16)
ylabel('Speed (\mum/min)','FontSize',16)
legend('Total','Radial','Rotational','Location','Northwest','Orientation','Horizontal')
title([imname 10 'Edge Speed vs Time'],'FontSize',16,'Interpreter','none')
set(gca,'FontSize',16)
xlim([min(hour_vec) max(hour_vec)])
% Save figure
savename = [directory 'Analytics Figures\' imname 'Radial Speed vs Time'];
saveas(gcf,[savename '.fig'],'fig')
saveas(gcf,[savename '.tif'],'tif')

figure(2)
plot(hour_vec,ang_dev_t(2,:),'ro','LineWidth',2)
hold on
plot(hour_vec,ang_dev_t(3,:),'bo','LineWidth',2)
plot([min(hour_vec) max(hour_vec)],sqrt(2)*[1 1],'--k')
hold off
xlabel('Hour','FontSize',16)
ylabel('Angular Deviation','FontSize',16)
legend('Edge','Center','Upper Bound','Location','Southwest','Orientation','Horizontal')
title([imname 10 'Angular Deviation'],'FontSize',16,'Interpreter','none')
set(gca,'FontSize',16)
xlim([min(hour_vec) max(hour_vec)])
ylim([0 1.5])
% Save figure
savename = [directory 'Analytics Figures\' imname 'Angular Deviation vs Time'];
saveas(gcf,[savename '.fig'],'fig')
saveas(gcf,[savename '.tif'],'tif')

end %This ends the loops that skips everything if files are missing