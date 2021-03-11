% Function to provide some automated analysis of FTLE data.
% Please do not use this function without thinking!
% - Watch the original movie
% - Check the hard coded parameters
% - Watch that the results seem consistent!
%
% Usage: [mean_ftle, SEM_ftle, per_pos, firstframe, lastframe] = ftle_analytics(directory, imname, r_scale, t_scale, firstframe, lastframe)

function [mean_ftle, SEM_ftle, per_pos, firstframe, lastframe] = ftle_analytics(directory, imname, ftleLRC, r_scale, t_scale, firstframe, lastframe)
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%% ****PARAMETERS**** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These parameters are often constant between experiments, but should be
% confirmed for all data sets.
param.edge_cut = 0.75; % Only include vectors with r/R greater than this in edge metrics (set to 0 to analyze all)
param.center_cut = 0.5; % Only include vectors with r/R less than this in center metrics (set to 1 to analzye all)

param.ftle_bins = -5:.01:1; % Bins (in hr^-1) for FTLE histograms

param.R_bin_size = 0.05; % bin metrics in units of r/R
param.R_bin_max = 1; % Maximum value of r/R to plot
param.size_limit = 200; % Number of vectors needed to keep a bin for r/R plots

param.r_step = 50; % microns for non-scaled bins
param.numel_limit = 42; % Number of vectors needed to keep a bin for heatmap
param.heat_limits = [-1 .5]; % Range for ftle kymograph colormap, in hr^-1


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
    load([directory imname 'L_edgedat.mat'])
    DotL = edgedat;
    clear edgedat

    load([directory imname 'R_edgedat.mat'])
    DotR = edgedat;
    clear edgedat

    % All outputs in um
    [X,Y,R,xL,xR,yL,yR,xC,yC] = do_edgedat_to_circle(DotL, DotR, param.pos_list, directory, imname,r_scale);

    clear DotL DotR
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
[ftles,r_scaled,beginning, ending] = do_scale_shift_ftle(ftleLRC,dot_piv,xL,yL,xR,yR,xC,yC,X,Y,R,r_scale,t_scale);

% Update frame range to not exceed rational ranges
firstframe = max(firstframe,beginning);
lastframe = min(lastframe,ending);
ftles = ftles(:,:,firstframe:lastframe);
r_scaled = r_scaled(:,:,firstframe:lastframe);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot Kymograph %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load ftle_cmap.mat
title_text = 'FTLE Values';
bar_text = 'FTLE (hr^-^1)';
savename = [directory 'FTLE Analytics\' imname 'heat_map'];
% plot_kymo(ftles,r_scaled(:,:,1)*R(firstframe),t_scale,param.r_step,param.numel_limit,1,lastframe-firstframe+1,param.heat_limits,title_text,savename,bar_text,imname,ftle_cmap)
plot_kymo(ftles,r_scaled(:,:,1)*R(firstframe),t_scale,param.r_step,param.numel_limit,1,lastframe-firstframe+1,param.heat_limits,title_text,savename,bar_text,imname)

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Calc Outputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mean_ftle(1) = mean(ftles(~isnan(ftles)));
SEM_ftle(1) = std(ftles(~isnan(ftles)))/numel(find(~isnan(ftles)));
per_pos(1) = numel(find(ftles > 0))/numel(find(~isnan(ftles)))*100;

ftle_slice = ftles(r_scaled > param.edge_cut);
mean_ftle(2) = mean(ftle_slice(~isnan(ftle_slice)));
SEM_ftle(2) = std(ftle_slice(~isnan(ftle_slice)))/numel(find(~isnan(ftle_slice)));
per_pos(2) = numel(find(ftle_slice > 0))/numel(find(~isnan(ftle_slice)))*100;

ftle_slice = ftles(r_scaled < param.center_cut);
mean_ftle(3) = mean(ftle_slice(~isnan(ftle_slice)));
SEM_ftle(3) = std(ftle_slice(~isnan(ftle_slice)))/numel(find(~isnan(ftle_slice)));
per_pos(3) = numel(find(ftle_slice > 0))/numel(find(~isnan(ftle_slice)))*100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Histograms %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All vectors
savename = [directory 'FTLE Analytics\' imname 'hist_total'];
title_text = ['Total Positive: ' num2str(per_pos(1),'%0.0f') ' %'];
plot_ftle_pdf(ftles(~isnan(ftles)),param.ftle_bins,savename,title_text,imname)

% Edge vectors only
savename = [directory 'FTLE Analytics\' imname 'hist_edge'];
title_text = ['Edge Positive: ' num2str(per_pos(2),'%0.0f') ' %'];
ftle_slice = ftles(r_scaled > param.edge_cut);
plot_ftle_pdf(ftle_slice(~isnan(ftle_slice)),param.ftle_bins,savename,title_text,imname)

% Center vectors only
savename = [directory 'FTLE Analytics\' imname 'hist_center'];
title_text = ['Center Positive: ' num2str(per_pos(3),'%0.0f') ' %'];
ftle_slice = ftles(r_scaled < param.center_cut);
plot_ftle_pdf(ftle_slice(~isnan(ftle_slice)),param.ftle_bins,savename,title_text,imname)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Radial Values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bins = 0:param.R_bin_size:param.R_bin_max;

% Preallocate storage
bin_mid = NaN*ones(1,length(bins) - 1);
N = bin_mid;
r_mean = bin_mid;
r_std = bin_mid;
ftle_mean = bin_mid;
ftle_std = bin_mid;
ftle_pos = bin_mid;

% Perform binning
for j = 1:(length(bins) - 1)

    now1 = r_scaled > bins(j);
    now2 = r_scaled < bins(j+1);
    now3 = now1.*now2;
    now = find(now3);
    
    bin_mid(j) = (bins(j)+bins(j+1))/2;
    
    ftle_now = ftles(now);    
    ftle_now = ftle_now(~isnan(ftle_now));
    N(j) = numel(ftle_now);
    
    if N(j) > param.size_limit
    
        r_mean(j) = mean(r_scaled(now));
        r_std(j) = std(r_scaled(now));
        
        ftle_mean(j) = mean(ftle_now);
        ftle_std(j) = std(ftle_now);
        ftle_pos(j) = 100*numel(find(ftle_now > 0))/N(j);

    end

end

% Save FTLE binned data
save([directory 'FTLE Analytics\' imname '_ftle_vs_r.mat'],...
   'bin_mid','N','r_mean','r_std','ftle_mean','ftle_pos','ftle_std',...
   'firstframe','lastframe')

errorbar(bin_mid,ftle_mean,ftle_std./sqrt(N),'k-s','LineWidth',2,'MarkerFaceColor','k')
hold on
plot([0 1],[0 0],'--k')
hold off
xlabel('r/R','FontSize',16)
ylabel('<FTLE> (hr^-^1)','FontSize',16)
title([imname 10 'FTLE vs r/R'],'FontSize',16,'Interpreter','none')
set(gca,'FontSize',16)
xlim([0 max(bins)+param.R_bin_size])
% Save figure
savename = [directory 'FTLE Analytics\' imname 'Radial FTLE'];
saveas(gcf,[savename '.fig'],'fig')
saveas(gcf,[savename '.tif'],'tif')


plot(bin_mid,ftle_pos,'k-s','LineWidth',2,'MarkerFaceColor','k')
xlabel('r/R','FontSize',16)
ylabel('Positive FTLE Values (%)','FontSize',16)
title([imname 10 'Positive FTLE Values'],'FontSize',16,'Interpreter','none')
set(gca,'FontSize',16)
xlim([0 max(bins)+param.R_bin_size])
% Save figure
savename = [directory 'FTLE Analytics\' imname 'Radial Positives'];
saveas(gcf,[savename '.fig'],'fig')
saveas(gcf,[savename '.tif'],'tif')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Versus Time %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_t = NaN*ones(size(ftles,3),1);
mean_ftle_t = N_t;
std_ftle_t = N_t;
per_pos_t = N_t;
N_t_edge = N_t;
mean_ftle_t_edge = N_t;
std_ftle_t_edge = N_t;
per_pos_t_edge = N_t;
N_t_bulk = N_t;
mean_ftle_t_bulk = N_t;
std_ftle_t_bulk = N_t;
per_pos_t_bulk = N_t;

for kk = 1:size(ftles,3)
    
    slice = ftles(:,:,kk);
    N_t(kk) = numel(find(~isnan(slice)));
    mean_ftle_t(kk) = mean(slice(~isnan(slice)));
    std_ftle_t(kk) = std(slice(~isnan(slice)));
    per_pos_t(kk) = 100*numel(find(slice > 0))/N_t(kk);
    
    slice2 = slice(r_scaled(:,:,kk) > param.edge_cut);
    N_t_edge(kk) = numel(find(~isnan(slice2)));
    mean_ftle_t_edge(kk) = mean(slice2(~isnan(slice2)));
    std_ftle_t_edge(kk) = std(slice2(~isnan(slice2)));
    per_pos_t_edge(kk) = 100*numel(find(slice2 > 0))/N_t_edge(kk);
    
    slice2 = slice(r_scaled(:,:,kk) < param.center_cut);
    N_t_bulk(kk) = numel(find(~isnan(slice2)));
    mean_ftle_t_bulk(kk) = mean(slice2(~isnan(slice2)));
    std_ftle_t_bulk(kk) = std(slice2(~isnan(slice2)));
    per_pos_t_bulk(kk) = 100*numel(find(slice2 > 0))/N_t_bulk(kk);
    
end

hour_vec = (firstframe:lastframe)*t_scale/60;

errorbar(hour_vec,mean_ftle_t,std_ftle_t./sqrt(N_t),'k-s','LineWidth',2,'MarkerFaceColor','k')
hold on
errorbar(hour_vec,mean_ftle_t_edge,std_ftle_t_edge./sqrt(N_t_edge),'b-o','LineWidth',2,'MarkerFaceColor','b')
errorbar(hour_vec,mean_ftle_t_bulk,std_ftle_t_bulk./sqrt(N_t_bulk),'r-^','LineWidth',2,'MarkerFaceColor','r')
plot([0 1],[0 0],'--k')
hold off
xlabel('Time (hours)','FontSize',16)
ylabel('<FTLE> (hr^-^1)','FontSize',16)
title([imname 10 'FTLE vs Time'],'FontSize',16,'Interpreter','none')
legend('Total','Edge','Bulk')
set(gca,'FontSize',16)
xlim([0 max(hour_vec)])
% Save figure
savename = [directory 'FTLE Analytics\' imname 'Time FTLE'];
saveas(gcf,[savename '.fig'],'fig')
saveas(gcf,[savename '.tif'],'tif')

plot(hour_vec,per_pos_t,'k-s','LineWidth',2,'MarkerFaceColor','k')
hold on
plot(hour_vec,per_pos_t_edge,'b-o','LineWidth',2,'MarkerFaceColor','b')
plot(hour_vec,per_pos_t_bulk,'r-^','LineWidth',2,'MarkerFaceColor','r')
hold off
xlabel('Time (hours)','FontSize',16)
ylabel('Positive FTLE Values (%)','FontSize',16)
title([imname 10 'Positive FTLE Values'],'FontSize',16,'Interpreter','none')
legend('Total','Edge','Bulk')
set(gca,'FontSize',16)
xlim([0 max(hour_vec)])
% Save figure
savename = [directory 'FTLE Analytics\' imname 'Time Positives'];
saveas(gcf,[savename '.fig'],'fig')
saveas(gcf,[savename '.tif'],'tif')

save([directory 'FTLE Analytics\' imname 'FTLE_time_Analytics.mat'],...
    'hour_vec','mean_ftle_t','per_pos_t','N_t','std_ftle_t'...
    ,'mean_ftle_t_edge','per_pos_t_edge','N_t_edge','std_ftle_t_edge'...
    ,'mean_ftle_t_bulk','per_pos_t_bulk','N_t_bulk','std_ftle_t_bulk')

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Calc Outputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save([directory 'FTLE Analytics\' imname 'FTLEAnalytics.mat'],'mean_ftle'...
    ,'SEM_ftle', 'per_pos','firstframe','lastframe','param')
