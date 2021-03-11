% Function to provide some automated analysis of FTLE data.
% Please do not use this function without thinking!
% - Watch the original movie
% - Check the hard coded parameters
% - Watch that the results seem consistent!
%
% Usage: [mean_ftle, SEM_ftle, per_pos, firstframe, lastframe] = ftle_analytics(directory, imname, r_scale, t_scale, firstframe, lastframe)

function plot_ftle_time_rPairedTif(directory, imnames, ftleLRC, r_scale, t_scale, firstframe, lastframe,ylimits)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% INPUTS %%%
% directory     = folder of images
% imname        = file names (missing L, R, or C at the end)
% r_scale       = magnification of images in microns per pixel
% t_scale       = time scale of movie, in minutes per frame
% firstframe    = Optional. First frame to analyze. Defaults to 1.
% lastframe     = Optional. Last frame to analyze. Defaults to full image series.
% ylimits       = Optional. Set the y axis scale. 2 x 2 - first row sets
%                   FTLE graph scale, second row sets pos percent
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
% 2015/06/29 RML Update for multiC, update colormap, clean up figures
% 2015/12/14 edgedat file was loaded incorrectly!!! Fixed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%% ****PARAMETERS**** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These parameters are often constant between experiments, but should be
% confirmed for all data sets.
param.R_bin_size = 0.05; % bin metrics in units of r/R
param.R_bin_max = 1; % Maximum value of r/R to plot
param.size_limit = 5000; % Number of vectors needed to keep a bin for r/R plots
param.time_bin = 20; %Number of frames per time bin

param.pos_list = 1; % 1 = excel file exists with locations, 0 = find from zvi files (slower)

if ~exist('firstframe','var') || isempty(firstframe)
    firstframe = 1;
end

if ~exist('lastframe','var') || isempty(lastframe)
    lastframe = inf;
end

%%%%%%%%%%%%%%%%%%%%%%% Make Folder for Saving %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist([directory 'FTLE Time and R\'],'file')
    mkdir(directory,'FTLE Time and R')
    
end
savedir = [directory 'FTLE Time and R\'];

if exist('ylimits','var') && ~isempty(ylimits)
    if ~exist([directory 'FTLE Time and R\Same Scale\'],'file')
        mkdir([directory 'FTLE Time and R\'],'Same Scale')
    end
    savedir = [directory 'FTLE Time and R\Same Scale\'];
end

imname = [imnames{1} '_' imnames{2}];
%%%%%%%%%%%%%%%%%%%%%%%%%%% Circle Fitting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if exist([directory 'Analytics Data\' imname '_dot_RvsT.mat'],'file')
    
    load([directory 'Analytics Data\' imname '_dot_RvsT.mat'],'X','Y','R','x_pos','y_pos')
    if exist('x_pos','var')
        type = 'new_version';
    else
        type = 'old_version';
    end
    
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
    
    type = 'new_version';
    
end

% disp(type)

%%%%%%%%%%%%%%%%%%%%%%% Convert FTLE to Radial %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Shift PIV locations and set to appropriate scale
% Function assumes that all three PIV files (L, R, C), have the same
% positions, so only loads one file to get x and y values
% load([directory imname 'Cdot_piv.mat'])
% Other options for files:
load([directory imnames{1} 'dot_piv_filtered.mat'],'dot_piv')
% load([directory imname 'Rdot_piv_filtered.mat'])

% ftles = hr^-1; r_scaled = unitless; beginning,ending = frames
if strcmp(type,'new_version')
    [ftles,r_scaled,beginning, ending] = do_scale_shift_ftleMultiC(ftleLRC,dot_piv,x_pos,y_pos,X,Y,R,r_scale,t_scale);
elseif strcmp(type,'old_version')
    % 2016/10/04
    %old version expect LRC, but new script reads in CLR if standard naming
    % manually reorganize:
    ftleCLR = ftleLRC;
    ftleLRC{1} = ftleCLR{2};
    ftleLRC{2} = ftleCLR{3};
    ftleLRC{3} = ftleCLR{1};
    clear ftleCLR
    [ftles,r_scaled,beginning, ending] = do_scale_shift_ftle(ftleLRC,dot_piv,xL,yL,xR,yR,xC,yC,X,Y,R,r_scale,t_scale);
else
    error('Something is wrong with your dot_RvsT.mat file')
end

% Update frame range to not exceed rational ranges
firstframe = max(firstframe,beginning);
lastframe = min(lastframe,ending);
ftles = ftles(:,:,firstframe:lastframe);
r_scaled = r_scaled(:,:,firstframe:lastframe);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Radial Values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bins = 0:param.R_bin_size:param.R_bin_max;
time_pieces = floor(size(ftles,3)/param.time_bin);

% Preallocate storage
bin_mid = NaN*ones(time_pieces,length(bins) - 1);
N = bin_mid;
r_mean = bin_mid;
r_std = bin_mid;
ftle_mean = bin_mid;
ftle_std = bin_mid;
ftle_pos = bin_mid;


% Perform binning
for kk = 1:time_pieces
    
    time_chunk = (param.time_bin*(kk-1)+1):(param.time_bin*kk);
    ftles_piece = ftles(:,:,time_chunk);
    r_piece = r_scaled(:,:,time_chunk);
    
    for j = 1:(length(bins) - 1)

        now1 = r_piece > bins(j);
        now2 = r_piece < bins(j+1);
        now3 = now1.*now2;
        now = find(now3);

        bin_mid(kk,j) = (bins(j)+bins(j+1))/2;

        ftle_now = ftles_piece(now);    
        ftle_now = ftle_now(~isinf(ftle_now));
        N(kk,j) = numel(ftle_now);

        if N(kk,j) > param.size_limit

            r_mean(kk,j) = mean(r_piece(now));
            r_std(kk,j) = std(r_piece(now));

            ftle_mean(kk,j) = mean(ftle_now);
            ftle_std(kk,j) = std(ftle_now);
            ftle_pos(kk,j) = 100*numel(find(ftle_now > 0))/N(kk,j);

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
    color_label{kk} = num2str(time_vec(kk));
end
x_tick_places = 0:(1/time_pieces):((time_pieces-1)/time_pieces);
x_tick_places = x_tick_places + mean(diff(x_tick_places))/2;

figure(1)
hold on
for kk = 1:time_pieces
    plot(bin_mid(kk,:),ftle_mean(kk,:),'-','LineWidth',2,'Color',cmap(kk,:))
end
hold off
xlabel('r/R','FontSize',16)
ylabel('<FTLE> (hr^-^1)','FontSize',16)
title([imname 10 'FTLE vs r/R'],'FontSize',16,'Interpreter','none')
set(gca,'FontSize',16)
xlim([0 max(bins)+param.R_bin_size])
if exist('ylimits','var') && ~isempty(ylimits)
    ylim(ylimits(1,:))
% else
%     a = get(gca,'YLim');
%     disp(a)
end
colormap(cmap)
h = colorbar;
set(h,'Ticks',x_tick_places,'TickLabels',color_label,'TickLength',0)
set(get(h,'Label'),'String','Hour','FontSize',16)
% Save figure
savename = [savedir imname 'FTLE Value'];
saveas(gcf,[savename '.fig'],'fig')
saveas(gcf,[savename '.tif'],'tif')

figure(2)
hold on
for kk = 1:time_pieces
    plot(bin_mid(kk,:),ftle_pos(kk,:),'-','LineWidth',2,'Color',cmap(kk,:))
end
hold off
xlabel('r/R','FontSize',16)
ylabel('Positive FTLE Values (%)','FontSize',16)
title([imname 10 'Positive FTLE Values'],'FontSize',16,'Interpreter','none')
set(gca,'FontSize',16)
xlim([0 max(bins)+param.R_bin_size])
if exist('ylimits','var') && ~isempty(ylimits)
    ylim(ylimits(2,:))
% else
%     a = get(gca,'YLim');
%     disp(a)
end
colormap(cmap)
h = colorbar;
set(h,'Ticks',x_tick_places,'TickLabels',color_label,'TickLength',0)
set(get(h,'Label'),'String','Hour','FontSize',16)
% % Save figure
savename = [savedir imname 'Positive FTLEs'];
saveas(gcf,[savename '.fig'],'fig')
saveas(gcf,[savename '.tif'],'tif')

