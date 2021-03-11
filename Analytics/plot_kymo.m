% Plot a colormap coded kymograph of a variable over both a length and time
% scale.  Note that time always starts from zero, regardless of first frame
% analyzed.
%
% Usage: plot_kymo(heat_var,len_var,t_scale,r_step,numel_limit,firstframe,lastframe,heat_limits,title_text,savename,bar_text,y_text)

function plot_kymo(heat_var,len_var,t_scale,r_step,numel_limit,firstframe,lastframe,heat_limits,title_text,savename,bar_text,y_text)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS
% heat_var      = MxNxP matrix to be plotted with color representing magnitude
% len_var       = MxN matrix giving the location of each MxN slice of heat_var
% t_scale       = time scale between MxN slices of heat_var (in minutes)
% r_step        = bin size (in microns)
% numel_limit   = number of data points necessary to keep a bin
% firstframe    = firstframe to include in kymograph
% lastframe     = lastframe to include in kymograph
% heat_limits   = limits for the colormap/colorbar
% title_text    = Optional. text to include in the title
% savename      = Optional. string to use for saving the figure.
% bar_text      = Optional. string to label colorbar
% y_text        = Optional. string to indentify graph. No tex interpeter to
%                 allow for Colibri file names
%
% OUTPUTS
% No outputs are sent to the workspace, but a figure is generated and
% optionally saved.
%
% Created by Rachel Lee 2014/02/12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('title_text','var')
    title_text = [];
end

if ~exist('bar_text','var')
    bar_text = [];
end

if exist('y_text','var') && ~isempty(y_text)
    y_text = [y_text 10 10 'Radial Location (mm)'];
else
    y_text = 'Radial Location (mm)';
end

% Set up bins
bins = 0:r_step:(max(len_var(:)));
% Find time limits
total_time = lastframe - firstframe + 1;
% Preallocate storage
to_plot = zeros(length(bins)-1,total_time);

% Bin heat_var using len_var
for kk = 1:(length(bins)-1)
    
    good1 = len_var > bins(kk);
    good2 = len_var < bins(kk+1);
    good = good1 & good2;
    
    for jj = 1:total_time
        
        slice = heat_var(:,:,jj+firstframe-1);
        slice = slice(good);
        
        if numel(slice(~isnan(slice))) >= numel_limit
            to_plot(kk,jj) = mean(slice(~isnan(slice)));
        else
            to_plot(kk,jj) = NaN;
        end
        
    end
    
end

% Set up axis vectors
time_vec = (0:(total_time-1))*t_scale/60; %time in hours
len_vec = (diff(bins)/2 + bins(1:end-1))/1000;

% Plot the heatmap
pcolor(time_vec,len_vec,to_plot), shading flat
xlabel('Time (hours)','FontSize',16)
ylabel(y_text,'FontSize',16,'Interpreter','none')
title(title_text,'FontSize',16)
set(gca,'FontSize',16)
caxis manual
caxis(heat_limits)
set(gcf,'renderer','zbuffer');
h = colorbar;
set(h,'FontSize',16)
set(get(h,'YLabel'),'String',bar_text,'FontSize',16)
ylim([0 (max(len_vec) + r_step/1000)])

% Optionally, save the figure
if exist('savename','var') && ~isempty(savename)
    saveas(gcf,[savename '.fig'],'fig')
    saveas(gcf,[savename '.tif'],'tif')
end
