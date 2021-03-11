% Plot a histogram normalized such that area under the curve is 1.
%
% Usage: plot_pdf(values,bins,savename,title_text,y_text)

function plot_ftle_pdf_time(values,bins,time_bin,t_scale,savename,title_text,y_text)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS
% values      = FTLE values to create histogram (hr^-1) (M x N x P)
% bins        = bins for histogram (hr^-1)
% savename    = Optional. Save figure to file.
% title_text  = Optional. Text for title.
% y_text      = Optional. Label for graph (placed in ylabel). No tex
%               interpretaton, to allow for Colibri file names.
%
% Created by Rachel Lee, 2014/02/20
% Adapated to add a time component (different colored lines) on 2014/04/01
% Added ability to take different sections of dot 2014/04/02
% 2015/06/15 changed to parula colormap,removed r_scaling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

time_pieces = floor(size(values,3)/time_bin);
time_step = time_bin*t_scale/60;
time_vec = 1:time_step:time_step*time_pieces;
color_label = cell(1,length(time_vec));
for kk = 1:length(time_vec)
    color_label{kk} = num2str(time_vec(kk));
end
x_tick_places = 0:(1/time_pieces):((time_pieces-1)/time_pieces);
x_tick_places = x_tick_places + mean(diff(x_tick_places))/2;

% Take histogram and normalize
n = NaN*ones(length(bins),time_pieces);
for kk = 1:time_pieces
    
    frames = (time_bin*(kk-1)+1):(time_bin*kk);
    
    slice = values(:,:,frames);
    slice = slice(~isinf(slice));
    
    n(:,kk) = hist(slice,bins);
    n(:,kk) = n(:,kk)/trapz(bins,n(:,kk));
    
end

cmap = parula(size(n,2));
hold on
for kk = 1:size(n,2)
    plot(bins,n(:,kk),'LineWidth',2,'Color',cmap(kk,:))
end
xlabel('FTLE Value (hr^-^1)','FontSize',16)
set(gca,'FontSize',16)
xlim([min(bins) max(bins)])
colormap(cmap)
% caxis manual
% caxis([1 time_pieces])
h = colorbar;
set(h,'Ticks',x_tick_places,'TickLabels',color_label,'TickLength',0)
set(get(h,'Label'),'String','Hour','FontSize',16)

% Optional label with image name
if exist('y_text','var') && ~isempty(y_text)
    ylabel([y_text 10 10 'Probability Density'],'FontSize',16,'Interpreter','none')
else
    ylabel('Probability Density','FontSize',16)
end

% Optionally title
if exist('title_text','var') && ~isempty(title_text)
    title(title_text,'FontSize',16)
end

% Save if there is a file name given
if exist('savename','var') && ~isempty(savename)
    saveas(gcf,[savename '.fig'],'fig')
    saveas(gcf,[savename '.tif'],'tif')
end