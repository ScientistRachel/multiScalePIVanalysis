% Plot a histogram normalized such that area under the curve is 1.
%
% Usage: plot_pdf(values,bins,savename,title_text,y_text)

function plot_ftle_pdf(values,bins,savename,title_text,y_text)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS
% values      = FTLE values to create histogram (hr^-1)
% bins        = bins for histogram (hr^-1)
% savename    = Optional. Save figure to file.
% title_text  = Optional. Text for title.
% y_text      = Optional. Label for graph (placed in ylabel). No tex
%               interpretaton, to allow for Colibri file names.
%
% Created by Rachel Lee, 2014/02/20
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Take histogram and normalize
n = hist(values,bins);
n = n/trapz(bins,n);

% Plot
plot(bins,n,'LineWidth',2)
xlabel('FTLE Value (hr^-^1)','FontSize',16)
set(gca,'FontSize',16)

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