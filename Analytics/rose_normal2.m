% Make a rose plot normalized by total count
% Patch does not work correctly in 2014b, color the lines instead
% Usage: rose_normal(angles,angle_bins,rose_limit,savename,line_color,title_text)

function rose_normal2(angles,angle_bins,rose_limit,savename,line_color,title_text,y_text)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS
% angles      = angles to create histogram
% bins        = number of bins for histogram
% rose_limit  = Optional.  Set scale of rose plot.
% savename    = Optional. Save figure to file.
% line_color = Optional. Color for rose plot lines
% title_text  = Optional. Text for title (placed in xlabel).
% y_text      = Optional. Label for graph (placed in ylabel). No tex
%               interpretaton, to allow for Colibri file names
%
% Created by Rachel Lee, 2014/02/06
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('line_color','var') || isempty(line_color)
    line_color = [0 0.4470 0.7410];
end

% bins = (0:angle_bins-1)*2*pi/angle_bins + pi/angle_bins;
[t,r] = rose(angles(:),angle_bins);

% Plot the initial rose plot
if exist('rose_limit','var') || ~isempty(rose_limit)
    h_fake = polar([0 pi],rose_limit*[1 1]);
    hold on
    h = polar(t,2*r/sum(r),'k'); % Normalize by total count (i.e. output is percent)
    hold off
    set(h_fake,'Visible','off');
else
    h = polar(t,2*r/sum(r)); % Normalize by total count (i.e. output is percent)
end

% Format rose plot
set(findall(gcf, 'String', '30', '-or','String','60',...
    '-or','String','120', '-or','String','150',...
    '-or','String','210', '-or','String','240',...
    '-or','String','300', '-or','String','330'),'String', '  ');
set(findall(gcf,'String','  0.05','-or','String','  0.15','-or','String','  0.25'),'String',' ');
set(findall(gcf, 'type','text'),'FontSize',16);

% Make it prettier
set(h,'LineWidth',2,'Color',line_color)

if exist('title_text','var') && ~isempty(title_text)
    xlabel(title_text,'FontSize',16)
end

if exist('y_text','var') && ~isempty(y_text)
    ylabel(y_text,'FontSize',16,'Interpreter','none')
end

% Save if there is a file name given
if exist('savename','var') && ~isempty(savename)
    saveas(gcf,[savename '.fig'],'fig')
    saveas(gcf,[savename '.tif'],'tif')
end