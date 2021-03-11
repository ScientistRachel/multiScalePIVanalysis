% Creates images comparing edge detection to PIV filtering.
%
% Usage: filter_check(directory, imname, frame_skip)

function filter_check(directory, imname, frame_skip)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS
% directory = folder containing PIV and edgedat.mat
% imname = image name
% frame_skip = how often to plot overlays
% OUTPUTS
% No variables are sent to the workspace, but images are saved in the
% directory in a subfolder 'FilterCheck'.
%
% Created by Rachel Lee, 2014/02/10
% CHANGE LOG
% 2014/05/27 RML Indexing error if edge or PIV did not start on frame 1
% 2015/12/07 RML flip the axis in the right direction for user's ease of
% use in comparing to edge detection images.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If necessary, create a folder for saving
if ~exist([directory 'FilterCheck\' imname '\'],'file')
    mkdir([directory 'FilterCheck\'],imname)
end

% Load in the revelant data
load([directory imname '_edgedat.mat'])
load([directory imname 'dot_piv_filtered.mat'])

% Make sure there is no frame shift between PIV and edge detection
piv1 = dot_piv.firstframe;
piv2 = size(dot_piv.u,3) + piv1 - 1;

% Make sure that we only analzye frames that have both an edge and PIV
if isfield(edgedat,'firstframe')
    edge1 = edgedat.firstframe;
else
    edge1 = 1;
end
edge2 = size(edgedat.points,3) + edge1 - 1;

max_frame = min(piv2,edge2);
min_frame = max(piv1,edge1);

% Make the images to save
for kk = min_frame:frame_skip:max_frame
    
    % Get current frame
    piv_slice = dot_piv.u_fi(:,:,kk-piv1+1);
    edge_slice = edgedat.points(:,:,kk-edge1+1);
    
    % Black and white PIV
    piv_slice(~isnan(piv_slice)) = 1;
    piv_slice(isnan(piv_slice)) = 0;
    
    % Plot it
    pcolor(dot_piv.x,dot_piv.y,piv_slice), shading flat
    colormap gray
    hold on
    plot(edge_slice(:,2),edge_slice(:,1),'g','LineWidth',2)
    hold off
    title(imname,'FontSize',16,'Interpreter', 'none')
    set(gca,'DataAspectRatio',[1 1 1],'ydir','reverse')
    axis off
    
    % Save it
    savename = [directory 'FilterCheck\' imname '\Frame' num2str(kk)];
    saveas(gcf,[savename '.tif'],'tif')    
    
end