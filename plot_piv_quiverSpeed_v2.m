% Create plots with quiver overlays on speed - uses dot_piv structure
% Usage: plot_piv_quiverSpeed(dot_piv, [bot], [top], [space], [frames],[qscale],[units],[savedir])

function plot_piv_quiverSpeed_v2(dot_piv, bot, top, space, frames,qscale,unitSpace,unitTime,savedir)

%%%% Inputs
% dot_piv:      structure output of dot_matpiv code suites
% bot:          (Optional) Top value for colorbar. Default = 1.
% top:          (Optional) Bottom value for colorbar. Default = 0.
% space:        (Optional) Conversion from pixels to unit of interest.
%                   Default value = 1.
% frames:       (Optional) Time between frames. Defaults = 1.
% qscale:       (Optional) Manipulates length of vectors. Default = 12. 
%                       Zero results in no scaling
% unitSpace:    (Optional) String which specifies units of space (e.g. um)
%                   Default = 'pixels'
% unitTime:     (Optional) String which specifies units of time (e.g. min)
%                   Default = '\Deltat'
% savedir:      (Optional) Location to save images.  Default is a folder
%                   created in the PIV directory with name dot_piv.imname
%
%%%% Outputs
% imname_quiverSpeed_0000.tif: Image sequence saved in savedir
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rachel Lee 2016/10/07 Based on plot_piv_pcolor
%                       Main difference: arrows in center of squares using
%                       imagesc instead of pcolor % more sensible defaults
% Initially tested in 2014b -- figure commands might not be compatible with
% other versions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fontsize = 16;
%%%% DEFAULT PARAMETERS
spaced = 1;
framesd = 1;
qscaled = 12;
topd = 1;
botd = 0;
unitSpaced = 'pixels';
unitTimed = '\Deltat';

if ~exist('space','var') || isempty(space)
    space = spaced;
end
if ~exist('frames','var') || isempty(frames)
    frames = framesd;
end
if ~exist('qscale','var') || isempty(qscale)
    qscale = qscaled;
end
if ~exist('bot','var') || isempty(bot)
    bot = botd;
end
if ~exist('top','var') || isempty(top)
    top = topd;
end
if ~exist('unitSpace','var') || isempty(unitSpace)
    unitSpace = unitSpaced;
end
if ~exist('unitTime','var') || isempty(unitTime)
    unitTime = unitTimed;
end
if ~exist('savedir','var') || isempty(savedir)
    % Check for correct directory formatting
    if ~strcmp(dot_piv.directory(end),filesep)
        dot_piv.directory = [dot_piv.directory filesep];
    end
    % Set up a directory based on image name
    if ~isempty(dot_piv.imname)
        savedir = [dot_piv.directory 'quiverSpeed_' dot_piv.imname filesep];
    else
        savedir = [dot_piv.directory 'quiverSpeed' filesep];
    end
end

% If the save folder doesn't exist, create it
if ~exist(savedir,'file')
    mkdir(savedir)
end

%%%%%%%%

% Calculate speed and convert units
speed = sqrt(dot_piv.u_fi.^2 + dot_piv.v_fi.^2);
speed = speed*space/frames;

% Get current frames of velocity
dot_piv.u_fi = space/frames*dot_piv.u_fi;
dot_piv.v_fi = space/frames*dot_piv.v_fi;
dot_piv.x = dot_piv.x*space;
dot_piv.y = dot_piv.y*space;

% Axis limits
x_plot = [dot_piv.x(1) dot_piv.x(end)]; % imagesc allows you to specify the corners of the matrix
y_plot = [dot_piv.y(1) dot_piv.y(end)];
limits = [x_plot y_plot];

% Find time range
lastframe = dot_piv.lastframe;
firstframe = dot_piv.firstframe;
if ~firstframe
    firstframe = firstframe + 1;
    lastframe = lastframe + 1;
end

% Make the images larger
scrsz = get(0,'ScreenSize');
figure('Position',[50 100 3*scrsz(3)/4 3*scrsz(4)/4]);
pause(1) % Allows user to change figure size

% Create frames from each image
for p = firstframe:(lastframe-1)
    
    % imagesc plots the speed as colors
    h = imagesc(x_plot,y_plot,speed(:,:,p));
    set(h,'alphadata',~isnan(speed(:,:,p))); % don't show the NaN regions
    % Set up the color range
    caxis manual
    caxis([bot top]);
    h = colorbar;
    set(get(h,'Label'),'String',['Speed (' unitSpace '/' unitTime ')'],'FontSize',fontsize)
    
    % quiver plots the flow field arrows
    hold on
    j = quiver(dot_piv.x,dot_piv.y,dot_piv.u_fi(:,:,p),dot_piv.v_fi(:,:,p),0,'k');
    if qscale ~= 0
        hU = get(j,'UData') ;
        hV = get(j,'VData') ;
        set(j,'UData',qscale*hU,'VData',qscale*hV)
    end
    hold off
    
    %Organize the plot
    axis(limits)
    xlabel(unitSpace,'FontSize',fontsize)
    set(gca,'DataAspectRatio',[1 1 1],'FontSize',fontsize)
    
    % Save the figure
    saveas(gcf,[savedir dot_piv.imname '_quiverSpeed_' num2str(p,'%04u') '.tif'],'tif')
    
%     imw = getframe;
%     imwrite(imw.cdata(:,:,:),[savedir '\PIV_' num2str(p) '.tif'],'tif', 'Compression','none');


end
