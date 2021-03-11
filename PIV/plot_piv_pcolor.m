% Create pcolor plots with quiver overlays of the information in dot_piv
% Usage: plot_piv_pcolor(dot_piv, [bot], [top], [space], [frames],[qscale])

function plot_piv_pcolor(dot_piv, bot, top, space, frames,qscale)

%%%% Inputs
% dot_piv:  structure output of dot_matpiv (can also be output of trimming)
% bot:      (Optional) Top value for colorbar. Defaults to 1.
% top:      (Optional) Bottom value for colorbar. Defaults to 0.
% space:    (Optional) Conversion from pixels to unit of interest. Defaults
%               to 0.65
% frames:   (Optional) Time between frames. Defaults to 5.
% qscale:   (Optional) Manipulates length of vectors. Defaults to 10.
%
%%%% Outputs
%   - PIV_0.tif - Image sequence saved in a subfolder labeled dot_piv.imname

% Rachel Lee 2013/04/11 Based on previously used code, but cleaned up
% 2013-12-19 RML added ability to make a folder even if imname = [];

%%%% DEFAULT PARAMETERS
spaced = 0.65;
framesd = 5;
qscaled = 12;
topd = 1;
botd = 0;

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

%%%%%%%%

Z = sqrt(dot_piv.u_fi.^2 + dot_piv.v_fi.^2);
Z = Z*space/frames;

limits = [min(dot_piv.x(:)) max(dot_piv.x(:)) min(dot_piv.y(:)) max(dot_piv.y(:))];

if ~isempty(dot_piv.imname)
    mkdir(dot_piv.directory, dot_piv.imname);
    long_dir = [dot_piv.directory '\' dot_piv.imname];
else
    mkdir(dot_piv.directory,'PlotPColor')
    long_dir = [dot_piv.directory '\PlotPColor'];
end

% Start up the waitbar with an initial guess of 2 seconds per fame
h = waitbar(0,'Writing Frame: 1');
t_ave = 2;

lastframe = dot_piv.lastframe;
firstframe = dot_piv.firstframe;

if ~firstframe
    firstframe = firstframe + 1;
    lastframe = lastframe + 1;
end


% Make the images larger
scrsz = get(0,'ScreenSize');
figure('Position',[50 100 3*scrsz(3)/4 3*scrsz(4)/4]);
pause(1) % Allows user to change figure size (and move waitbar if desired)

if ispc

% Create frames from each image
for p = firstframe:(lastframe-1)
 
    % begin recording the average time of computation for this frame
    tic
    
    % calculate the estimated time of finishing the sequence in minues
    % and seconds, and display the estimate in the waitbar
    mins = floor(t_ave*(lastframe-p+1)/60);
    secs = round(floor((t_ave*(lastframe-p+1))) - mins*60);
    waitbar((p-1)/(lastframe),h,sprintf('%s%u%s%u%s%u',...
        'Writing frame: ',...
        p,...
        '      Est time remaining: ',...
        mins,...
        'mins ',...
        secs,...
        'secs'));
    
    % Plot the interpolated velocity with quiver and contourf
    u_slice = space/frames*dot_piv.u_fi(:,:,p);
    v_slice = space/frames*dot_piv.v_fi(:,:,p);
    Z_slice = Z(:,:,p);
    

    pcolor(dot_piv.x,dot_piv.y,Z_slice), shading flat
    caxis manual
    caxis([bot top]);
    colorbar
    hold on
    j = quiver(dot_piv.x,dot_piv.y,u_slice,v_slice,0,'k','LineWidth',2);
    hU = get(j,'UData') ;
    hV = get(j,'VData') ;
    set(j,'UData',qscale*hU,'VData',qscale*hV)
    axis(limits)
    set(gca,'YDir','reverse','DataAspectRatio',[1 1 1])
    hold off
    
    imw = getframe;
    imwrite(imw.cdata(:,:,:),[long_dir '\PIV_' num2str(p) '.tif'],'tif', 'Compression','none');

    % record running average computation time per frame
    time = toc;
    t_ave = (t_ave+time)/2;
end

else

% Create frames from each image
for p = firstframe:(lastframe-1)
 
    % begin recording the average time of computation for this frame
    tic
    
    % calculate the estimated time of finishing the sequence in minues
    % and seconds, and display the estimate in the waitbar
    mins = floor(t_ave*(lastframe-p+1)/60);
    secs = round(floor((t_ave*(lastframe-p+1))) - mins*60);
    waitbar((p-1)/(lastframe),h,sprintf('%s%u%s%u%s%u',...
        'Writing frame: ',...
        p,...
        '      Est time remaining: ',...
        mins,...
        'mins ',...
        secs,...
        'secs'));
    
    % Plot the interpolated velocity with quiver and contourf
    u_slice = space/frames*dot_piv.u_fi(:,:,p);
    v_slice = space/frames*dot_piv.v_fi(:,:,p);
    Z_slice = Z(:,:,p);
    
    pcolor(dot_piv.x,dot_piv.y,Z_slice), shading flat
    caxis manual
    caxis([bot top]);
    colorbar
    hold on
    j = quiver(dot_piv.x,dot_piv.y,u_slice,v_slice,0,'k','LineWidth',2);
    hU = get(j,'UData') ;
    hV = get(j,'VData') ;
    set(j,'UData',qscale*hU,'VData',qscale*hV)
    axis(limits)
    set(gca,'YDir','reverse','DataAspectRatio',[1 1 1])
    hold off
    
    imw = getframe;
    imwrite(imw.cdata(:,:,:),[dot_piv.directory '/' dot_piv.imname '/PIV_' num2str(p) '.tif'],'tif', 'Compression','none');

    % record running average computation time per frame
    time = toc;
    t_ave = (t_ave+time)/2;
end

end

% delete the waitbar
delete(h);

end
