function [b,da]=stretchg(u,v,x,y,savefiles,noisy)
% Usage: [b,da]=stretchg(u,v,x,y,[savefiles],[noisy])
% Working from the velocity field specified by "u", "v", "x", and "y",
% stretchg calculates its stretching field, returning the results in "b".
% An estimate of the angular error is also returned in "da". If "savefiles"
% is a string, results at each time step are saved as .mat files in a
% subdirectory with that name. Otherwise, if savefiles~=0, results are
% saved in the current directory. The array "u" should specify the
% x-direction velocity on a three-dimensional (x,y,time) grid. Likewise
% "v" should specify the y-direction velocity on a grid, and "x" and "y" 
% should give the locations of the grid points. If noisy==0, no status
% udpates are displayed. See also stretch.m. 

% Written 14 February 2012 by Doug Kelley, based largely on stretch.m. 

savefilesdefault=0;
noisydefault=1;

if nargin<4
    error(['Usage: [b,da] = ' mfilename '(u,v,x,y,[savefiles],[noisy])'])
end
if ~exist('savefiles','var') || isempty(savefiles)
    savefiles=savefilesdefault;
elseif ischar(savefiles)
    savefiles=['./' savefiles]; % to avoid seeing stuff elsewhere on path with same name
    if exist(savefiles,'file')~=7
        mkdir(savefiles)
    end
elseif savefiles % anything but 0
    savefiles=pwd;
end
if ~exist('noisy','var') || isempty(noisy)
    noisy=noisydefault;
end

if isvector(x)
    [x,y]=meshgrid(x,y);
end
sx=size(x);
sy=size(y);
su=size(u);
sv=size(v);
if any(sx(1:2)~=su(1:2)) || any(sy(1:2)~=su(1:2)) || ...
    any(sx(1:2)~=sv(1:2)) || any(sy(1:2)~=sv(1:2))
    error('Sorry, x and y must match the size of u and v.')
end

T=size(u,3);
xp=x; % coordinates of virtual tracers (start on grid)
yp=y;
U2=u(:,:,1); % initial velocity field
V2=v(:,:,1);
dx=x(1,2)-x(1,1); % x grid size
dy=y(2,1)-y(1,1); % y grid size
b=NaN(size(x)); % stretching
da=NaN(size(x)); % error estimate

% -=- Integrate fluid elements -=-
for ii=1:T-1
    if noisy && ~mod(ii-1,100)
        disp(['Now calculating frame ' num2str(ii) ' of ' num2str(T) '.']);
    end
    U1=U2; % velocity on grid
    V1=V2;
    u1=interp2(x,y,U1,xp,yp,'*cubic'); % velocity at virtual tracers
    v1=interp2(x,y,V1,xp,yp,'*cubic');
    u1(isnan(u1))=0;
    v1(isnan(v1))=0;
    xp1=xp+u1; % new tracer locations, first guess
    yp1=yp+v1;

    U2=u(:,:,ii+1); % later velocity on grid
    V2=v(:,:,ii+1);
    u2=interp2(x,y,U2,xp1,yp1,'*cubic'); % later velocity at virtual tracers
    v2=interp2(x,y,V2,xp1,yp1,'*cubic');
    u2(isnan(u2))=0;
    v2(isnan(v2))=0;
    xp=xp+(u1+u2)/2; % new tracer locations, improved guess
    yp=yp+(v1+v2)/2;

    drx=xp-x; % offset from original location
    dry=yp-y;

    % -=- Finite differences for gradients -=-
    sd = size(drx);
    drx_dx=(drx(3:sd(1),2:sd(2)-1)-drx(1:sd(1)-2,2:sd(2)-1))/2/dx;
    drx_dy=(drx(2:sd(1)-1,3:sd(2))-drx(2:sd(1)-1,1:sd(2)-2))/2/dy;
    dry_dx=(dry(3:sd(1),2:sd(2)-1)-dry(1:sd(1)-2,2:sd(2)-1))/2/dx;
    dry_dy=(dry(2:sd(1)-1,3:sd(2))-dry(2:sd(1)-1,1:sd(2)-2))/2/dy;

    % -=- Compute eigenvalues -=-
    P = drx_dx.^2 + drx_dy.^2 + dry_dx.^2 + dry_dy.^2;
    Q = (drx_dx.*dry_dy).^2 + (drx_dy.*dry_dx).^2 - ...
        2.0 * drx_dx.*drx_dy.*dry_dx.*dry_dy;
    L1 = P/2.0 + sqrt((P/2.0).^2 - Q);
    L2 = P/2.0 - sqrt((P/2.0).^2 - Q);
    b(2:end-1,2:end-1)=sqrt(max(abs(L1),abs(L2))); % stretching
    b(1,:)=b(2,:); % fill in the edges
    b(end,:)=b(end-1,:);
    b(:,1)=b(:,2);
    b(:,end)=b(:,end-1); 
    b(1,1)=(b(1,2)+b(2,1))/2; % and corners
    b(1,end)=(b(2,end)+b(1,end-1))/2;
    b(end,1)=(b(end-1,1)+b(end,2))/2;
    b(end,end)=(b(end,end-1)+b(end-1,end))/2;
    da(2:end-1,2:end-1)=asin(sqrt(min(abs(L1),abs(L2))) ./ ...
        b(2:end-1,2:end-1)); % error estimate
    da(1,:)=da(2,:); % fill in the edges
    da(end,:)=da(end-1,:);
    da(:,1)=da(:,2);
    da(:,end)=da(:,end-1); 
    da(1,1)=(da(1,2)+da(2,1))/2; % and corners
    da(1,end)=(da(2,end)+da(1,end-1))/2;
    da(end,1)=(da(end-1,1)+da(end,2))/2;
    da(end,end)=(da(end,end-1)+da(end-1,end))/2;

    if savefiles
        save(fullfile(savefiles,[num2str(ii,'%05.0f') '.mat']), ...
            'b','x','y','da');
    end
    
end % for ii=1:T
if noisy
    disp('Done.')
end

