% Function to shift multiple FTLE files into the frame of the microscope
% stage.  Analysis is limited by the extend of edge detection data.
%
% Usage: [ftles,r_scaled,beginning, ending] = do_scale_shift_ftle(ftleLRC,piv,xL,yL,xR,yR,xC,yC,X,Y,R,r_scale,t_scale)

function [ftles,r_scaled,beginning, ending] = do_scale_shift_ftleMultiC(ftleLRC,piv,x_pos,y_pos,X,Y,R,r_scale,t_scale)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS
% ftleLRC           = cell with {1} = FTLE for left, {2} = R, {3} = C
% piv               = any PIV file with the correct x and y positons
% xL,yL,xR,yR,xC,yC = Microscope stage postions (um)
% X,Y               = Dot center locations (um)
% R                 = radius of dot(um)
% r_scale           = microns per pixel
% t_scale           = minutes per frame
%
% OUTPUTS
% ftles             = combined ftle values (hr^-1)
% r_scaled          = radial location of PIV vectors, scaled by dot radius
% beginning         = first frame valid for analysis
% ending            = last frame valid for analysis
%
% Created by Rachel Lee, 2014/02/19, based on do_scale_shift_piv.m
% 2015/03/18 - Reinstanted lines that turn inf into NaN
% 2015/12/14 RML - fixed a problem where R isn't interpolated properly if
% more than one frame in a row are bad.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Remove edge problems from FTLE data
ftles = [];
for kk = 1:length(ftleLRC)
    slice = ftleLRC{kk};
    slice = slice*60/t_scale; %60/t_scale puts it in hr^-1
    %Get rid of edge effects
    slice(1:2,:,:) = [];
    slice(end-1:end,:,:) = [];
    slice(:,1:2,:) = [];
    slice(:,end-1:end,:) = [];
    % Inf -> NaN
    slice(isinf(slice)) = NaN;
    % Save the tweaked version
    ftles = [ftles slice];
end

% Don't use more frames than available in FTLE data
max_frame = min([size(ftles,3),length(R)]);
R = R(1:max_frame);

%%%%%%% PIV coordinates
x = r_scale*piv.x;
y = r_scale*piv.y;
clear dot_piv

x([1:2,end-1:end],:) = [];
x(:,[1:2,end-1:end]) = [];
y([1:2,end-1:end],:) = [];
y(:,[1:2,end-1:end]) = [];

% Dot Center
X = mean(X(~isnan(X)));
Y = mean(Y(~isnan(Y)));

x_old = x;
y_old = y;
x = [];
y = [];
% Group and use appropriate scale
for kk = 1:length(x_pos)
    x = [x x_old+x_pos{kk}-X];
    y = [y y_old+y_pos{kk}-Y];
end
x = repmat(x,[1 1 max_frame]);
y = repmat(y,[1 1 max_frame]);

r = sqrt( x.^2 + y.^2 );

%%%%%%%%%%% Scale r for each time point
%%%% First smooth over NaN gaps in R
% Need to fill in any isolated R values, but not those before/after
% complete edge finding failure
beginning = find(~isnan(R),1);
ending = find(~isnan(R),1,'last');
R_good = R(beginning:ending);
bad_values = find(isnan(R_good));
real_values = find(~isnan(R_good));
x_temp = 1:length(R_good);
R_good_fill = NaN*bad_values;
for kk = 1:length(bad_values)
    R_good_fill(kk) = interp1(x_temp(real_values),R_good(real_values),x_temp(bad_values(kk)));
end
R_good(bad_values) = R_good_fill;
R(beginning:ending) = R_good;

% disp(['     FTLE analysis is only valid for frames ' num2str(beginning) ' to ' num2str(ending) '.'])

% Scale radial locations by the size of the dot
r_scaled = zeros(size(r,1),size(r,2),size(r,3)-1);
for k = beginning:ending
    
    r_scaled(:,:,k) = r(:,:,1)./R(k);
    
end
