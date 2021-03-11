% Function to shift multiple FTLE files into the frame of the microscope
% stage.  Analysis is limited by the extend of edge detection data.
%
% Usage: [ftles,r_scaled,beginning, ending] = do_scale_shift_ftle(ftleLRC,piv,xL,yL,xR,yR,xC,yC,X,Y,R,r_scale,t_scale)

function [ftles,r_scaled,beginning, ending] = do_scale_shift_ftle(ftleLRC,piv,xL,yL,xR,yR,xC,yC,X,Y,R,r_scale,t_scale)

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Remove edge problems from FTLE data
ftleL = ftleLRC{1}*60/t_scale; %60/t_scale puts it in hr^-1
ftleR = ftleLRC{2}*60/t_scale;
ftleC = ftleLRC{3}*60/t_scale;

ftleL(1:2,:,:) = [];
ftleL(end-1:end,:,:) = [];
ftleL(:,1:2,:) = [];
ftleL(:,end-1:end,:) = [];
ftleL(isinf(ftleL)) = NaN;

ftleR(1:2,:,:) = [];
ftleR(end-1:end,:,:) = [];
ftleR(:,1:2,:) = [];
ftleR(:,end-1:end,:) = [];
ftleR(isinf(ftleR)) = NaN;

ftleC(1:2,:,:) = [];
ftleC(end-1:end,:,:) = [];
ftleC(:,1:2,:) = [];
ftleC(:,end-1:end,:) = [];
ftleC(isinf(ftleC)) = NaN;

ftles = [ftleL ftleC ftleR]; 

% Don't use more frames than available in FTLE data
max_frame = min([size(ftleL,3),size(ftleR,3),size(ftleC,3),length(R)]);
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

% Group and use appropriate scale
x = [(x+xL) (x+xC) (x+xR)] - X;
x = repmat(x,[1 1 max_frame]);

y = [(y+yL) (y+yC) (y+yR)] - Y;
y = repmat(y,[1 1 max_frame]);

r = sqrt( x.^2 + y.^2 );

%%%%%%%%%%% Scale r for each time point
%%%% First smooth over NaN gaps in R
% Need to fill in any isolated R values, but not those before/after
% complete edge finding failure
beginning = find(~isnan(R),1);
ending = find(~isnan(R),1,'last');
R_good = R(beginning:ending);
real_values = find(~isnan(R_good));
x_temp = 1:length(R_good);
R_good_fill = interp1(x_temp(real_values),R_good(real_values),x_temp(~real_values));
R_good(~real_values) = R_good_fill;
R(beginning:ending) = R_good;

disp(['     FTLE analysis is only valid for frames ' num2str(beginning) ' to ' num2str(ending) '.'])

% Scale radial locations by the size of the dot
r_scaled = zeros(size(r,1),size(r,2),size(r,3)-1);
for k = beginning:ending
    
    r_scaled(:,:,k) = r(:,:,1)./R(k);
    
end
