% Function to shift multiple PIV files into the frame of the microscope
% stage.  Analysis is limited by the extend of edge detection data.
%
% Usage: [x,y,u,v,r_scaled,beginning, ending] = do_scale_shift_piv(pivC,pivL,pivR,xL,yL,xR,yR,xC,yC,X,Y,R,r_scale,t_scale)

function [x_all,y_all,u,v,r_scaled,beginning, ending] = do_scale_shift_edges(u,v,x,y,x_pos,y_pos,X,Y,R,r_scale,t_scale)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS
% pivC, pivL, pivR  = The three PIV files corresponding to center, left, and right images
% xL,yL,xR,yR,xC,yC = Microscope stage postions (um)
% X,Y               = Dot center locations (um)
% R                 = radius of dot(um)
% r_scale           = microns per pixel
% t_scale           = minutes per frame
%
% OUTPUTS
% x                 = x locations for PIV vectors in dot frame of reference (um)
% y                 = y locations for PIV vectors in dot frame of reference (um)
% u,v               = combined velocity vectors from PIV data (um/min)
% r_scaled          = radial location of PIV vectors, scaled by dot radius
% beginning         = first frame valid for analysis
% ending            = last frame valid for analysis
%
% Created by Rachel Lee, 2014/02/04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Adjust coordinates
for kk = 1:length(x_pos)
    x{kk} = r_scale*x{kk} + x_pos{kk};
    y{kk} = r_scale*y{kk} + y_pos{kk};
end

X = mean(X(~isnan(X)));
Y = mean(Y(~isnan(Y)));

% Don't use more frames than available in PIV data
max_frame = min([length(R) size(u,3) size(v,3)]);
R = R(1:max_frame);

x_all = [];
y_all = x_all;
for kk = 1:length(x_pos)
    x_all = [x_all x{kk}];
    y_all = [y_all y{kk}];
end

x_all = x_all - X;
x_all = repmat(x_all,[1 1 max_frame]);

y_all = y_all - Y;
y_all = repmat(y_all,[1 1 max_frame]);

r = sqrt( x_all.^2 + y_all.^2 );


u = r_scale/t_scale*u;
u = u(:,:,1:max_frame);
v = r_scale/t_scale*v;
v = v(:,:,1:max_frame);

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

disp(['     Radial analysis is only valid for frames ' num2str(beginning) ' to ' num2str(ending) '.'])

% Scale radial locations by the size of the dot
r_scaled = zeros(size(r,1),size(r,2),size(r,3)-1);
for k = beginning:ending
    
    r_scaled(:,:,k) = r(:,:,1)./R(k);
    
end
