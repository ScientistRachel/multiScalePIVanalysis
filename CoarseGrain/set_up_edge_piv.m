% Function that loads PIV, rescales, and puts into good shape
% For coarse_grain_r_and_t_frame
% 2015-06-01 RML - THIS ISN'T COMPATIBLE WITH MUTLIC
function [u_out,v_out,x_out,y_out,a,final_frame] = set_up_edge_piv(directory, imname,final_frame_user,edge_cut,r_scale,t_scale)


% Load the PIV data
load([directory imname 'Rdot_piv_filtered.mat'])
pivR = dot_piv;

load([directory imname 'Ldot_piv_filtered.mat'])
pivL = dot_piv;

load([directory imname 'Cdot_piv.mat'])
pivC = dot_piv;
clear dot_piv

xL = NaN; yL = NaN; xR = NaN; yR = NaN; % just to get rid of warning, not important
load([directory 'Analytics Data\' imname '_dot_RvsT.mat'])

%%%% PIV Processing
% Convert to radial coordinates
[x,y,u,v,r_scaled,beginning, ending] = do_scale_shift_piv(pivC,pivL,pivR,xL,yL,xR,yR,xC,yC,X,Y,R,r_scale,t_scale);
if beginning ~=1
    error('need to adapt for data that does not start at frame 1')
end
% Don't keep unecessary frames
final_frame = min(ending,final_frame_user);
x = x(:,:,1:final_frame);
y = y(:,:,1:final_frame);
u = u(:,:,1:final_frame);
v = v(:,:,1:final_frame);
r_scaled = r_scaled(:,:,1:final_frame);

% Only keep 'edges' -- first delete center to decrease processing time
[~,a,~] = size(x);
x(:,(a/3+1):2/3*a,:) = [];
y(:,(a/3+1):2/3*a,:) = [];
u(:,(a/3+1):2/3*a,:) = [];
v(:,(a/3+1):2/3*a,:) = [];
r_scaled(:,(a/3+1):2/3*a,:) = [];
% Then only keep those values > edge_cut
u(r_scaled < edge_cut) = NaN;
v(r_scaled < edge_cut) = NaN;
x(r_scaled < edge_cut) = NaN;
y(r_scaled < edge_cut) = NaN;

u_out = u;
v_out = v;
y_out = y;
x_out = x;