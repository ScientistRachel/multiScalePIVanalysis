% Function that loads PIV, rescales, and puts into good shape
% For coarse_grain_r_and_t_frame

function [u_out,v_out,x_out,y_out,a,final_frame, num_pieces] = set_up_edge_piv_MultiC(directory, imname,final_frame_user,edge_cut,r_scale,t_scale)


% % Load the PIV data
% load([directory imname 'Rdot_piv_filtered.mat'])
% pivR = dot_piv;
% 
% load([directory imname 'Ldot_piv_filtered.mat'])
% pivL = dot_piv;
% 
% load([directory imname 'Cdot_piv.mat'])
% pivC = dot_piv;
% clear dot_piv
% 
% xL = NaN; yL = NaN; xR = NaN; yR = NaN; % just to get rid of warning, not important
% load([directory 'Analytics Data\' imname '_dot_RvsT.mat'])
% 
% %%%% PIV Processing
% % Convert to radial coordinates
% [x,y,u,v,r_scaled,beginning, ending] = do_scale_shift_piv(pivC,pivL,pivR,xL,yL,xR,yR,xC,yC,X,Y,R,r_scale,t_scale);

% Shift PIV locations and set to appropriate scale
list = [dir([directory imname 'C*dot_piv.mat']) ; dir([directory imname '*dot_piv_filtered.mat'])];
u_raw = []; v_raw = u_raw;
x = cell(numel(list),1); y = x;
for kk = 1:numel(list)
    load([directory list(kk).name])
%     u_raw = [u_raw dot_piv.u_fi];
%     v_raw = [v_raw dot_piv.v_fi];
%     x{kk} = dot_piv.x;
%     y{kk} = dot_piv.y;
    c = size(dot_piv.u_fi,3); %2016/03/02 In case PIV is differently sized...
%     disp(c)
    if size(u_raw,3) > c
        u_raw = u_raw(:,:,1:c);
        v_raw = v_raw(:,:,1:c);
    end
    if kk == 1
        b = c;
    else
        b = size(u_raw,3);
    end
    u_raw = [u_raw dot_piv.u_fi(:,:,1:b)];
    v_raw = [v_raw dot_piv.v_fi(:,:,1:b)];
    x{kk} = dot_piv.x;
    y{kk} = dot_piv.y;
end
clear dot_piv

load([directory 'Analytics Data\' imname '_dot_RvsT.mat'])
num_pieces = length(x_pos);

% x,y = um; u,v = um/min; r_scaled = unitless; beginning,ending = frames
[x,y,u,v,r_scaled, beginning, ending] = do_scale_shift_MultiC(u_raw,v_raw,x,y,x_pos,y_pos,X,Y,R,r_scale,t_scale);
[~,a,~] = size(x);

% Need to work on making this more general
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

% Then only keep those values > edge_cut
u(r_scaled < edge_cut) = NaN;
v(r_scaled < edge_cut) = NaN;
x(r_scaled < edge_cut) = NaN;
y(r_scaled < edge_cut) = NaN;

u_out = u;
v_out = v;
y_out = y;
x_out = x;