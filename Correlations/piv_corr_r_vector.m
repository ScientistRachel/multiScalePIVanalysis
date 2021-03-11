% Function to provide some automated analysis of PIV data.
% Please do not use this function without thinking!
% - Watch the original movie
% - Check the hard coded parameters
% - Watch that the results seem consistent!
%
% Usage: piv_corr_r_vector(directory, imname, r_scale, t_scale, firstframe, lastframe,pos_list)
function piv_corr_r_vector(directory, imname, r_scale, t_scale, firstframe, lastframe,pos_list)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% INPUTS %%%
% directory     = folder of images
% imname        = file names (missing L, R, or C at the end)
% r_scale       = magnification of images in microns per pixel
% t_scale       = time scale of movie, in minutes per frame
% firstframe    = Optional. First frame to analyze. Defaults to 1.
% lastframe     = Optional. Last frame to analyze. Defaults to full image series.
% pos_list      = Specify how to calculation image locations.  Defaults to 1.
%                   0 = Use zvi files to extract position (robust but slow)
%                   1 = Use 'PosList.csv' to find postions (must be in same directory)
%                   2 = 'dot_RvsT.mat' file already exists, use previous analysis
%
% Note: frame range will be limited by edge detection results.
%
%%% OUTPUTS %%%
% Several figures and .mat files are saved to the directory
%
%%% CHANGE LOG %%%
% 2015/07/09 RML - created piv_corr_t based on piv_analyticsMultiC
% 2015/12/14 edgedat file was loaded incorrectly!!! Fixed
% 2016/06/21 RML adapted piv_corr_r_NaNnorm to be a simple vector corr func
% 2016/06/22 now accepts input from either MultiC or old version of code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%% ****PARAMETERS**** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These parameters are often constant between experiments, but should be
% confirmed for all data sets.
param.edge_cut = 0.75; % Only include vectors with r/R greater than this in edge metrics (set to 0 to analyze all)
param.center_cut = 0.5; % Only include vectors with r/R less than this in center metrics (set to 1 to analzye all)
max_dist = 400; % in um, put inf to use as much as possible based on image size

if exist('pos_list','var') && ~isempty(pos_list)
    param.pos_list = pos_list;
elseif exist([directory 'Analytics Data\' imname '_dot_RvsT.mat'],'file')
    param.pos_list = 2;
else
    param.pos_list = 1; % 1 = excel file exists with locations, 0 = find from zvi files (slower)
end

if ~exist('firstframe','var') || isempty(firstframe)
    firstframe = 1;
end

if ~exist('lastframe','var') || isempty(lastframe)
    lastframe = inf;
end

%%%%%%%%%%%%%%%%%%%%%%% Make Folder for Saving %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist([directory 'Correlation Data\'],'file')
    mkdir(directory,'Correlation Data')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% Circle Fitting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if param.pos_list == 2
    
    load([directory 'Analytics Data\' imname '_dot_RvsT.mat'])
    if exist('x_pos','var')
        type = 'new_version';
    else
        type = 'old_version';
    end
    
else

    edgelist = dir([directory '*' imname '*edgedat.mat']);
    DotE = cell(2,1);
    for kk = 1:2
        load([directory edgelist(kk).name])
        DotE{kk} = edgedat;
        clear edgedat
    end

    % All outputs in um
    [X,Y,R,x_pos,y_pos] = do_edgedat_to_circleMultiC(DotE{1}, DotE{2}, pos_list, directory, imname,r_scale);
    clear DotE edgelist
    
end

%%%%%%%%%%%%%%%%%%%%%%% Convert PIV to Radial %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(type,'new_version')
    % Shift PIV locations and set to appropriate scale
    list = [dir([directory imname 'C*dot_piv.mat']) ; dir([directory imname '*dot_piv_filtered.mat'])];
    u_raw = []; v_raw = u_raw;
    x = cell(numel(list),1); y = x;
    for kk = 1:numel(list)
        load([directory list(kk).name])
        c = size(dot_piv.u_fi,3); %2016/03/02 In case PIV is differently sized...
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

    % x,y = um; u,v = um/min; r_scaled = unitless; beginning,ending = frames
    [x,y,u,v,~, beginning, ending] = do_scale_shift_MultiC(u_raw,v_raw,x,y,x_pos,y_pos,X,Y,R,r_scale,t_scale);
elseif strcmp(type,'old_version')
    % Shift PIV locations and set to appropriate scale
    load([directory imname 'Cdot_piv.mat'])
    dot_piv.directory = directory;
    dot_piv.imname = [imname 'C'];
    pivC = dot_piv;
    clear dot_piv

    load([directory imname 'Ldot_piv_filtered.mat'])
    dot_piv.directory = directory;
    dot_piv.imname = [imname 'L'];
    pivL = dot_piv;
    clear dot_piv

    load([directory imname 'Rdot_piv_filtered.mat'])
    dot_piv.directory = directory;
    dot_piv.imname = [imname 'R'];
    pivR = dot_piv;
    clear dot_piv

    % x,y = um; u,v = um/min; r_scaled = unitless; beginning,ending = frames
    [x,y,u,v,~, beginning, ending] = do_scale_shift_piv(pivC,pivL,pivR,xL,yL,xR,yR,xC,yC,X,Y,R,r_scale,t_scale);
else
    error('Something is wrong with your dot_RvsT.mat file')
end

% Update frame range to not exceed rational ranges
firstframe = max(firstframe,beginning);
lastframe = min(lastframe,ending);

%%%%%%%%%%%%%%%%%%%%%%%%% Useful quantities %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u = u(:,:,firstframe:lastframe);
v = v(:,:,firstframe:lastframe);
x = x(:,:,firstframe:lastframe);
y = y(:,:,firstframe:lastframe);

%%%%%%%%%%%%%%%%%%%%%%%%% Time correlations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
savedir = [directory 'Correlation Data\' imname '_vector_radial_corr'];
% correlatePIV_r_vector_NaNnorm(u,v,x,y,t_scale,r_scale,0,lastframe,max_dist,savedir,'Velocity Vector Autocorrelation')
correlatePIV_r_vector_NaNnorm(u,v,x,y,t_scale,0,lastframe,max_dist,savedir,'Velocity Vector Autocorrelation')

