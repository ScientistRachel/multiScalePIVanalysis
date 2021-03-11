% do_edgedat_to_cirle fits the information in edgedat.mat to a circle
% Usage: [X,Y,R,xL,xR,yL,yR,xC,yC] = do_edgedat_to_circle(DotL, DotR, pos_list, directory, imname,r_scale)

function [X,Y,R,xL,xR,yL,yR,xC,yC] = do_edgedat_to_circle(DotL, DotR, pos_list, directory, imname,r_scale)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% INPUTS %%%%%%
% DotL      = edgedat data for left edge
% DotR      = edgedat data for right edge
% pos_list  = 1 -> Use excel file to find image positions
%             0 -> Use zvi file to find image positions (slower)
% directory = main directory for data set, includes edgedat.mat
% imname    = zvi file name -- excluding L or R
% rscale    = um per pixel

%%%%%% OUTPUTS %%%%%%
% All outputs are in um!
% X      = x position of the dot center
% Y      = y position of the dot center
%               -- both in the microscope reference frame
% R      = radius of the dot in pixels
% xL/yL  = x/y location of left image
% xR/yR  = x/y location of right image
% xC/yC  = x/y location of center image
%
%
% Created by Rachel Lee, 2014/02/04
% 2014/04/08 RML - first .csv file in the directory is used to get
%                   positions (changed from requiring naming convention of 
%                   PosList.csv)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Image Locations %%%%%%%%%%%%%%%%%%%%%%%%
if pos_list %Use an excel file
    
    % Assume a standardized file in the same directory (may need to update)
%     poslist_file = [directory 'PosList.csv'];
    poslist_file = dir([directory '*.csv']);
    poslist_file = [directory poslist_file(1).name];

    % Load in position list from microscope (expecting .csv in standard format)
    [~, ~, M] = xlsread(poslist_file);
    
    % Get appropriate name length
    found = M{10,1};
    cut = length(found) - 2;

    % Position for L edge
    found = 0;
    k = 10;  % should be junk at top of excel file, so start with k = 10
    while ~found && k <= size(M,1)
        found = strcmp(M{k,1},[imname(end-cut:end) 'L']);
        k = k+1;
    end
    if found
        xL = M{k-1,2};
        yL = M{k-1,3};
    else
        error(['Microscope stage position not found for ' imname 'L.'])
    end

    % Position for R edge
    found = 0;
    k = 10;  % should be junk at top of excel file, so start with k = 10
    while ~found && k <= size(M,1)
        found = strcmp(M{k,1},[imname(end-cut:end) 'R']);
        k = k+1;
    end
    if found
        xR = M{k-1,2};
        yR = M{k-1,3};
    else
        error(['Microscope stage position not found for ' imname 'R.'])
    end

    % Position for Center
    found = 0;
    k = 10;  % should be junk at top of excel file, so start with k = 10
    while ~found && k <= size(M,1)
        found = strcmp(M{k,1},[imname(end-cut:end) 'C']);
        k = k+1;
    end
    if found
        xC = M{k-1,2};
        yC = M{k-1,3};
    else
        error(['Microscope stage position not found for ' imname 'C.'])
    end
    
else % Find the stage locations using the zvi files -- should more generally work, but is slower
    
    filename = [directory imname 'L.zvi'];
    [xL, yL] = get_zvi_location(filename);
    
    filename = [directory imname 'R.zvi'];
    [xR, yR] = get_zvi_location(filename);
    
    filename = [directory imname 'C.zvi'];
    [xC, yC] = get_zvi_location(filename);
    
    disp(' ')
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Frames to analyze %%%%%%%%%%%%%%%%%%%%%%%%
if isfield(DotL,'firstframe')
    startL = DotL.firstframe;
else
    startL = 1;
end

if isfield(DotR,'firstframe')
    startR = DotR.firstframe;
else
    startR = 1;
end

framesL = size(DotL.points,3) + startL - 1;
framesR = size(DotR.points,3) + startR - 1;

minstart = max(startL,startR);
maxend = min(framesL,framesR);

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Use circle fitting %%%%%%%%%%%%%%%%%%%%%%%%

% Preallocate storage
X = NaN*zeros(maxend,1);
Y = NaN*zeros(maxend,1);
R = NaN*zeros(maxend,1);

% Loop through applicable frames
for k = minstart:maxend
        
    sliceL = r_scale*double(cutter2(DotL.points(:,:,k-startL+1)));
    sliceR = r_scale*double(cutter2(DotR.points(:,:,k-startR+1)));

    if ~isempty(sliceL) && ~isempty(sliceR) %If no edge found, no center found
        
        % Combine the two edges in the microscope reference frame
        points = [(sliceL(:,1) + yL) , (sliceL(:,2) + xL) ; (sliceR(:,1) + yR) , (sliceR(:,2) + xR)];

        % Use circfit to fit the edges to a circle and calculate output
        [X(k), Y(k), R(k)] = circfit(points(:,2),points(:,1));
        
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% SAVE DATA %%%%%%%%%%%%%%%%%%%%%%%%%%

save([directory 'Analytics Data\' imname '_dot_RvsT.mat'],'X','Y','R','xL','xR','yL','yR','xC','yC')
