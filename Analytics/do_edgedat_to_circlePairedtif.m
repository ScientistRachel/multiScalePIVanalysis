% do_edgedat_to_cirle fits the information in edgedat.mat to a circle
% Usage: [X,Y,R,xL,xR,yL,yR,xC,yC] = do_edgedat_to_circle(DotL, DotR, pos_list, directory, imname,r_scale)

function [X,Y,R,x_pos,y_pos] = do_edgedat_to_circlePairedtif(Dot1, Dot2, pos_list, directory, imnames,r_scale)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% INPUTS %%%%%%
% DotL      = edgedat data for left edge
% DotR      = edgedat data for right edge
% pos_list  = 1 -> Use excel file to find image positions
%             0 -> Use zvi file to find image positions (slower)
% directory = main directory for data set, includes edgedat.mat
% imnames    = names of the two tif  files of interest
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
% 2016/01/06 RML - when reading from zvi files, get file names a simpler way
% 2019/08/29 RML - Pairedtif version to parse spinning disk data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
list = [dir([directory imnames{1} '*dot_piv_filtered.mat']) ; dir([directory imnames{2} '*dot_piv_filtered.mat'])];
piece = cell(numel(list),1);
for kk = 1:numel(list)
    piece{kk} = list(kk).name(1:end-20); %Edge pieces
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Image Locations %%%%%%%%%%%%%%%%%%%%%%%%
if pos_list %Use an excel file
    
    % Assume a standardized file in the same directory (may need to update)
%     poslist_file = [directory 'PosList.csv'];
%     poslist_file = dir([directory '*.csv']);
    
    % assuming the files were set up with convertMVD2, this is correct
    C = strsplit(directory,'\');    
    poslist_file = [C{1} '\' C{2} '\' C{3} '\' 'PosList.xls'];

    % Load in position list from microscope (expecting .csv in standard format)
    [~, ~, M] = xlsread(poslist_file);
    
    x_pos = cell(length(piece),1);
    y_pos = x_pos;
    for kk = 1:length(piece)
        found = 0;
        k = 1; %start at the top
        while ~found && k <= size(M,1)
            found = strcmp(M{k,1},piece{kk});
            k = k+1;
        end
        if found
            x_pos{kk} = M{k-1,2};
            y_pos{kk} = M{k-1,3};
        else
            error(['Microscope stage position not found for ' piece{kk}])
        end        
        
    end

    
else % Find the stage locations using the zvi files -- should more generally work, but is slower
    
    error('This code is not set up to read mvd2 metadata')
    
%     for kk = 1:length(piece)
%         filename = [directory imname piece{kk} '.zvi'];
%         [xtemp,ytemp] = get_zvi_location(filename);
%         x_pos{kk} = xtemp;
%         y_pos{kk} = ytemp;
%         
%         if strcmp(piece{kk}(end),'L')
%             xL = xtemp;
%             yL = ytemp;
%         elseif strcmp(piece{kk}(end),'R')
%             xR = xtemp;
%             yR = ytemp;
%         end
%         
%     end

%     for kk = 1:length(list)
%         if strcmp(list(kk).name(end-19:end),'dot_piv_filtered.mat')
%             filename_part = list(kk).name(1:end-20); %Edge pieces
%         else
%             filename_part = list(kk).name(1:end-11); %Center pieces
%         end
%         filename = [directory filename_part '.zvi'];
%         [xtemp,ytemp] = get_zvi_location(filename);
%         x_pos{kk} = xtemp;
%         y_pos{kk} = ytemp;
%         
%         if strcmp(piece{kk}(end),'L')
%             xL = xtemp;
%             yL = ytemp;
%         elseif strcmp(piece{kk}(end),'R')
%             xR = xtemp;
%             yR = ytemp;
%         end
%         
%     end
%     
%     disp(' ')
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Frames to analyze %%%%%%%%%%%%%%%%%%%%%%%%
if isfield(Dot1,'firstframe')
    start1 = Dot1.firstframe;
else
    start1 = 1;
end

if isfield(Dot2,'firstframe')
    start2 = Dot2.firstframe;
else
    start2 = 1;
end

frames1 = size(Dot1.points,3) + start1 - 1;
frames2 = size(Dot2.points,3) + start2 - 1;

minstart = max(start1,start2);
maxend = min(frames1,frames2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Use circle fitting %%%%%%%%%%%%%%%%%%%%%%%%

% Preallocate storage
X = NaN*zeros(maxend,1);
Y = NaN*zeros(maxend,1);
R = NaN*zeros(maxend,1);

% Loop through applicable frames
for k = minstart:maxend
        
    sliceL = r_scale*double(cutter2(Dot1.points(:,:,k-start1+1)));
    sliceR = r_scale*double(cutter2(Dot2.points(:,:,k-start2+1)));

    if ~isempty(sliceL) && ~isempty(sliceR) %If no edge found, no center found
        
        % Combine the two edges in the microscope reference frame
        points = [(sliceL(:,1) + y_pos{1}) , (sliceL(:,2) + x_pos{1}) ; (sliceR(:,1) + y_pos{2}) , (sliceR(:,2) + x_pos{2})];

        % Use circfit to fit the edges to a circle and calculate output
        [X(k), Y(k), R(k)] = circfit(points(:,2),points(:,1));
        
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% SAVE DATA %%%%%%%%%%%%%%%%%%%%%%%%%%

save([directory 'Analytics Data\' imnames{1} '_' imnames{2} '_dot_RvsT.mat'],'X','Y','R','x_pos','y_pos')
