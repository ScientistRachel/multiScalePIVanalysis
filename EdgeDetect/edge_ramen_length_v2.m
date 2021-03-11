function edge_ramen_length_v2(directory,imname,edgedat,timestep, diststep, saveimpose,firstframe,lastframe,color_map,high,low)
% 2019/09/27 RML _v2 allows for vertical edges

% Check all the inputs
if ~exist('timestep','var') || isempty(timestep)
    timestep = 1;
end

if ~exist('diststep','var') || isempty(diststep)
    diststep = 1;
end

% % defines a 3 by 'numcolors' matrix that represents the jet colormap.
if ~exist('color_map','var') || isempty(color_map)
    warning('using default colormap')
    color_map = uint8(jet(2^8)*((2^8)-1));
end

% determines the normalization constant of for 'paramater' that will map
% all values to be of the range [-1,1]
if ~exist('high','var') || isempty(high)
    high = 3; %mm, based on rough graphs
end
if ~exist('low','var') || isempty(low)
    low = 1; % mm, based on rough graphs
end
para_size = abs(high);

if ~exist('saveimpose','var') || isempty(saveimpose)
    saveimpose = directory;
end

if ~exist('firstframe','var') || isempty(firstframe)
    firstframe = 1;
end
if ~exist('lastframe','var') || isempty(lastframe)
    lastframe = size(edgedat.points,3);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c = size(edgedat.points,3);
save_warning = [];
if lastframe > c
    disp(['WARNING: ' num2str(lastframe) ' requested, but only ' num2str(c) ' frames available'])
    save_warning = ['- only ' num2str(c) ' frames'];
    lastframe = c;
end

% edgedat.points = edgedat.points(:,:,firstframe:lastframe); %2016-03-02
% this line is unecessary if the for loop includes the frame contraints and
% will cause an error if the firstframe ~= 1.

% measure size of input arrays and imintialize memory for output
c = size(edgedat.points,3);

% 2019/08/27 left and right used to be opposite of this, but I've changed
% the naming convention in sheet_edger_combo_v2 - RML
if strcmp(edgedat.side, 'left') || strcmp(edgedat.side,'top')
    big = max(edgedat.points(:));
elseif strcmp(edgedat.side, 'right') || strcmp(edgedat.side,'bottom')
    big = min(edgedat.points(:)); 
    c = 0;
    diststep = -diststep;
else
    error('Input structure does not contain the correct .side field');
end

if strcmp(edgedat.side, 'left') || strcmp(edgedat.side,'right')
    
    stack_im = zeros(edgedat.imsize(1),big+diststep*c,3,'uint8');

    for i = firstframe:timestep:lastframe

        % eleminates formatting zeros from the 3d element 'i' 'from points_seq'
        points_seq_slice = cutter2(edgedat.points(:,:,i));

        x = diff(double(points_seq_slice(:,2)));
        y = diff(double(points_seq_slice(:,1)));
        L = sum(sqrt(x.^2 + y.^2))*0.65/1000; %mm

        % adds the distance step to the horrizontal postion in each xy points
        % in 'points_seq_slice'
        points_seq_slice(:,2) = points_seq_slice(:,2) + diststep*(c-i);

        g = size(points_seq_slice,1);

        parameter_slice = L*ones(g,1);
        para_color = round(255*((parameter_slice - low)/para_size));
        %Allow for saturation
        para_color(para_color < 1) = 1;
        para_color(para_color > 2^8) = 2^8;

        % stamps the vector in color_map that corresponds to the frame
        % each xy point to each xy point on the image 'stack_im'
        for j = 1:g
            stack_im(points_seq_slice(j,1), points_seq_slice(j,2),1) = color_map(para_color(j), 1);
            stack_im(points_seq_slice(j,1), points_seq_slice(j,2),2) = color_map(para_color(j), 2);
            stack_im(points_seq_slice(j,1), points_seq_slice(j,2),3) = color_map(para_color(j), 3);
        end

    end
    
elseif strcmp(edgedat.side, 'bottom') || strcmp(edgedat.side,'top')
    
    stack_im = zeros(big+diststep*c,edgedat.imsize(2),3,'uint8');

    for i = firstframe:timestep:lastframe

        % eleminates formatting zeros from the 3d element 'i' 'from points_seq'
        points_seq_slice = cutter2(edgedat.points(:,:,i));

        x = diff(double(points_seq_slice(:,2)));
        y = diff(double(points_seq_slice(:,1)));
        L = sum(sqrt(x.^2 + y.^2))*0.65/1000; %mm

        % adds the distance step to the horrizontal postion in each xy points
        % in 'points_seq_slice'
        points_seq_slice(:,1) = points_seq_slice(:,1) + diststep*(c-i);

        g = size(points_seq_slice,1);

        parameter_slice = L*ones(g,1);
        para_color = round(255*((parameter_slice - low)/para_size));
        %Allow for saturation
        para_color(para_color < 1) = 1;
        para_color(para_color > 2^8) = 2^8;

        % stamps the vector in color_map that corresponds to the frame
        % each xy point to each xy point on the image 'stack_im'
        for j = 1:g
            stack_im(points_seq_slice(j,1), points_seq_slice(j,2),1) = color_map(para_color(j), 1);
            stack_im(points_seq_slice(j,1), points_seq_slice(j,2),2) = color_map(para_color(j), 2);
            stack_im(points_seq_slice(j,1), points_seq_slice(j,2),3) = color_map(para_color(j), 3);
        end

    end   
    
end


imwrite(stack_im, [saveimpose imname 'LengthColorStack' save_warning '.bmp']);
