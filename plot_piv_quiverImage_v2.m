% Create plots with quiver overlays on speed - uses dot_piv structure
% Usage: plot_piv_quiverImage(dot_piv, [qscale], [savedir],[color_arrow])

function plot_piv_quiverImage(dot_piv, qscale,savedir,color_arrow)

%%%% Inputs
% dot_piv:      structure output of dot_matpiv code suites
% qscale:       (Optional) Manipulates length of vectors. Default = 12. 
%                       Zero results in no scaling
% savedir:      (Optional) Location to save images.  Default is a folder
%                   created in the PIV directory with name dot_piv.imname
% color_arrow:  (Optional) Color of the arrows.  Defaults to 'g'.
%                   Accepts RGB or MATLAB color strings
%
%%%% Outputs
% imname_quiverImage_0000.tif: Image sequence saved in savedir
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rachel Lee 2016/10/07 Based on plot_piv_pcolor
% Initially tested in 2014b -- figure commands might not be compatible with
% other versions
% WARNING: NOT ALL FILE TYPES HAVE BEEN TESTED YET
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Take care of optional inputs
% Quiver scaling:
qscale_default = 12;
if ~exist('qscale','var') || isempty(qscale)
    qscale = qscale_default;
end
% Save directory:
if ~exist('savedir','var') || isempty(savedir)
    % Check for correct directory formatting
    if ~strmcp(dot_piv.directory(end),filesep)
        dot_piv.directory = [dot_piv.directory filesep];
    end
    % Set up a directory based on image name
    if ~isempty(dot_piv.imname)
        savedir = [dot_piv.directory 'quiverImage_' dot_piv.imname filesep];
    else
        savedir = [dot_piv.directory 'quiverImage' filesep];
    end
end
color_arrow_default = 'g';
if ~exist('color_arrow','var') || isempty(color_arrow)
    color_arrow = color_arrow_default;
end

% If the save folder doesn't exist, create it
if ~exist(savedir,'file')
    mkdir(savedir)
end

%%%%%%%%

figure(1)

% Load in images based on type
if strcmp(dot_piv.type,'zvi') || strcmp(dot_piv.type,'czi')
    
    top_cut = 50;% zvi files often have problems....
    
    % Load the images
    h = waitbar(0,'Opening Bio-Formats Reader...');
    images = bfopen4PIV([dot_piv.directory dot_piv.imname '.' dot_piv.type],h);
    images = images{1};
    delete(h); % Delete the waitbar
    
    % Sort through images to only keep desired set
    in_channel = dot_piv.channel:dot_piv.num_channel:size(images,1);
    images = images(in_channel,1);
    
    firstframe = dot_piv.firstframe;
    lastframe = min([dot_piv.lastframe,length(images),size(dot_piv.u_fi,3)]);

    % Create frames from each image
    for p = firstframe:lastframe
        
        % Show the image
        imshow(imadjust(images{p}(top_cut+1:end,:)))

        % quiver plots the flow field arrows
        hold on
        j = quiver(dot_piv.x,dot_piv.y,dot_piv.u_fi(:,:,p),dot_piv.v_fi(:,:,p),0,'Color',color_arrow);
        if qscale ~= 0
            hU = get(j,'UData') ;
            hV = get(j,'VData') ;
            set(j,'UData',qscale*hU,'VData',qscale*hV)
        end
        hold off    

        % Save the figure
        saveas(gcf,[savedir dot_piv.imname '_quiverImage_PIV' num2str(p,'%04u') '_Frame' num2str(p,'%04u') '.tif'],'tif')
        
        % Show the next image but same quiver
        imshow(imadjust(images{p+1}(top_cut+1:end,:)))

        % quiver plots the flow field arrows
        hold on
        j = quiver(dot_piv.x,dot_piv.y,dot_piv.u_fi(:,:,p),dot_piv.v_fi(:,:,p),0,'Color',color_arrow);
        if qscale ~= 0
            hU = get(j,'UData') ;
            hV = get(j,'VData') ;
            set(j,'UData',qscale*hU,'VData',qscale*hV)
        end
        hold off    

        % Save the figure
        saveas(gcf,[savedir dot_piv.imname '_quiverImage_PIV' num2str(p,'%04u') '_Frame' num2str(p+1,'%04u') '.tif'],'tif')

    end
    
elseif strcmp(dot_piv.type,'multiff')  %%% Multipage tif files
    
    % Find correct file names
    image_info = imfinfo([dot_piv.directory dot_piv.imname '.tif']);
       
    % Pick reasonable last frame
    firstframe = dot_piv.firstframe;
    lastframe = min([dot_piv.lastframe length(image_info)/dot_piv.num_channel]);  
    % Sort through images to only keep desired set
    in_channel = firstframe*dot_piv.channel:dot_piv.num_channel:lastframe*dot_piv.num_channel;
    
    % Create frames from each image
    for p = firstframe:lastframe
        
        % Read the image
        image = imread([dot_piv.directory dot_piv.imname '.tif'],in_channel(p));
        
        % Show the image
        imshow(imadjust(image))

        % quiver plots the flow field arrows
        hold on
        j = quiver(dot_piv.x,dot_piv.y,dot_piv.u_fi(:,:,p),dot_piv.v_fi(:,:,p),0,'Color',color_arrow);
        if qscale ~= 0
            hU = get(j,'UData') ;
            hV = get(j,'VData') ;
            set(j,'UData',qscale*hU,'VData',qscale*hV)
        end
        hold off    

        % Save the figure
        saveas(gcf,[savedir dot_piv.imname '_quiverImage_PIV' num2str(p,'%04u') '_Frame' num2str(p,'%04u') '.tif'],'tif')
        
        % Show the next image but same quiver
        image = imread([dot_piv.directory dot_piv.imname '.tif'],in_channel(p+1));
        
        % Show the image
        imshow(imadjust(image))

        % quiver plots the flow field arrows
        hold on
        j = quiver(dot_piv.x,dot_piv.y,dot_piv.u_fi(:,:,p),dot_piv.v_fi(:,:,p),0,'Color',color_arrow);
        if qscale ~= 0
            hU = get(j,'UData') ;
            hV = get(j,'VData') ;
            set(j,'UData',qscale*hU,'VData',qscale*hV)
        end
        hold off    

        % Save the figure
        saveas(gcf,[savedir dot_piv.imname '_quiverImage_PIV' num2str(p,'%04u') '_Frame' num2str(p+1,'%04u') '.tif'],'tif')

    end
    
elseif strcmp(dot_piv.type,'mat')  %%% .mat file -- will assume first variable is a 3D matrix of images!
    
    m = matfile([directory imname '.mat']);
    varlist = who(m);
    images = m.(varlist{1}); %first attempt, try loading all images in at once -- might be better to load as needed??
    
    % Pick reasonable last frame
    firstframe = dot_piv.firstframe;
    lastframe = min([dot_piv.lastframe size(images,3)/dot_piv.num_channel]);  
    % Sort through images to only keep desired set
    in_channel = firstframe*dot_piv.channel:dot_piv.num_channel:lastframe*dot_piv.num_channel;
    
    % Create frames from each image
    for p = firstframe:lastframe
        
        % Show the image
        imshow(imadjust(images(:,:,in_channel(p))))

        % quiver plots the flow field arrows
        hold on
        j = quiver(dot_piv.x,dot_piv.y,dot_piv.u_fi(:,:,p),dot_piv.v_fi(:,:,p),0,'Color',color_arrow);
        if qscale ~= 0
            hU = get(j,'UData') ;
            hV = get(j,'VData') ;
            set(j,'UData',qscale*hU,'VData',qscale*hV)
        end
        hold off    

        % Save the figure
        saveas(gcf,[savedir dot_piv.imname '_quiverImage_PIV' num2str(p,'%04u') '_Frame' num2str(p,'%04u') '.tif'],'tif')
        
        % Show the next image same quiver
        imshow(imadjust(images(:,:,in_channel(p+1))))

        % quiver plots the flow field arrows
        hold on
        j = quiver(dot_piv.x,dot_piv.y,dot_piv.u_fi(:,:,p),dot_piv.v_fi(:,:,p),0,'Color',color_arrow);
        if qscale ~= 0
            hU = get(j,'UData') ;
            hV = get(j,'VData') ;
            set(j,'UData',qscale*hU,'VData',qscale*hV)
        end
        hold off    

        % Save the figure
        saveas(gcf,[savedir dot_piv.imname '_quiverImage_PIV' num2str(p,'%04u') '_Frame' num2str(p+1,'%04u') '.tif'],'tif')
    end
    
else  %%% Image sequences that are not in .zvi format, individual images
    
    % Find correct file names
    list = dir([directory imname '*' imname2 '.' type]);
    
    % Pick reasonable last frame
    firstframe = dot_piv.firstframe;
    lastframe = min([dot_piv.lastframe length(list)]);  
    % Sort through images to only keep desired set
    in_channel = firstframe*dot_piv.channel:dot_piv.num_channel:lastframe*dot_piv.num_channel;
    list = list(in_channel);
    
    % Create frames from each image
    for p = firstframe:lastframe
        
        % Read the image
        image = imread([dot_piv.directory list(p).name],dot_piv.type);
        
        % Show the image
        imshow(imadjust(image))

        % quiver plots the flow field arrows
        hold on
        j = quiver(dot_piv.x,dot_piv.y,dot_piv.u_fi(:,:,p),dot_piv.v_fi(:,:,p),0,'Color',color_arrow);
        if qscale ~= 0
            hU = get(j,'UData') ;
            hV = get(j,'VData') ;
            set(j,'UData',qscale*hU,'VData',qscale*hV)
        end
        hold off    

        % Save the figure
        saveas(gcf,[savedir dot_piv.imname '_quiverImage_PIV' num2str(p,'%04u') '_Frame' num2str(p,'%04u') '.tif'],'tif')
        
        % Read the image
        image = imread([dot_piv.directory list(p+1).name],dot_piv.type);
        
        % Show the image
        imshow(imadjust(image))

        % quiver plots the flow field arrows
        hold on
        j = quiver(dot_piv.x,dot_piv.y,dot_piv.u_fi(:,:,p),dot_piv.v_fi(:,:,p),0,'Color',color_arrow);
        if qscale ~= 0
            hU = get(j,'UData') ;
            hV = get(j,'VData') ;
            set(j,'UData',qscale*hU,'VData',qscale*hV)
        end
        hold off    

        % Save the figure
        saveas(gcf,[savedir dot_piv.imname '_quiverImage_PIV' num2str(p,'%04u') '_Frame' num2str(p+1,'%04u') '.tif'],'tif')

    end
end
