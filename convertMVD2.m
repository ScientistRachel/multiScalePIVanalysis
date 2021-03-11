clc
clear
close all

%%%%%% !!!!!!!! EXCEL NAMING CONVENTION MAKES ASSUMPTIONS !!!!!!!!!!!
% to dos: display imname, add overwrite option

directories = {'E:\2020-12-15_10Amutants_serumStarvation\2020-12-15_10Amutants_serumStarvation_Migration\'};

r_scale = 0.582; % microns/pixel
t_scale = 3; % minutes/frame
scalebar = 100; % um

for kk = 1:numel(directories)  % Run each directory, then each analysis
        
        % Format the directory string correctly
        if ~strcmp(directories{kk}(end),'\')
            directory = [directories{kk} '\'];
        else
            directory = directories{kk};
        end
        disp(['Working on the folder: ' directory ' (' num2str(kk) ' of ' num2str(numel(directories)) ')'])
        
        filename = dir([directory '*.mvd2']);
        if length(filename) > 1
            error('Two files found')
        end
        filename = filename(1).name(1:end-5);
        
        savedir1 = [directory 'aviFormat_scalebar' num2str(scalebar) '_HHMMclock' filesep];
        if ~exist(savedir1,'file')
            mkdir(savedir1)
        end
        
        savedir2 = [directory 'tifFormat' filesep];
        if ~exist(savedir2,'file')
            mkdir(savedir2)
        end
        
        %%% Make an excel file
        filenameXLS = [directory 'PosList.xls'];
        
        im = bfopen([directory filename '.mvd2']);
        
        x_pos = NaN*ones(size(im,1),1);
        y_pos = x_pos;
        for jj = 1:size(im,1)

            x_pos(jj) = im{jj,2}.get('X Location');
            y_pos(jj) = im{jj,2}.get('Y Location');
            
            pointName = ['XYpoint' num2str(jj,'%04u')];
            disp(['     ' pointName])
            
            xlswrite(filenameXLS,{pointName},'Sheet1',['A' num2str(jj+9)]) % The + 9 is to start data on row 10 to match other formats
            
            imPoint = im{jj,1}(:,1);
                     
            imwrite(imPoint{1},[savedir2 pointName '.tif'],'tif')
            
            imNow = rescale_image(imPoint{1},8); % 8 bit for movies
            % Scale Bar
            [a,~] = size(imNow);
            imNow = insertShape(imNow,'filledrectangle',[20 a-30 20+scalebar/r_scale 15],'color','white','opacity',1);
            % Clock
            clock = (1-1)*3; % minutes
            h = floor(clock/60);
            m = clock - h*60;
            clock = [num2str(h,'%02u') ':' num2str(m,'%02u')];
            imNow = insertText(imNow,[15 5],clock,'BoxOpacity',0,'TextColor','White','FontSize',60);
            imwrite(imNow,[savedir1 pointName '.tif'],'tif')
            
            % Set up ability to also make an avi file for 
            v = VideoWriter([savedir1 pointName '.mp4'],'MPEG-4');
            v.FrameRate = 20;
            open(v)
            writeVideo(v,imNow)
            
            for ii = 2:length(imPoint)
                
                imwrite(imPoint{ii},[savedir2 pointName '.tif'],'tif','WriteMode','append')
                
                imNow = rescale_image(imPoint{ii},8); % 8 bit for movies
                % Scale Bar
                [a,~] = size(imNow);
                imNow = insertShape(imNow,'filledrectangle',[20 a-30 20+scalebar/r_scale 15],'color','white','opacity',1);
                % Clock
                clock = (ii-1)*3; % minutes
                h = floor(clock/60);
                m = clock - h*60;
                clock = [num2str(h,'%02u') ':' num2str(m,'%02u')];
                imNow = insertText(imNow,[15 5],clock,'BoxOpacity',0,'TextColor','White','FontSize',60);
                imwrite(imNow,[savedir1 pointName '.tif'],'tif','WriteMode','append')
                
                writeVideo(v,imNow)
                
            end            

            close(v)
            
        end

        xlswrite(filenameXLS,x_pos,'Sheet1','B10')
        xlswrite(filenameXLS,y_pos,'Sheet1','C10')
    
end

disp(' ')
disp('Conversion Complete')
