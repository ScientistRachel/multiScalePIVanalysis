%%%%%%% I NEED TO ADD THE TIME CORRELATIONS TO THIS STILL

clc
clear
close all

% tic
% % 
% % load('J:\2014-11-15-EcadshM1\2014-11-15-EcadshM1-0023\2014-11-15-EcadshM1-0023_A2-Ldot_piv_filtered.mat')
% % lastframe = 200;
% % max_dist = 400; % in um, put inf to use as much as possible based on image size
% % r_scale = 0.65;
% % t_scale = 3;
% % PIVtopixel = dot_piv.x(1,2) - dot_piv.x(1,1);
% % 
% % speed = sqrt(dot_piv.u_fi.^2 + dot_piv.v_fi.^2);
% % 
% % correlatePIV_r(speed,t_scale,r_scale,0,lastframe,PIVtopixel,max_dist)
% % 
% % toc
% 
% directory = 'J:\2014-11-26-EcadshM1\2014-11-26-EcadshM1-0005\';
% imname = '2014-11-26-EcadshM1-0005_B1-';
% 
% % directory = 'J:\2015-04-04 M4 Ecadsh\2015-04-04 Ecadsh M4-0008\';
% % imname = '2015-04-04 Ecadsh M4-0008_B1-';
% 
% lastframe = 200;
% r_scale = 0.65;
% t_scale = 3;
% speed = piv_corr_r(directory, imname, r_scale, t_scale, 1, lastframe);
% 
% toc
% 
% %%
% Please set all of the following correctly before running this file
% load('E:\Documents\2015-04 Ecad Story\EcadStory_data_names.mat')
% directories = data_names(:,1);

directories = {'J:\2014-M1-M4-Data\2014-04-28-M1M3M4-00011\',...
    'J:\2014-M1-M4-Data\2014-05-05-M1M3M4-0065\',...
    'J:\2014-M1-M4-Data\2014-06-01-M1M3M4-0003\',...
    'J:\2014-M1-M4-Data\2014-09-07 HCl vs Borate-0001\'};

% load('E:\Documents\2015-08 Long Term Dynamics\Ecad_shRNA_data_names.mat')
% directories = data_names(:,1);

r_scale = 0.65; % microns/pixel
t_scale = 3; % minutes/frame
firstframe = 1; % First frame to analyze
lastframe = 200;  % Set lastframe to inf to analyze all frames
over_write = 0;

disp(datestr(now)) % Display what time the batch began
    
for k = 1:numel(directories)

    % Check directory for proper format
    directory = directories{k};
    if ~strcmp(directory(end),'\');
        directory = [directories{k} '\'];
    end

    % Update user on progress
    disp(['Working on the folder: ' directory ' (' num2str(k) ' of ' num2str(numel(directories)) ')'])

    % Look for possible files to analyze
    list = dir([directory '*L.zvi']);
    names_list = cell(numel(list),1);
    for i = 1:numel(list)
        names_list{i} = list(i).name(1:end-5);
    end
    names_list = unique(names_list);


    % Loop through the files in this directory
    for jj = 1:numel(names_list)
        imname = names_list{jj};
        if ~exist([directory 'Correlation Data\' imname '_vth_radial_corr.mat'],'file') || over_write
            disp(['  Analyzing ' imname ' (' num2str(jj) ' of ' num2str(numel(names_list)) ')'])
            piv_corr_r(directory, imname, r_scale, t_scale, firstframe, lastframe);
        end
    end

end

disp(' ')
disp('Batch Complete')
disp(' ')

disp(datestr(now)) % Display what time the batch ended

