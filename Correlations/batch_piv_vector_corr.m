clc
clear
close all

%%%%%%%%%%%%%%%%%%%%%%% USER PARAMETERS
directories ={'H:\2016-05 LPA Data Analysis\Ki16425 expt\140507\'};

r_scale = 0.65; % microns/pixel
t_scale = 3; % minutes/frame
firstframe = 1; % First frame to analyze
lastframe = 180;  % Set lastframe to inf to analyze all frames
over_write = 0;
%%%%%%%%%%%%%%%%%%%%%%%
% Change Log
% File organized by RML 2016/06/21

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
        if ~exist([directory 'Correlation Data\' imname '_vth_radial_corrNaNnorm.mat'],'file') || over_write
            disp(['  Analyzing ' imname ' (' num2str(jj) ' of ' num2str(numel(names_list)) ')'])
            piv_corr_r_vector(directory, imname, r_scale, t_scale, firstframe, lastframe); % 2016/06/22 taking ~5 min to run per imname
        end
    end

end

disp(' ')
disp('Batch Complete')
disp(' ')

disp(datestr(now)) % Display what time the batch ended
close all

