% Runs piv_analytics on a batch of files
% Sends e-mail progress updates and error messages if necessary
%
%%%%%%% WARNING:
% Do not use this batch blindly! Various parameters will need to be 
% adjusted for new data.
%
% Created 2014/02/07 by Rachel Lee
% 2015/06/29 RML updated for MultiC, etc.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc % Clear your workspace and command window

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% USER INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Please set all of the following correctly before running this file
directories = {'J:\2014-M1-M4-Data\2014-04-28-M1M3M4-00011\',...
    'J:\2014-M1-M4-Data\2014-05-05-M1M3M4-0065\',...
    'J:\2014-M1-M4-Data\2014-06-01-M1M3M4-0003\',...
    'J:\2014-M1-M4-Data\2014-09-07 HCl vs Borate-0001\'};

r_scale = 0.65; % microns/pixel
t_scale = 3; % minutes/frame
firstframe = 1; % First frame to analyze
user_lastframe = 200;  % Set lastframe to inf to analyze all frames

disp(datestr(now)) % Display what time the batch began

p = path; %Save the original path to reset later
p2 = genpath('PIV');
p3 = genpath('EdgeDetect');
p4 = genpath('Analytics');
p5 = genpath('FTLE');
p6 = genpath('CoarseGrain');
path(p,[p2,p3,p4,p5,p6]);
clear p2 p3 p4 p5 p6

    
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
        disp(['  Analyzing ' imname ' (' num2str(jj) ' of ' num2str(numel(names_list)) ')'])
        plot_piv_time_hist(directory, imname, r_scale, t_scale, firstframe, user_lastframe,2);
        close all
    end

end

% Display and e-mail completion of batch
disp(' ')
disp('Batch Complete')
disp(' ')

disp(datestr(now)) % Display what time the batch ended

% Reset the path
path(p)