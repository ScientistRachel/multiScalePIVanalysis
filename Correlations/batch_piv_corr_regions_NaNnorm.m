%2016/01/12 Based on batch_piv_corr
% This take a long time to run -- comment out sections of
% piv_corr_r_regions to correlate less combinations of means and regions...

clc
clear
close all

directories = {'F:\2012_Data\2012-03-15-M1-0001\',...
    'F:\2012_Data\2012-06-02-M1-Dots-0005\',...
    'F:\2012_Data\2012-06-16-M1-Dots-0007\',...
    'F:\2012_Data\2012-08-17-M1-Dots-0004\'};

r_scale = 0.65; % microns/pixel
t_scale = 2; % minutes/frame
firstframe = 1; % First frame to analyze
lastframe = 500;  % Set lastframe to inf to analyze all frames
over_write = 0;

my_email = 'rmlee@umd.edu';

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
    list = dir([directory '*-L.zvi']);
    names_list = cell(numel(list),1);
    for i = 1:numel(list)
        names_list{i} = list(i).name(1:end-5);
    end
    names_list = unique(names_list);


    % Loop through the files in this directory
    for jj = 1:numel(names_list)
        imname = names_list{jj};
        disp(['  Analyzing ' imname ' (' num2str(jj) ' of ' num2str(numel(names_list)) ')'])
        piv_corr_r_regions_NaNnorm(directory, imname, r_scale, t_scale, firstframe, lastframe,[],over_write);
    end

end

disp(' ')
disp('Batch Complete')
disp(' ')

close all
disp(datestr(now)) % Display what time the batch ended

%%%%%%%% Only necessary if you are changing the e-mail settings
% These lines change the e-mail address that e-mail is sent FROM
myaddress = 'matlabLCMB@gmail.com';
mypassword = 'GuineaPig!';
% 
setpref('Internet','E_mail',myaddress);
setpref('Internet','SMTP_Server','smtp.gmail.com');
setpref('Internet','SMTP_Username',myaddress);
setpref('Internet','SMTP_Password',mypassword);

props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.socketFactory.class', ...
                  'javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port','465');

sendmail(my_email, 'MATLAB Progress Update', 'Correlation Batch Complete :)');
