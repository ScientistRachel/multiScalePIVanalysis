function basic_analyticsMultiC(directory,r_scale, t_scale, firstframe, user_lastframe,data_type,over_write)

% Look for possible files to analyze
list = dir([directory '*-R.zvi']);
names_list = cell(numel(list),1);
for i = 1:numel(list)
    names_list{i} = list(i).name(1:end-5);
end
% names_list = unique(names_list); 
% if strcmp(names_list{end}(end-4:end),'repea') %delete the jitter file
%     names_list = names_list(1:end-1);
% end
     
%Set up variables to write into excel
colB = NaN*ones(numel(names_list),3);
colE = NaN*ones(numel(names_list),3);
colH = NaN*ones(numel(names_list),3);
colK = NaN*ones(numel(names_list),1);
colL = NaN*ones(numel(names_list),1);
colM = NaN*ones(numel(names_list),1);
        
% Loop through the files in this directory
for jj = 1:numel(names_list)
    imname = names_list{jj};
    if ~exist([directory 'Analytics Data\' imname 'Analytics.mat'],'file') || over_write
        %Note: previous analysis checked for files here, but this
        %should not be a problem with the main_batch
            disp(['    Analyzing ' imname ' (' num2str(jj) ' of ' num2str(numel(names_list)) ')'])
            [colB(jj,:), colE(jj,:), colH(jj,:), colK(jj), colL(jj), colM(jj)] = piv_analyticsMultiC(directory, imname, r_scale, t_scale, firstframe, user_lastframe,data_type);

    else
        disp(['    Analysis already exists for ' imname ' (' num2str(jj) ' of ' num2str(numel(names_list)) ')'])
        load([directory 'Analytics Data\' imname 'Analytics.mat'])
        colB(jj,:) = ang_dev;
        colE(jj,:) = mean_speed;
        colH(jj,:) = SEM_speed;
        colK(jj) = deltaR;
        colL(jj) = firstframe;
        colM(jj) = lastframe;
    end            
end


% Create a useful file name for the excel file
C = strsplit(directory,'\');
filename = [directory C{end-1} '_PIV_Analytics_' datestr(now,'yymmdd') '_' num2str(firstframe) '-' num2str(round(mean(colM))) '.xlsx']; 
% Put outputs in excel file....
xlswrite(filename,{'File Name','Angular Deviation','','','Mean Speed (um/min)','','','SEM Speed (um/min)','','','Delta R (um)','First Frame Analyzed','Last Frame Analyzed'},'Sheet1','A1')
xlswrite(filename,{'','Total','Edge','Center','Total','Edge','Center','Total','Edge','Center'},'Sheet1','A2')
xlswrite(filename,names_list,'Sheet1','A3')
xlswrite(filename,colB,'Sheet1','B3')
xlswrite(filename,colE,'Sheet1','E3')
xlswrite(filename,colH,'Sheet1','H3')
xlswrite(filename,colK,'Sheet1','K3')
xlswrite(filename,colL,'Sheet1','L3')
xlswrite(filename,colM,'Sheet1','M3')
        
% Save collective outputs in .mat format
save([directory 'PIVanalytics.mat'],'names_list','colB','colE','colH','colK','colL','colM')        