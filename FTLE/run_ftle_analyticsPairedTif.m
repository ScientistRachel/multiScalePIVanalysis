function run_ftle_analyticsPairedTif(directory,r_scale, t_scale, firstframe, lastframe,over_write)

% Look for possible files to analyze
list = dir([directory '*tif']);
names_list = cell(2,numel(list)/2);
for i = 1:numel(list)
    names_list{i} = list(i).name(1:end-4);
end
names_list = names_list'; % Each row is a well now

% Load FTLE values
load([directory 'FTLEall.mat'],'pivfiles','ftle_all')
     
%Set up variables to write into excel
colB = NaN*ones(size(names_list,1),3);
colE = NaN*ones(size(names_list,1),3);
colH = NaN*ones(size(names_list,1),3);
colK = NaN*ones(size(names_list,1),1);
colL = NaN*ones(size(names_list,1),1);
        
% Loop through the files in this directory
for jj = 1:size(names_list,1)
    imnames = names_list(jj,:);
    imname = [imnames{1} '_' imnames{2}];
    if ~exist([directory 'FTLE Analytics\' imname 'FTLEAnalytics.mat'],'file') || over_write

        disp(['    FTLE analytics running on ' imname ' (' num2str(jj) ' of ' num2str(size(names_list,1)) ')'])
        
        %Find the correct ftle files
        list = [dir([directory imnames{1} '*dot_piv_filtered.mat']) ; dir([directory imnames{2} '*dot_piv_filtered.mat'])];
        ftleLRC = cell(length(list),1);
        for kk = 1:length(list)
            found_it = strcmp(pivfiles,list(kk).name);
            ftleLRC{kk} = ftle_all{found_it};
        end

        [colB(jj,:), colE(jj,:), colH(jj,:), colK(jj), colL(jj)] = ftle_analyticsPairedTif(directory, imnames, ftleLRC, r_scale, t_scale, firstframe, lastframe);

    else
        disp(['    FTLE analytics already exist for ' imname ' (' num2str(jj) ' of ' num2str(size(names_list,1)) ')'])
        load([directory 'FTLE Analytics\' imname 'FTLEAnalytics.mat'],...
            'mean_ftle','SEM_ftle', 'per_pos','firstframe','lastframe')
        colB(jj,:) = mean_ftle;
        colE(jj,:) = SEM_ftle;
        colH(jj,:) = per_pos;
        colK(jj) = firstframe;
        colL(jj) = lastframe;
    end
end
% Save collective outputs in .mat format
save([directory 'FTLEanalytics.mat'],'names_list','colB','colE','colH','colK','colL')        