function run_ftle_time_rPairedTif(directory,r_scale, t_scale, firstframe, lastframe,over_write)

% Look for possible files to analyze
list = dir([directory '*tif']);
names_list = cell(2,numel(list)/2);
for i = 1:numel(list)
    names_list{i} = list(i).name(1:end-4);
end
names_list = names_list'; % Each row is a well now

% Load FTLE values
load([directory 'FTLEall.mat'],'pivfiles','ftle_all')
          
% Loop through the files in this directory
for jj = 1:size(names_list,1)
    imnames = names_list(jj,:);
    imname = [imnames{1} '_' imnames{2}];
    if ~exist([directory 'FTLE Time and R\' imname 'Positive FTLEs.tif'],'file') || over_write

        disp(['    FTLE r/R vs time running on ' imname ' (' num2str(jj) ' of ' num2str(size(names_list,1)) ')'])
        
        %Find the correct ftle files
        list = [dir([directory imnames{1} '*dot_piv_filtered.mat']) ; dir([directory imnames{2} '*dot_piv_filtered.mat'])];
        ftleLRC = cell(length(list),1);
        for kk = 1:length(list)
            found_it = strcmp(pivfiles,list(kk).name);
            ftleLRC{kk} = ftle_all{found_it};
        end

        plot_ftle_time_rPairedTif(directory, imnames, ftleLRC, r_scale, t_scale, firstframe, lastframe)
        close all

    else
        disp(['    FTLE histograms already exist for ' imname ' (' num2str(jj) ' of ' num2str(size(names_list,1)) ')'])
    end
end
