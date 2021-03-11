function run_coarse_grain_pairedTif(directory,user_lastframe,r_scale,t_scale,overwrite,parpoolName)

if ~exist('parpoolName','var') || isempty(parpoolName)
    parpoolName = 'local';
end

% Runs coarse_grain_r_and_t on a batch of files
% 2016/04/25 Changed time coarse grain to have sqrt jitter timescale
% 2019/08/29 RML 

% Look for possible files to analyze
list = dir([directory '*tif']);
names_list = cell(2,numel(list)/2);
for i = 1:numel(list)
    names_list{i} = list(i).name(1:end-4);
end
names_list = names_list'; % Each row is a well now

%Preallocate storage    
lc = NaN*ones(size(names_list,1),1);
tc = NaN*ones(size(names_list,1),1);
fitR = NaN*ones(size(names_list,1),3);
fitT = NaN*ones(size(names_list,1),4);
               
% Loop through the files in this directory
for jj = 1:size(names_list,1)
    imnames = names_list(jj,:);
    imname = [imnames{1} '_' imnames{2}];

    if ~exist([directory 'CoarseGrain\' imname '_coarse_data.mat'],'file') || overwrite
        
        if isempty(gcp('nocreate'))
            parpool(parpoolName);
        end
        
        [lc(jj), tc(jj),fitR(jj,:),fitT(jj,:)] = coarse_grain_pairedTif_sqrt(directory,imnames,user_lastframe,r_scale,t_scale);
    else
        load([directory 'CoarseGrain\' imname '_coarse_data.mat'],'l_c','t_c','pR','pT')
        lc(jj) = l_c;
        tc(jj) = t_c;
        fitR(jj,:) = pR;
        fitT(jj,:) = pT;
    end
end

% Save collective outputs in .mat format
save([directory 'CoarseGrain\coarse_scales.mat'],'names_list','lc','tc','fitR','fitT')