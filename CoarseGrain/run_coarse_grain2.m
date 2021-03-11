function run_coarse_grain2(directory,user_lastframe,r_scale,t_scale)

% Runs coarse_grain_r_and_t on a batch of files


% Look for possible files to analyze
list = dir([directory '*dot_piv_filtered.mat']);
names_list = cell(numel(list),1);
for i = 1:numel(list)
    names_list{i} = list(i).name(1:end-21);
end
names_list = unique(names_list);

%Preallocate storage    
lc = NaN*ones(numel(names_list),1);
tc = NaN*ones(numel(names_list),1);
fitR = NaN*ones(numel(names_list),3);
fitT = NaN*ones(numel(names_list),4);
               
% Loop through the files in this directory
for jj = 1:numel(names_list)
    imname = names_list{jj};
    [lc(jj), tc(jj),fitR(jj,:),fitT(jj,:)] = coarse_grain_r_and_t_frame2(directory,imname,user_lastframe,r_scale,t_scale);
end

% Save collective outputs in .mat format
save([directory 'Coarse Test\Forced\coarse_scales.mat'],'names_list','lc','tc','fitR','fitT')