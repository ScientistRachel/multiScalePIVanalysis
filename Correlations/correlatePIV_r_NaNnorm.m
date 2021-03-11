% Performs spatial correlations on data in a PIV grid format
%
% Created by Rachel Lee, 2015/07/22
% 2015/07/23 Added x and y for more accurate answers when combining
% matrices (such as left and right edge speeds)
% 2015/08/27 Save more information for better plotting later
% 2016/06/22 May have previously been calculating distances in wrong
% units?? dists is now calcuated from x and y assumed to be in correct
% units -- r_scale is no longer an input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function correlatePIV_r_NaNnorm(piv_values,x,y,t_scale,sub_mean,lastframe,max_dist,savedir,title_text)
function correlatePIV_r_NaNnorm(piv_values,x,y,t_scale,sub_mean,lastframe,max_dist,savedir,title_text)
% function correlatePIV_r_NaNnorm(piv_values,x,y,t_scale,r_scale,sub_mean,lastframe,max_dist,savedir,title_text)

% Hard coded parameters
ylim_to_use = [0 1];
r_bin_size = 15;
time_bin = floor(60/t_scale); %2016/01/31 make it 1 hour instead of fixed frame number

% Default Parameters
if ~exist('t_scale','var') || isempty(t_scale)
    t_scale = 1; %min per frame
end
% if ~exist('r_scale','var') || isempty(r_scale)
%     r_scale = 1; %um/pixel
% end
if ~exist('sub_mean','var') || isempty(sub_mean)
    sub_mean = 0; % subtract the mean before correlating? 0 = no
end
if ~exist('lastframe','var') || isempty(lastframe)
    lastframe = inf; % use all frames as default
end
% if ~exist('PIVtopixel','var') || isempty(PIVtopixel)
%     PIVtopixel = 1; %pixels/PIV grid size
% end
if ~exist('max_dist','var') || isempty(max_dist)
    max_dist = inf; % use full image size as default
end

[a, b, c] = size(piv_values);
c = min(c,lastframe);
% If requested, subtract the mean
if sub_mean
    for kk = 1:c
        slice = piv_values(:,:,kk);
        slice = slice - mean(slice(~isnan(slice)));
        piv_values(:,:,kk) = slice;
    end
end

% Calculate all possible distances
dists = NaN*ones(a,b,a,b);
for i = 1:a
    for j = 1:b
        for k = i:a
            for m = j:b
                dists(i,j,k,m) = sqrt( (x(i,j) - x(k,m))^2 + (y(i,j) - y(k,m)).^2);
%                 dists(i,j,k,m) = sqrt((i-k)^2+(j-m)^2);
            end
        end
    end
end

% Convert distances to useful units
% dists = dists*r_scale*PIVtopixel; %now in um
% dists = dists*r_scale; %now in um
max_dist = min(nanmax(dists(:)), max_dist);

% Bins for averaging
r_bins = 0:r_bin_size:max_dist;
bin_mids = r_bins(1:end-1) + r_bin_size/2;
bin_mids = [0 bin_mids];

% Preallocate storage
Cr = NaN*ones(length(r_bins),c);
N = Cr;
for jj = 1:length(r_bins)   

    if jj == 1
        s = find(dists == r_bins(jj));
    else
        s = find((dists > r_bins(jj-1)) & (dists < r_bins(jj)));
    end
    
    [r, s, t, u] = ind2sub(size(dists),s);
    subs1 = sub2ind(size(piv_values),r,s);
    subs2 = sub2ind(size(piv_values),t,u);
    
    for kk = 1:c
        slice = piv_values(:,:,kk);
        slice1 = slice(subs1);
        slice2 = slice(subs2);
        
        % Especially at large distances, some of the indices might be in
        % the non-cell region  -- the NaN mask makes sure the terms have
        % the same number of ~isNaN values for summing
        NaNmask = slice1.*slice2./(slice1.*slice2);
        slice1 = NaNmask.*slice1;
        slice2 = NaNmask.*slice2;
        
        corr_temp = slice1.*slice2;        
        term1 = nansum(nansum(slice1.^2));
        term2 = nansum(nansum(slice2.^2));
        Cr(jj,kk) = nansum(nansum(corr_temp))./(sqrt(term1.*term2));
        N(jj,kk) = sum(sum(~isnan(corr_temp)));
    end   
    
end

% Average over time
Cr_mean = mean(Cr,2);

% Plot time average
figure(1)
plot(bin_mids,Cr_mean,'s','LineWidth',2)
ylim(ylim_to_use)
xlabel('\Deltar (\mum)','FontSize',16)
ylabel('C(\Deltar)','FontSize',16)
set(gca,'FontSize',16)
xlim([min(bin_mids) max(bin_mids)])
ylim(ylim_to_use)
% Optionally title
if exist('title_text','var') && ~isempty(title_text)
    title(title_text,'FontSize',16)
end
% Save if there is a file name given
if exist('savedir','var') && ~isempty(savedir)
    savename = [savedir ' Spatial Correlation - NaNnorm '];
    saveas(gcf,[savename '.fig'],'fig')
    saveas(gcf,[savename '.png'],'png')
end

% Group by time
time_pieces = floor(c/time_bin);
time_step = time_bin*t_scale/60;
time_vec = 1:time_step:time_step*time_pieces;
color_label = cell(1,length(time_vec));
for kk = 1:length(time_vec)
    color_label{kk} = num2str(time_vec(kk));
end
x_tick_places = 0:(1/time_pieces):((time_pieces-1)/time_pieces);
x_tick_places = x_tick_places + mean(diff(x_tick_places))/2;

% Take average over time_bin sizes
n = NaN*ones(size(Cr,1),time_pieces);
for kk = 1:time_pieces
    
    frames = (time_bin*(kk-1)+1):(time_bin*kk);    
    slice = Cr(:,frames);
    n(:,kk) = mean(slice,2);
    
end

figure(2)
cmap = parula(size(n,2));
plot(bin_mids,n(:,1),'LineWidth',2,'Color',cmap(1,:))
hold on
for kk = 2:size(n,2)
    plot(bin_mids,n(:,kk),'LineWidth',2,'Color',cmap(kk,:))
end
hold off
xlabel('\Deltar (\mum)','FontSize',16)
ylabel('C(\Deltar)','FontSize',16)
set(gca,'FontSize',16)
xlim([min(bin_mids) max(bin_mids)])
ylim(ylim_to_use)
colormap(cmap)
% caxis manual
% caxis([1 time_pieces])
h = colorbar;
set(h,'Ticks',x_tick_places,'TickLabels',color_label,'TickLength',0)
set(get(h,'Label'),'String','Hour','FontSize',16)

% Optionally title
if exist('title_text','var') && ~isempty(title_text)
    title(title_text,'FontSize',16)
end


% Save if there is a file name given
if exist('savedir','var') && ~isempty(savedir)
    savename = [savedir ' Spatial Correlation over Time - NaNnorm'];
    saveas(gcf,[savename '.fig'],'fig')
    saveas(gcf,[savename '.png'],'png')
end

save([savedir 'NaNnorm.mat'],'r_bin_size','time_bin','Cr','N','Cr_mean',...
    'time_vec','n','max_dist','bin_mids')

end