function [l_c, t_c, pR, pT] = coarse_grain_MultiC(directory,imname,final_frame_user,r_scale,t_scale)
% coarse_grain_r_t but with std calculated per frame compared to overall
% also deleted plotting parts for simplicity
% 2015-06-01 RML - THIS ISN'T REALLY COMPATIBLE WITH MUTLIC IN THAT IT ONLY
% HANDLES ONE CENTER, BUT WILL TAKE THE X Y COORDS FROM MULTIC
% 2015-12-14 RML -- If final_frame from edge detection is too short, output
% NaN for all time variables
% 2016-01-06 RML -- fix a problem with finding dx if the sides are messed
% up and added a different way of finding the edges that doesn't fail if a
% tiny piece of the center has r_scaled > edge_cut
% 2016/04/25 fix error in calculating r (was missing sqrt!)

%%%%% PARAMETERS %%%%%%
r_sizes = [2 3 4 6 8 9]; %in grid points (not um or pixels)
t_sizes = [2 4 10 20 40 60 80 90 100 110 120 130]; %must be even, in frames
edge_cut = 0.75; % only keep r/R > edge_cut
options = optimset('Display','none');

% Set up a save dir if necessary
if ~exist([directory 'Coarse Test'],'file')
    mkdir(directory,'Coarse Test')
end

if ~exist([directory 'Coarse Test\Forced'],'file')
    mkdir([directory 'Coarse Test'],'Forced')
end

[u,v,x,y,a,final_frame,num_pieces] = set_up_edge_piv_MultiC(directory, imname,final_frame_user,edge_cut,r_scale,t_scale);

count = 0;
find_edge = NaN*ones(num_pieces,2);
for kk = 1:num_pieces
    u_slice = u(:,((kk-1)*a/num_pieces + 1):(kk*a/num_pieces),:);
    u_slice = ~isnan(u_slice);
    find_edge(kk,:) = [sum(u_slice(:)),kk];
end
find_edge = sortrows(find_edge,-1);
    
for kk = find_edge(1:2,2)'
%     if sum(u_slice(:))
        if count
            xR = x(:,((kk-1)*a/num_pieces + 1):(kk*a/num_pieces),:);
            yR = y(:,((kk-1)*a/num_pieces + 1):(kk*a/num_pieces),:);
            uR = u(:,((kk-1)*a/num_pieces + 1):(kk*a/num_pieces),:);
            vR = v(:,((kk-1)*a/num_pieces + 1):(kk*a/num_pieces),:);
        else
            xL = x(:,((kk-1)*a/num_pieces + 1):(kk*a/num_pieces),:);
            yL = y(:,((kk-1)*a/num_pieces + 1):(kk*a/num_pieces),:);
            uL = u(:,((kk-1)*a/num_pieces + 1):(kk*a/num_pieces),:);
            vL = v(:,((kk-1)*a/num_pieces + 1):(kk*a/num_pieces),:);
            count = count + 1;
        end
%     end
end

% pcolor(xR(:,:,1)), shading flat, axis equal
% figure;
% pcolor(xL(:,:,1)), shading flat, axis equal
% error('why')

% % Now split into two pieces to avoid weird spatial averaging of disconnected pieces
% xR = x(:,(a/3+1):2/3*a,:);
% yR = y(:,(a/3+1):2/3*a,:);
% uR = u(:,(a/3+1):2/3*a,:);
% vR = v(:,(a/3+1):2/3*a,:);
% 
% xL = x(:,1:a/3,:);
% yL = y(:,1:a/3,:);
% uL = u(:,1:a/3,:);
% vL = v(:,1:a/3,:);

% Calculate original grid size
dx = xR(1,2) - xR(1,1);
if isnan(dx)
    dx = xR(end,2) - xR(end,1);
    if isnan(dx)
        dx = xR(1,end) - xR(1,end-1);
        if isnan(dx)
            error('Check the PIV x coordinates')
        end
    end
end

%%%% Values for non-coarse-grained grid
% Calculate radial velocities
    rL = sqrt(xL.^2 + yL.^2);  %%%%% 2016/04/25 ERROR -- wasn't taking sqrt before...
    rR = sqrt(xR.^2 + yR.^2);
% v_r = (u.*x + v.*y)./r; % Dot product of velocity and r_hat (polar coordinates)
v_thL = (-uL.*yL + vL.*xL)./rL; % Dot product of velocity and theta_hat (polar coordinates)
v_thR = (-uR.*yR + vR.*xR)./rR; % Dot product of velocity and theta_hat (polar coordinates)
v_th = [v_thR v_thL];

% Preallocate storeage
% mean_vth_r = NaN*ones(1,length(r_sizes)+1);
% std_vth_r = NaN*ones(1,length(r_sizes)+1);
std_vth_r_f = NaN*ones(1,length(r_sizes)+1);
std_vth_r_fs = NaN*ones(1,length(r_sizes)+1); %Deviation over frames; weighted fit

% mean_vth_t = NaN*ones(1,length(t_sizes)+1);
% std_vth_t = NaN*ones(1,length(t_sizes)+1);
std_vth_t_f = NaN*ones(1,length(t_sizes)+1);
std_vth_t_fs = NaN*ones(1,length(t_sizes)+1);
% Calculate values for original grid
% mean_vth_r(1) = mean(v_th(~isnan(v_th)));
% std_vth_r(1) = std(v_th(~isnan(v_th)));
% By frame for r
std_temp = NaN*ones(size(v_th,3),1);
for kk = 1:size(v_th,3)
    slice = v_th(:,:,kk);
    std_temp(kk) = std(slice(~isnan(slice)));
end
std_vth_r_f(1) = mean(std_temp);
std_vth_r_fs(1) = std(std_temp);

% mean_vth_t(1) = mean_vth_r(1); % r and t have the same starting point
% std_vth_t(1) = std_vth_r(1);
std_vth_t_f(1) = std_vth_r_f(1);
std_vth_t_fs(1) = std_vth_r_fs(1);

%%%%%%%%%%%%%%%%% Spatial Coarse Graining
% Set up counter for later calculations
count = 2;
% Loop through different grid sizes
for kk = r_sizes
    
    x_new = coarse_grain_grid_r(xL,kk);
    y_new = coarse_grain_grid_r(yL,kk);
    u_new = coarse_grain_grid_r(uL,kk);
    v_new = coarse_grain_grid_r(vL,kk);
    
    r = sqrt(x_new.^2 + y_new.^2); %%%% 2016/04/25 error -- previously not taking sqrt

    v_thL = (-u_new.*y_new + v_new.*x_new)./r; % Dot product of velocity and theta_hat (polar coordinates)
    
    x_new = coarse_grain_grid_r(xR,kk);
    y_new = coarse_grain_grid_r(yR,kk);
    u_new = coarse_grain_grid_r(uR,kk);
    v_new = coarse_grain_grid_r(vR,kk);
    
    r = sqrt(x_new.^2 + y_new.^2); %%%% 2016/04/25 error -- previously not taking sqrt

    v_thR = (-u_new.*y_new + v_new.*x_new)./r; % Dot product of velocity and theta_hat (polar coordinates)
    v_th = [v_thR v_thL]; %Combine for averaging
    
%     mean_vth_r(count) = mean(v_th(~isnan(v_th)));
%     std_vth_r(count) = std(v_th(~isnan(v_th)));
    
    std_temp = NaN*ones(size(v_th,3),1);
    for jj = 1:size(v_th,3)
        slice = v_th(:,:,jj);
        std_temp(jj) = std(slice(~isnan(slice)));
    end
    std_vth_r_f(count) = mean(std_temp);
    std_vth_r_fs(count) = std(std_temp);
    
    count = count+1;    

end

r_sizes = [1 r_sizes];
r_to_plot = (1:0.1:19)*dx;

% plot(dx*r_sizes,std_vth_r_f,'sk','LineWidth',2)

% Start with all frames data
p = polyfit(dx*r_sizes,log(std_vth_r_f),1); % Get initial guesses
x0 = [-1/p(1) real(log(p(2))) min(std_vth_r_f)];

W = 1./std_vth_r_fs;
[pA, resA] = lsqnonlin(@expweightO,x0,[0 -Inf -Inf],[],options,dx*r_sizes,std_vth_r_f);

x0 = pA; x0(1) = 35; %force it down to what has been previously seen
[pB, resB] = lsqnonlin(@expweightO,x0,[0 -Inf -Inf],[],options,dx*r_sizes,std_vth_r_f);

% Keep the better fit
if resA < resB
    pR = pA;
else
    pR = pB;
end
l_c = pR(1);
fvalR = pR(3) + pR(2)*exp(-r_to_plot/pR(1));


figure(1)
plot(r_to_plot,fvalR,'k',dx*r_sizes,std_vth_r_f,'sk','LineWidth',2)
legend(['l_c = ' num2str(l_c,'%0.0f') ' \mum'],'Calculated \sigma(v_\theta)')
xlabel('Grid Size (\mum)','FontSize',16)
ylabel('\sigma(v_\theta)','FontSize',16)
set(gca,'FontSize',16)
xlim([0 200])

savename = [directory 'Coarse Test\Forced\' imname '_r fit each frame'];
saveas(gcf,[savename '.tif'],'tif')
saveas(gcf,[savename '.fig'],'fig')
 

%%%%%%%%%%%%%%%% Time Coarse Graining
if final_frame > max(t_sizes) %If edge detection is enough, keep going -- otherwise everything = NaN;
    % Set up counter for later calculations
    count = 2;
    % Loop through different frame groupings
    for kk = t_sizes

        x_new = coarse_grain_grid_t(xL,kk);
        y_new = coarse_grain_grid_t(yL,kk);
        u_new = coarse_grain_grid_t(uL,kk);
        v_new = coarse_grain_grid_t(vL,kk);

        r = sqrt(x_new.^2 + y_new.^2); %%%% 2016/04/25 error -- previously not taking sqrt

        v_thL = (-u_new.*y_new + v_new.*x_new)./r; % Dot product of velocity and theta_hat (polar coordinates)

        x_new = coarse_grain_grid_t(xR,kk);
        y_new = coarse_grain_grid_t(yR,kk);
        u_new = coarse_grain_grid_t(uR,kk);
        v_new = coarse_grain_grid_t(vR,kk);

        r = sqrt(x_new.^2 + y_new.^2); %%%% 2016/04/25 error -- previously not taking sqrt

        v_thR = (-u_new.*y_new + v_new.*x_new)./r; % Dot product of velocity and theta_hat (polar coordinates)
        v_th = [v_thR v_thL]; %Combine for averaging

    %     mean_vth_t(count) = mean(v_th(~isnan(v_th)));
    %     std_vth_t(count) = std(v_th(~isnan(v_th)));

        std_temp = NaN*ones(size(v_th,3),1);
        for jj = 1:size(v_th,3)
            slice = v_th(:,:,jj);
            std_temp(jj) = std(slice(~isnan(slice)));
        end
        std_vth_t_f(count) = mean(std_temp);
        std_vth_t_fs(count) = std(std_temp);  

        count = count+1;    

    end

    t_sizes = [1 t_sizes];
    t_to_plot = (1:final_frame)*t_scale/60;
    % 
    % disp(std_vth_t_f)

    % Start with all frames data
    pt = polyfit(t_scale/60*t_sizes,log(std_vth_t_f),1);
    x0 = [-1/pt(1) real(log(pt(2))) min(std_vth_t_f)];

    W = 1./std_vth_t_fs;

    [pT,~] = lsqnonlin(@expweightO3,[x0 x0(2)/2],[0 -Inf -Inf -Inf],[],options,t_scale/60*t_sizes,std_vth_t_f);
    fvalT = pT(3) + pT(2)*exp(-t_to_plot/pT(1))+pT(4)*exp(-t_to_plot/(3/60));
    t_c = pT(1);
    % display(res1)

    % [pft,res2] = lsqnonlin(@expweightO,x0,[0 -Inf -Inf],[],options,t_scale/60*t_sizes,std_vth_t_f);
    % fvalft = pft(3) + pft(2)*exp(-t_to_plot/pft(1));
    % t_cf = pft(1);
    % display(res2)

    figure(2)
    plot(t_to_plot,fvalT,'k',t_scale/60*t_sizes,std_vth_t_f,'sk','LineWidth',2)
    legend(['t_c = ' num2str(t_c,'%0.0f') ' h'],'Calculated \sigma(v_\theta)')
    xlabel('Frames (h)','FontSize',16)
    ylabel('\sigma(v_\theta)','FontSize',16)
    set(gca,'FontSize',16)
    xlim([0 final_frame*t_scale/60])

    savename = [directory 'Coarse Test\Forced\' imname '_t fit each frame'];
    saveas(gcf,[savename '.tif'],'tif')
    saveas(gcf,[savename '.fig'],'fig')

else
    t_to_plot = NaN;
    pT = NaN;
    t_c = NaN;
    std_vth_t_f = NaN;
    fvalT = NaN;
end


%%%%% Save relevant data
save([directory 'Coarse Test\Forced\' imname '_coarse_data.mat'],'imname','directory',...
    'r_sizes','t_sizes','r_scale','t_scale','dx',...
    'r_to_plot','t_to_plot','pR','pT','t_c','l_c',...
    'std_vth_r_f','std_vth_t_f','fvalR','fvalT')

%%%%%%%%%% Exponential Fitting Functions

    % Fit to a single exponential with offset - 3 fit parameters
    function diff = expweightO(p,X,Y)
        % Y = p(3) + p(2)*exp(-X/p(1))
        fitY = p(2)*exp(-X/p(1))+p(3);
        raw_diff = fitY - Y;
        diff = raw_diff.*sqrt(W);
    end

    % Fit to a double exponential with offset and 3 min frame rate - 4 fit parameters
    function diff = expweightO3(p,X,Y)
        % Y = p(3) + p(2)*exp(-X/p(1))
        fitY = p(2)*exp(-X/p(1))+p(3)+p(4)*exp(-X/(t_scale/60));
        raw_diff = fitY - Y;
        diff = raw_diff.*sqrt(W);
    end

%     % Fit to a single exponential with offset, unweighted - 3 fit parameters
%     function diff = expO(p,X,Y)
%         % Y = p(3) + p(2)*exp(-X/p(1))
%         fitY = p(2)*exp(-X/p(1))+p(3);
%         diff = fitY - Y;
%     end


end
