% Make images of FTLE values over time
clc
clear all
close all

% directory = 'F:\2014-02-06-Ecad-shRNA\2014-02-06-Ecad-shRNA0047.zvhi_Files\2014-02-06-Ecad-shRNA0047_2\';
% directories = {'F:\2014-02-06-Ecad-shRNA\2014-02-06-Ecad-shRNA0047.zvhi_Files\2014-02-06-Ecad-shRNA0047_2\'...
%     'F:\2014-02-08-Ecad-shRNA\2014-02-08-Ecad-shRNA0048.zvhi_Files\2014-02-08-Ecad-shRNA0048_2\',...
%     'F:\2014-02-16-Ecad-shRNA\2014-02-16-Ecad-shRNA0049.zvhi_Files\2014-02-16-Ecad-shRNA0049_2\',...
%     'F:\2014-02-17-Ecad-shRNA\2014-02-17-Ecad-shRNA0050.zvhi_Files\2014-02-17-Ecad-shRNA0050_2\',...
%     'F:\140220 M1M4 shECad 1%HS\140220 M1M4 shE-Cad-0002\',...
%     'F:\140225 M1M4 shECad 1%HS\140225 M1M4 shECad -0005\'};
% directory = 'H:\140507 edge length\140507 M4 inh LPA-0006\';
directory = 'J:\2014-11-19-EcadshM1\2014-11-19-EcadshM1-0031\';
time_scale = 3;

load([directory 'FTLEall.mat'])
load('ftle_cmap.mat')

if ~exist([directory 'FTLE movie plots\'],'file')
    mkdir(directory, 'FTLE movie plots')
end

figure;
pause(2) %user's chance to make the figure whatever size they want

for file_use = 19:20%1:numel(pivfiles)

%     file_use = 19;
    imname = pivfiles{file_use}(1:end-13);
    if ~exist([directory 'FTLE movie plots\' imname],'file')
        mkdir([directory 'FTLE movie plots\'],imname)
    end

    ftle = ftle_all{file_use}*60/time_scale;


    for kk = 1:size(ftle,3)

        slice = ftle(:,:,kk);
        slice(isinf(slice)) = NaN;

        slice(1:2,:,:) = [];
        slice(end-1:end,:,:) = [];
        slice(:,1:2,:) = [];
        slice(:,end-1:end,:) = [];

        pcolor(slice), shading flat
        set(gcf,'renderer','zbuffer')
        set(gca,'DataAspectRatio',[1 1 1],'ydir','reverse')
        axis off
        caxis manual
        caxis([-2 .5])
        h = colorbar;
        set(h,'FontSize',16)
        colormap(ftle_cmap)

        imw = getframe;
        imwrite(imw.cdata(:,:,:),[directory 'FTLE movie plots\' imname '\FTLE_' num2str(kk) '.tif'],'tif', 'Compression','none');

    end
    
end