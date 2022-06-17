%{
- Computes & plots spectra in frequency domain
- Set up for vertical velocities, can also be used for TKE if given input
file
%}
%%
clear all
close all
clc
%% Load QC Data
load('/Users/lillienders/Desktop/ADCP/DEC-17/QC/FAST_3_1_QC.mat');
Nz=45; 
z=Data.IBurst_Range;
z=z(1:Nz)';
phase = Data.phase; 
wind_idx = Data.wind_idx;
%% Set Noise Params
fs = 4; % Sampling frequency in Hz
noise = 0.1421; % Beam precision in m/s
noisevar=noise^2;
%% Load Vertical Velocity and Time (CONCATENATED) 
Beam5=0; % Process beam 5?
load('/Users/lillienders/Desktop/ADCP/DEC-17/DATA/Concatenated/v5.mat') 
load('/Users/lillienders/Desktop/ADCP/DEC-17/DATA/Concatenated/time_all.mat')
vel=v5;
if Beam5==1
    time=timeallV5;
else
    time=timeall;
end
%% Get Ensembles Measurements (Raw Pings, Sorted into 5 Min Bins)
all_vel = get_ensembles(fs,timeall,vel);
%% Calculate Spectra for Subset of Data at Cell (from bottom)
cell = 10;
subset = 1:size(all_vel,3);
[fw,pxx] = get_spec(cell,all_vel,subset);
%% Average by Tide Phase 
mean_ebb_dec = mean(pxx(phase==0,:),1,'omitnan');
mean_flood_acc = mean(pxx(phase==1,:),1,'omitnan');
mean_flood_peak = mean(pxx(phase==2,:),1,'omitnan');
mean_flood_dec = mean(pxx(phase==3,:),1,'omitnan');
mean_ebb_acc = mean(pxx(phase==4,:),1,'omitnan');
mean_ebb_peak = mean(pxx(phase==5,:),1,'omitnan');

all_phases = [mean_ebb_dec(:,:); mean_flood_acc(:,:); mean_flood_peak(:,:); ...
    mean_flood_dec(:,:); mean_ebb_acc(:,:); mean_ebb_peak(:,:)];
%% PLOTS: Plot all spectra, coloured by tide stage 
f1 = figure(1);
c = ['#4098b7';'#e75948';'#a50844';'#fcb466';'#85cfa5';'#505faa'];
for i=1:length(subset)
    if ~isnan(phase(subset(i)))
       p = loglog(fw,pxx(i,:),'linewidth',1);
       p.Color = c(phase(subset(i))+1,:); 
       hold on
    end
end

xlabel('Frequency (Hz)')
ylabel('S (m^2s^-^2Hz^-^1)')
title('All Spectra')
ylim([0.5e-3 1])
set(gca,'FontSize',14,'fontname','times')
set(gcf,'PaperPositionMode','auto')
%% PLOTS: Spectra averaged by tide phase
f2 = figure(2);
c = ['#4098b7';'#e75948';'#a50844';'#fcb466';'#85cfa5';'#505faa'];
for i=1:size(all_phases,1)
    p = loglog(fw,all_phases(i,:),'linewidth',1); 
    p.Color = c(i,:);
    hold on
end

plot([.2 0.8],1e-2*[.2 0.8].^(-5/3),'--k')
text(0.35,1e-1,'~f ^{-5/3}','FontSize',14,'FontName','times')

xlabel('Frequency (Hz)')
ylabel('S (m^2s^-^2Hz^-^1)')
title('Phase Averaged Spectra')
ylim([0.5e-3 1])
legend('Decelerating Ebb','Accelerating Flood','Peak Flood',...
    'Decelerating Flood','Accelerating Ebb','Peak Ebb')
set(gca,'FontSize',14,'fontname','times')
set(gcf,'PaperPositionMode','auto')