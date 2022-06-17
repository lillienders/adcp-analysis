%{
- Adds a new metric to Data structure in QC file, which has an index 
- Do for true u,v,w, and along/across stream velocities (one at a time)
- Get data for each bin, and construct 5 minute ensembles

Last Edit: June 15 2022
Set Up For: March 2018 Sig 500 Deployment (FORCE)
%}
%%
clc
clear all
close all
%% Load Files
load('/Users/lillienders/Desktop/ADCP/DEC-17/QC/FAST_3_1_QC.mat');
load('/Users/lillienders/Desktop/ADCP/DEC-17/DATA/Ensembles/along_stream.mat')
path_to_QC = '/Users/lillienders/Desktop/ADCP/DEC-17/QC/FAST_3_1_QC.mat'; 
%% Organize Variables
Nz=45; 
z=Data.IBurst_Range;
z=z(1:Nz)';
for k=1:Nz
    var_mean(k,:)=Profile(k).vel_mean;
    var_std(k,:)=Profile(k).vel_std;
    t(k,:)=Profile(k).time_ens_mean;
end
%%
vel_mean = var_mean';
vel_err = var_std';

NEns=size(var_mean,2);
day0=datenum(2018,0,0,0,0,0);
yd=t-day0+1;
%% Find maximum (flood) and minimum (ebb) values for each tidal cycle
[pk_fl,locs_fl] = findpeaks(vel_mean(:,20),'MinPeakProminence',4); % Get MAX flood values and index
[pk_eb,locs_eb] = findpeaks(-vel_mean(:,20),'MinPeakProminence',4); % Get MAX ebb values and index

%% Find zero down and zero up crossings (accelerating/decelerating demarcations)
zero_down = [];
zero_up = [];
for k=1:length(vel_mean)-1
    if (sign(vel_mean(k,20))==1 & sign(vel_mean(k+1,20))==-1)
        zero_down = [zero_down,k];
    end
    if (sign(vel_mean(k,20))==-1 & sign(vel_mean(k+1,20))==1)
        zero_up = [zero_up,k];
    end
end
zero_down = zero_down';
zero_up = zero_up';

%% Assign phase values (0-5) to each ensemble
phase = NaN(1,size(var_mean,2));
for i=1:length(locs_fl)-1 % Looping through each tidal cycle
    phase(zero_up(i+1)+1:locs_fl(i)) = 1 ; % Flood Accelerating
    phase(locs_fl(i)+1:zero_down(i+1)) = 3; % Flood Decelerating
    phase(zero_down(i+1)+1:locs_eb(i+1)-1) = 4; % Ebb Accelerating
    phase(locs_eb(i+1):zero_up(i+2)) = 0; % Ebb Decelerating Redundant
end

%% Find Peak Flood and Peak Ebb
for k=1:length(locs_fl)
    std_peak = vel_err(locs_fl(k),20);
    i = 0;
    while abs(vel_mean(locs_fl(k)+i,20)-vel_mean(locs_fl(k),20)) < std_peak
        phase(locs_fl(k)+i) = 2;
        i = i+1;
    end
    i = 0;
    while abs(vel_mean(locs_fl(k)-i,20)-vel_mean(locs_fl(k),20)) < std_peak
        phase(locs_fl(k)-i) = 2;
        i = i+1;
    end
end
for k=1:length(locs_eb)
    std_peak = vel_err(locs_eb(k),20);
    i = 0;
    while abs(vel_mean(locs_eb(k)+i,20)-vel_mean(locs_eb(k),20)) < std_peak
        phase(locs_eb(k)+i) = 5;
        i = i+1;
    end
    i = 0;
    while abs(vel_mean(locs_eb(k)-i,20)-vel_mean(locs_eb(k),20)) < std_peak
        phase(locs_eb(k)-i) = 5;
        i = i+1;
    end
end

%% Plot to check things out
scatter(yd(1,:),vel_mean(:,20),50,phase,'filled');
colormap(flipud(brewermap(6,'Spectral')))
colorbar('Ticks',0:5);
xlim([16 17])
ylim([-4 4])
yline(0)
ylabel('Along Stream Velocity (m/s)')
xlabel('Year Day 2018')
set(gca,'FontSize',15,'fontname','times')
%%
Data.phase = phase;
save('/Users/lillienders/Desktop/ADCP/MAR-18/QC/MAR18_1_QC.mat','Data');