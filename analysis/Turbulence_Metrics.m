%{
- Calculates turbulence metrics (Beam 5 variance, TKE, TI) and makes colour plots 
- Set up to plot in sigma coordinates (converted using to_sigma function)
Last Edit: June 15 2022
Set Up For: March 2018 Sig 500 Deployment (FORCE)
%}
%%
close all
clear all
clc
%% Paths
path_to_QC = '~/Desktop/ADCP/MAR-18/QC/MAR18_1_QC.mat'; % Path to QC File
path_to_fs = '~/Desktop/ADCP/MAR-18/DATA/Attitude-Params/FreeSurfaceEnsembles.mat'; % Path to free surface measurement 
%% Noise Var
theta=25*pi/180;
load('~/Desktop/ADCP/MAR-18/Raw-Data/MAR18_1.mat'); % Load configuration file - has beam precision 

nb = Config.burst_beamPrecision/100; %Along Beam precision in m/s
nb_var=nb^2;
nv = Config.burst_verticalPrecision/100; %Vertical velocity precision in m/s
nv_var=nv^2;
nh = Config.burst_horizontalPrecision/100; %Horizontal velocity precision in m/s
nh_var=nh^2;
%% Load beam variances and convert to sigma
%  Beam 5
load('~/Desktop/ADCP/MAR-18/DATA/Ensembles/Ens_v5.mat');
Nz = 45;
for k=1:Nz
    v5_mean(k,:)=Profile(k).vel_mean;
    v5_var(k,:)=Profile(k).vel_var;
    t(k,:)=Profile(k).time_ens_mean;
end
clear Profile

day0=datenum(2018,0,0,0,0,0);
yd=t-day0+1;

[yd_sigma,v5_var_sigma,u_std,sigma_norm] = to_sigma(path_to_fs,path_to_QC,yd,v5_var);

% Beam 4
load('/Users/lillienders/Desktop/ADCP/DEC-17/DATA/Ensembles/Ens_v4.mat');
Nz = 45;
for k=1:Nz
    v4_var(k,:)=Profile(k).vel_var;
    t(k,:)=Profile(k).time_ens_mean;
end
clear Profile
day0=datenum(2018,0,0,0,0,0);
yd=t-day0+1;
[yd_sigma,v4_var_sigma,u_std,sigma_norm] = to_sigma(path_to_fs,path_to_QC,yd,v4_var);

% Beam 3
load('/Users/lillienders/Desktop/ADCP/DEC-17/DATA/Ensembles/Ens_v3.mat');
Nz = 45;
for k=1:Nz
    v3_var(k,:)=Profile(k).vel_var;
    t(k,:)=Profile(k).time_ens_mean;
end
clear Profile
day0=datenum(2018,0,0,0,0,0);
yd=t-day0+1;
[yd_sigma,v3_var_sigma,u_std,sigma_norm] = to_sigma(path_to_fs,path_to_QC,yd,v3_var);

% Beam 2
load('/Users/lillienders/Desktop/ADCP/DEC-17/DATA/Ensembles/Ens_v2.mat');
Nz = 45;
for k=1:Nz
    v2_var(k,:)=Profile(k).vel_var;
    t(k,:)=Profile(k).time_ens_mean;
end
clear Profile
day0=datenum(2018,0,0,0,0,0);
yd=t-day0+1;
[yd_sigma,v2_var_sigma,u_std,sigma_norm] = to_sigma(path_to_fs,path_to_QC,yd,v2_var);

% Beam 1
load('/Users/lillienders/Desktop/ADCP/DEC-17/DATA/Ensembles/Ens_v1.mat');
Nz = 45;
for k=1:Nz
    v1_var(k,:)=Profile(k).vel_var;
    t(k,:)=Profile(k).time_ens_mean;
end
clear Profile
day0=datenum(2018,0,0,0,0,0);
yd=t-day0+1;
[yd_sigma,v1_var_sigma,u_std,sigma_norm] = to_sigma(path_to_fs,path_to_QC,yd,v1_var);

%% Load Along-Stream Velocity (in case you want to do TI calcs)
load('/Users/lillienders/Desktop/ADCP/DEC-17/DATA/Ensembles/along_stream.mat');
Nz = 45;
for k=1:Nz
    u_mean(k,:)=Profile(k).vel_mean;
    u_var(k,:)=Profile(k).vel_var;
    t(k,:)=Profile(k).time_ens_mean;
end
clear Profile

day0=datenum(2018,0,0,0,0,0);
yd=t-day0+1;
% Convert to sigma coords
[yd_sigma,u_sigma,u_std,sigma_norm] = to_sigma(path_to_fs,path_to_QC,yd,u_mean);

% Remove Slack
for i = 1:length(u_sigma)
    if abs(u_sigma(20,i)) < 2
        u_sigma(20,i)=NaN;
    end
end

%% TKE CALCULATION
% Remove Noise Var
v1_var=v1_var_sigma-nb_var;
v2_var=v2_var_sigma-nb_var;
v3_var=v3_var_sigma-nb_var;
v4_var=v4_var_sigma-nb_var;
v5_var=v5_var_sigma-nb_var;

% Compute u/v/w variances from beam measurements using variance method
u_var=(v2_var+v4_var-2.*(v5_var)*cos(theta)^2)./(2*sin(theta).^2);
v_var=(v1_var+v3_var-2.*v5_var*cos(theta)^2)./(2*sin(theta).^2);
w_var=v5_var;

% Calc TKE
tke=1/2.*(u_var+v_var+w_var);

% Calc Horizontal Variance
hor_var = 1/2.*(u_var+v_var);

%% TI CALCULATION
TI_vertical=sqrt((v5_var_sigma))./abs(v5_mean_sigma); % Vertical
TI_along =sqrt((u_var_sigma))./abs(u_sigma); % Along Stream
%% PLOTS: Beam 5 Variance 
%  With Surface
f1 = figure(1);
set(f1,'Position',[230   459   956   300])
pcolor(yd_sigma,-sigma_norm,v5_var)
shading flat
colormap(flipud(brewermap([],'Spectral')))
cb = colorbar;
caxis([0 0.05])
ylim([-0.95 0.05])
yline(-0.15,'w','LineWidth',4)
yline(0,'k','LineWidth',4)
ylabel(['Sigma Coordinate'])
xlabel('Year day 2018')
ylabel(cb,'Vertical Variance (m^{2}/s^{2})')
cb.Label.FontName = 'times';
cb.Label.FontSize = 16;
xlim([5 10])
set(gca,'FontSize',16,'fontname','times')
set(gcf,'PaperPositionMode','auto')

%  Without Surface
f2 = figure(2);
set(f2,'Position',[230   459   956   300])
v5_var_ns = v5_var;
v5_var_ns(1:6,:)=NaN;
pcolor(yd_sigma,-sigma_norm,v5_var_ns)
shading flat
colormap(flipud(brewermap([],'Spectral')))
cb = colorbar;
caxis([0 0.05])
ylim([-0.95 -0.05])
yline(-0.15,'k','LineWidth',4)
ylabel(['Sigma Coordinate'])
xlabel('Year day 2018')
ylabel(cb,'Vertical Variance (m^{2}/s^{2})')
cb.Label.FontName = 'times';
cb.Label.FontSize = 16;
set(gca,'FontSize',16,'fontname','times')
set(gcf,'PaperPositionMode','auto')

%% PLOTS: TKE
%  With Surface
f3 = figure(3);
set(f3,'Position',[230   459   956   300])
pcolor(yd_sigma,-sigma_norm,tke)
shading flat
colormap(flipud(brewermap([],'Spectral')))
cb = colorbar;
caxis([0 0.06])
ylim([-0.95 -0.05])
yline(-0.15,'w','LineWidth',4)
ylabel(['Sigma Coordinate'])
xlabel('Year day 2018')
ylabel(cb,'TKE (m^{2}/s^{2})')
cb.Label.FontName = 'times';
cb.Label.FontSize = 16;
set(gca,'FontSize',16,'fontname','times')
set(gcf,'PaperPositionMode','auto')

%  Without Surface
f4 = figure(4);
set(f4,'Position',[230   459   956   300])
tke_ns = tke;
tke_ns(1:6,:)=NaN;
pcolor(yd_sigma,-sigma_norm,tke_ns)
shading flat
colormap(flipud(brewermap([],'Spectral')))
cb = colorbar;
caxis([0 0.1])
ylim([-0.95 -0.05])
xlim([5 12])
yline(-0.15,'k','LineWidth',4)
ylabel(['Sigma Coordinate'])
xlabel('Year day 2018')
ylabel(cb,'TKE (m^{2}/s^{2})')
cb.Label.FontName = 'times';
cb.Label.FontSize = 16;
set(gca,'FontSize',16,'fontname','times')
set(gcf,'PaperPositionMode','auto')

%% PLOTS: Horizontal Variance
%  With Surface
f5 = figure(5);
set(f5,'Position',[230   459   956   300])
pcolor(yd_sigma,-sigma_norm,hor_var)
shading flat
colormap(flipud(brewermap([],'Spectral')))
cb = colorbar;
caxis([0 0.06])
ylim([-0.95 -0.05])
yline(-0.15,'w','LineWidth',4)
ylabel(['Sigma Coordinate'])
xlabel('Year day 2018')
ylabel(cb,'TKE (m^{2}/s^{2})')
cb.Label.FontName = 'times';
cb.Label.FontSize = 16;
set(gca,'FontSize',16,'fontname','times')
set(gcf,'PaperPositionMode','auto')

%  Without Surface
f6 = figure(6);
hor_var_ns = hor_var;
hor_var_ns(1:6,:)=NaN;
set(f6,'Position',[230   459   956   300])
tke_ns = tke;
tke_ns(1:6,:)=NaN;
pcolor(yd_sigma,-sigma_norm,hor_var_ns)
shading flat
colormap(flipud(brewermap([],'Spectral')))
cb = colorbar;
caxis([0 0.1])
ylim([-0.95 -0.05])
xlim([5 10])
yline(-0.15,'k','LineWidth',4)
ylabel(['Sigma Coordinate'])
xlabel('Year day 2018')
ylabel(cb,'Horizontal Variance (m^{2}/s^{2})')
cb.Label.FontName = 'times';
cb.Label.FontSize = 16;
set(gca,'FontSize',16,'fontname','times')
set(gcf,'PaperPositionMode','auto')

%% BEAM 5 AMPLITUDE
%  Load beam 5 amplitude and convert to sigma layers
load('/Users/lillienders/Desktop/ADCP/DEC-17/DATA/Ensembles/Ens_v5Amp.mat');
Nz = 45;
for k=1:Nz
    v5_amp(k,:)=Profile(k).vel_var;
    t(k,:)=Profile(k).time_ens_mean;
end
clear Profile
day0=datenum(2018,0,0,0,0,0);
yd=t-day0+1;
[yd_sigma,v5_amp_sig,u_std,sigma_norm] = to_sigma(path_to_fs,path_to_QC,yd,v5_amp);
%% PLOTS: Beam 5 Amplitude 
%  With Surface
f7 = figure(7);
set(f7,'Position',[230   459   956   300])
pcolor(yd_sigma,-sigma_norm,abs(v5_amp_sig))
shading flat
colormap(flipud(brewermap([],'Spectral')))
cb = colorbar;
caxis([0 100])
ylim([-0.95 0.05])
yline(-0.15,'w','LineWidth',4)
yline(0,'k','LineWidth',4)
ylabel(['Sigma Coordinate'])
xlabel('Year day 2018')
ylabel(cb,'Amplitude (dB)')
cb.Label.FontName = 'times';
cb.Label.FontSize = 16;
xlim([5 48])
set(gca,'FontSize',16,'fontname','times')
set(gcf,'PaperPositionMode','auto')

%  Without Surface
f8 = figure(8);
set(f8,'Position',[230   459   956   300])
pcolor(yd_sigma,-sigma_norm,TI_vertical)
shading flat
colormap(flipud(brewermap([],'Spectral')))
cb = colorbar;
caxis([0 0.1])
ylim([-0.95 -0.05])
xlim([5 45])
yline(0,'k','Linewidth',4)
yline(-0.15,'w','LineWidth',4)
ylabel(['Sigma Coordinate'])
xlabel('Year day 2018')
ylabel(cb,'Turbulence Intensity')
cb.Label.FontName = 'times';
cb.Label.FontSize = 16;
set(gca,'FontSize',16,'fontname','times')
set(gcf,'PaperPositionMode','auto')