%{
- Computes & plots logarithmic fits 
- Includes options to average by tide stage or wind conditions
%}
%%
clc
clear all
close all
%% Load along stream velocity
load('/Users/lillienders/Desktop/ADCP/DEC-17/QC/FAST_3_1_QC.mat');
load('/Users/lillienders/Desktop/ADCP/DEC-17/DATA/Ensembles/along_stream.mat')

phase = Data.phase;
wind_idx = Data.wind_idx;
Nz=45; 
z=Data.IBurst_Range;
z=z(1:Nz)';
for k=1:Nz
    var_mean(k,:)=Profile(k).vel_mean;
    var_std(k,:)=Profile(k).vel_std;
    t(k,:)=Profile(k).time_ens_mean;
end

s_signed = var_mean';
s_err = var_std';
NEns=size(var_mean,2);
day0=datenum(2018,0,0,0,0,0);
yd=t-day0+1;
zplot=repmat(z,1,NEns);

%% Remove Slack Values
for i=1:length(s_signed)
    if abs(mean(s_signed(i,:),'omitnan'))<1
       s_signed(i,:)= NaN;
    end
end
%% Choose tide stage and wind flag to process 
tide_stage = 1;
wind_flag = 0;
s_temp = zeros(sum(phase == tide_stage & wind_idx == wind_flag),45);
err_temp = zeros(sum(phase == tide_stage & wind_idx == wind_flag),45);
j=1;
for i=1:length(phase)
    if phase(i) == tide_stage && wind_idx(i) == wind_flag
        s_temp(j,:) = s_signed(i,:);
        err_temp(j,:) = s_err(i,:);
        j=j+1;
    end
end
%% Averaging set up (choose bins)
umin = -6;
dumean = 0.25;
umax = 6;

edges = umin:dumean:umax;
centers = umin+dumean/2:dumean:umax-dumean/2;

[bincounts,edg,ind] = histcounts(s_temp(:,20),edges);
for i = 1:length(bincounts)
    av = find(ind==i);
    u_mean(i,:) = mean(s_temp(av,:),1,'omitnan');
    u_std(i,:) = sqrt(var(s_temp(av,:),1,'omitnan'));
    num_ens(i)=length(av);
end
%% Optional - If you want to save any parameters to look at/plot later
zos = zeros(length(u_mean),1);
u_stars = zeros(length(u_mean),1);
u_means = zeros(length(u_mean),1);
u_1m = zeros(length(u_mean),1);
u_da = zeros(length(u_mean),1);
%% Logarithmic Fit - Linear to get u* & zo 
kappa = 0.41;
fit = zeros(size(u_mean));
rsquare = zeros(size(u_mean));

f1 = figure(1);
set(gcf,'Position',[400 500 500 400])
for i=1:48
    ii=u_mean(i);
    semilogy(u_mean(i,:),z,'.-r','MarkerSize',10);
    hold on;
    errorbar(u_mean(i,:),z,u_std(i,:),'horizontal','r','MarkerSize',10);
    % Linear fit to find prams    
    for k = 5:length(z)
        x = log(z(1:k));
        y = u_mean(i,1:k);
        if isnan(y) < 1
            c = polyfit(x,y,1);
            y_est = polyval(c, x)';
            SStot = sum((y-mean(y)).^2); 
            SSres = sum((y-y_est).^2);
            rsquare(i,k) = 1-SSres/SStot;
        end
    end
    [M,I] = max(rsquare(i,:));
    x = log(z(1:I));
    y = u_mean(i,1:I)';
    if isnan(y) < 1
        c = polyfit(x,y,1);
        y_est = polyval(c, x);
        u_star = c(1)*kappa;
        zo = exp(c(2)/-c(1));
        fit(i,:) = c(1)*log(z/zo);
        zos(i) = zo;
        u_stars(i) = u_star;
        u_means(i) = mean(u_mean(i));
        u_1m(i) = u_mean(i,1);
        u_mean(u_mean==0) = NaN;
        u_da(i) = mean(u_mean(i,:),2,'omitnan');
        hold on
        semilogy(y_est,exp(x),'k','LineWidth',1.5)
    end
end
ylim([0 50])
xlabel('Along Stream Velocity (m/s)')
ylabel('z (m)')
yt = get(gca, 'YTick');
ytkvct = 10.^linspace(1, 10*size(yt,2), 10*size(yt,2));
set(gca, 'YTick', ytkvct);
set(gca, 'YMinorTick','on', 'YMinorGrid','on')
set(gca, 'YTick', [10.^0 10.^1 20.^1 30.^1 40.^1])
set(gca,'FontSize',16,'fontname','times')

%% Logarithmic Fit - Look at Fits
set(gcf,'Position',[400 500 500 400])
for i=1:45
    ii=i;
    he = plot(u_mean(ii,:),z,'.-b','MarkerSize',10);
    hold on;
    errorbar(u_mean(i,:),z,u_std(i,:),'horizontal','b','MarkerSize',10);
    hold on
    plot(fit(ii,1:34),z(1:34),'k','Linewidth',2);
end
ylim([0 35])
set(gca,'FontSize',14,'fontname','times')
set(gcf,'PaperPositionMode','auto')
xlabel('Along Stream Velocity (m/s)')
ylabel('Depth')
title('(c)')
set(gca,'FontSize',15,'fontname','times')
xlim([0 6])
