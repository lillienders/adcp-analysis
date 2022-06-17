%{
- Computes & plots wake fits for along-stream velocities
- Includes options to average by tide stage or wind conditions
%}
%%
clear all
close all
clc
%% Load along stream velocity and convert to sigma coordinates
load('/Users/lillienders/Desktop/ADCP/DEC-17/DATA/Ensembles/along_stream.mat')
path_to_QC = '/Users/lillienders/Desktop/ADCP/DEC-17/QC/FAST_3_1_QC.mat';
path_to_fs = '/Users/lillienders/Desktop/ADCP/DEC-17/DATA/Attitude-Params/FreeSurfaceEnsembles.mat';

Nz = 45;
for k=1:Nz
    vel_mean(k,:)=Profile(k).vel_mean;
    t(k,:)=Profile(k).time_ens_mean;
    vel_std(k,:) = Profile(k).vel_std;
end
day0=datenum(2018,0,0,0,0,0);
yd=t-day0+1;
[yd_sigma,u_sigma,u_sig_std,sigma_norm] = to_sigma(path_to_fs,path_to_QC,yd,vel_mean);
%% Format velocities, remove contaminated data (above sig = -0.15)
u_sigma = u_sigma';
phase = Data.phase;
wind_idx = Data.wind_idx;
sigma_range = [7:39];
eta = (1-sigma_norm(sigma_range,1))';
opts = optimset('MaxFunEvals',50000, 'MaxIter',10000);
u_sigma = u_sigma(:,6:39);
u_sigma = u_sigma(:,1:34);
u_sigma(:,34)=[];
%% Choose tide stage and wind flag to process 
tide_stage = 1;
wind_flag = 0;
s_temp = zeros(sum(phase == tide_stage & wind_idx == wind_flag),size(u_sigma,2));
j=1;
for i=1:length(phase)
    if phase(i) == tide_stage && wind_idx(i) == wind_flag
        s_temp(j,:) = u_sigma(i,:);
        j=j+1;
    end
end
%%
[bincounts,edg,ind] = histcounts(s_temp(:,20),edges);
for i = 1:length(bincounts)
    av = find(ind==i);
    u_mean(i,:) = mean(s_temp(av,:),1,'omitnan');
    u_std(i,:) = sqrt(var(s_temp(av,:),1,'omitnan'));
    num_ens(i)=length(av);
end

u_mean(any(isnan(u_mean), 2), :) = [];
%%
opts = optimset('MaxFunEvals',500000, 'MaxIter',100000);

ufit=u_sigma*0;
for ii=1:size(u_mean,1)
    u=abs(u_mean(ii,:));

    b0=[min(u)/10, 2, 0.1];
    
    OLS = @(b) sum((wake_function(eta,b)-u).^2); 
    [B(ii,:),error,exitflag(ii),output(ii)] = fminsearch(OLS, b0, opts);

    ufit(ii,:)=wake_function2(eta,B(ii,:));
    nrmse(ii) =sqrt(error)/sqrt(sum(u.^2));
    plot(u(2:33),(eta(2:33)-1), 'Color','#a50844','Linewidth',1.5)
    hold on
    scatter(u(2:33),(eta(2:33)-1),30,'MarkerEdgeColor','#a50844','MarkerFaceColor','#a50844')
	plot(ufit(ii,:),(eta-1), 'k','LineWidth',1.5)
    
    set(gca,'FontSize',15,'fontname','times')
    xlabel('|U_A_l_o_n_g|')
    ylabel('σ') 
    xlim([0 6])
end