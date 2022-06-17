function [all_vel] = get_ensembles(fs,timeall,vel)
% get_ensembles.m: Sorts raw pings into 5 minute ensemble bins, which can then be used to compute spectra
%inputs: fs - sampling frequency of adcp in hz, timeall - concatenated time
%(load in spec file), vel - concatenated velocity you want ensembles of

%outputs: all_vel - structure containing sorted velocities with dimensions
%of [fs*300 x no. depth cells x no. ensembles]
%%
MinEns=5;
dt=1/fs;
timeall_sec=(timeall-timeall(1))*3600*24;

Ens_sec=find(diff(timeall_sec)>100); % or Ens==find(diff(timeall)>1/24/60) %1 minute in days
Ens=find(diff(timeall)>1/24/60); % or 0.001 also works, I used plus 1 min
NEns=length(Ens); % How many ensembles

sizeens=diff(Ens); %Size of each ensembles
index_1=zeros(NEns,1);
index_end=zeros(NEns,1);
index_1(1,1)=1;

for j = 1:length(Ens)-1
    index_end(j,1)=Ens(j);
    index_1(j+1,1)=Ens(j)+1;
    size_ens(j,1)=index_end(j,1)-index_1(j,1)+1; 
end

lastEnssize=length(timeall)-index_end(end);
index_end(NEns,1)=length(timeall);
size_ens(NEns,1)=index_end(NEns,1)-index_1(NEns,1)+1;

k=1;
for j = 1:length(Ens)
    if size_ens(j)==fs*5*60
        time_ens(:,k)=timeall(index_1(j,1):index_end(j,1));
        k=k+1;
    end
end
GEns=size(time_ens,2); %How many good ensembles are there?
for j=1:GEns
    time_ens_sec(:,j)=(time_ens(:,j))-time_ens(1,j)*3600*24;
end
k=1;
for j=1:NEns   
    if size_ens(j)==fs*5*60 %Save data
        vel_ens=vel(index_1(j,1):index_end(j,1),:); %Here I grab all z values to make it faster
        all_vel(:,:,k)=vel_ens;
        k=k+1;
    end
end
end

