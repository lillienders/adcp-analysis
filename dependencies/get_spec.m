function [fw,pxx] = get_spec(cell,all_vel,subset)
% get_spec: uses pwelch to compute spectrum of subset (or all) ensembles 
% inputs: cell - depth at which you want to calculate spectra, all_vel -
% output matrix from get_ensembles function with sorted ensemble
% measurements, subset - subset of spectra that you want to look through
% outputs: fw - frequency (x-axis), pxx - spectral powers for each ensemble
% in the subset (y-axis)
%%
pxx = zeros(size(all_vel,3),129);
fw = zeros(1,129);
all_vel(isnan(all_vel))=0;
for i=1:length(subset)
    y = all_vel(:,cell,subset(i));
    nfft = length(y); 
    [pxx(i,:),fw] = pwelch(y,60,[],[]);
    hold on
end
end

