function [BI,sigBI,randBI] = compute_behavior_information(calcium_data,calcium_time,calcium_fps,...
                                                    behav_data,behav_time,behav_fps,...
                                                    nrand)

% This function calculates the 'behavior information' and its significance
% level for each neuron using mutual information

% INPUTS:   calcium_data, neurons x time matrix
%           calcium_time, 1 x time matrix
%           calcium_fps, acquisition rate of endoscope camera
%           behav_data, behaviors x time matrix
%           behav_time, 1 x time matrix
%           behav_fps, acquisition rate of behavior camera
%           nrand, number of permutation for significance level
%
% OUTPUTS:  BI, neurons x 1 matrix, observed value of behavior information (ie mutual information between fiuring and behavior time series)
%           sigBI, neurons x 1 matrix, significance level of BI
%           randBI, neurons x nrand matrix, shuffled distribution of BI values


%init output
MItp = zeros(size(calcium_data,1),size(behav_data,1));
MIrand = zeros(size(calcium_data,1),size(behav_data,1),nrand);

%select commun time bins between calcium recording and behavior recording
[~,ind_cal,ind_loco] = intersect(round(calcium_time*calcium_fps)/calcium_fps,round(behav_time*behav_fps)/behav_fps);
event = calcium_data(:,ind_cal);
behav = behav_data(:,ind_loco);

for bb=1:1:size(behav_data,1)
	rateQ = nanmean(event(:,behav(bb,:)==1),2)./nanmean(event(:,:),2);
	%rateQ(rateQ==0) = eps;
	MItp(:,bb) = sum(behav(bb,:))/size(behav,2)*rateQ.*log2(rateQ);
end
MItp(~isfinite(MItp)) = 0;
MI = sum(MItp,2);

%permutationd to compute significance level
nrand = 1000;
for rr=1:1:nrand
    irand = randperm(size(behav,2));
    rbehav = behav(:,irand);
    for bb=1:1:length(behaviors_list)
        rateQ = nanmean(event(:,rbehav(bb,:)==1),2)./nanmean(event(:,:),2);
        %rateQ(rateQ==0) = eps;
        MIrand(:,bb,rr) = sum(rbehav(bb,:))/size(rbehav,2)*rateQ.*log2(rateQ);
    end
end
MIrand(~isfinite(MIrand)) = 0;
MIrand2 = squeeze(nansum(MIrand,2));
SigTh =   (MI-mean(MIrand2,2))./std(MIrand2,[],2);

BI = MI;
sigBI = SigTh;
randBI = MIrand;