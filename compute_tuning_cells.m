function [observed_average_firing,shuffle_firing_threshold,shuffled_average_firing] = compute_tuning_cells(calcium_data,calcium_time,calcium_fps,...
                                                    behav_data,behav_time,behav_fps,...
                                                    nrand,percentile_threshold)

% This function evaluates if a cell is significantly more active during a
% given behavior compares to shuffled distribution of averaged firing

% INPUTS:   calcium_data, neurons x time matrix
%           calcium_time, 1 x time matrix
%           calcium_fps, acquisition rate of endoscope camera
%           behav_data, behaviors x time matrix
%           behav_time, 1 x time matrix
%           behav_fps, acquisition rate of behavior camera
%           nrand, number of permutation for significance level
%           percentile_threshold, percentile of the shuffled distribution to detect significant neuron
%
% OUTPUTS:  observed_average_firing, neurons x behaviors matrix, observed value of average firing for each behavior 
%           shuffle_firing_threshold, neurons x behaviors matrix, threshold value for significance of average firing for each behavior defined using percentile_threshold 
%           shuffled_average_firing, neuronsx behaviors x nrand matrix, shuffled distribution of average activity


%init output
observed_average_firing = NaN(size(event,1),size(behav_data,1));
shuffled_average_firing = NaN(size(event,1),size(behav_data,1),shuffle_nb);


%select commun time bins between calcium recording and behavior recording
[~,ind_cal,ind_loco] = intersect(round(calcium_time*calcium_fps)/calcium_fps,round(behav_time*behav_fps)/behav_fps);
event = calcium_data(:,ind_cal);
behav = behav_data(:,ind_loco);

%matrix of indices for shuffle
mat_shuffle = NaN(shuffle_nb,size(event,2));
for rr=1:1:nrand
    mat_shuffle(rr,:) = randperm(size(event,2));
    shuff_event = event(:,mat_shuffle(rr,:));
    for bb=1:1:size(behav_data,1)
        shuffled_average_firing(:,bb,rr) = sum(shuff_event(:,behav(bb,:)==1),2)/sum(behav(bb,:)==1);
    end
end
                        
for bb=1:1:size(behav_data,1)
    observed_average_firing(:,bb) = sum(event(:,behav(bb,:)==1),2)/sum(behav(bb,:)==1);
end

%calculate threshold values                                 
shuffle_firing_threshold = quantile(shuffled_average_firing,percentile_threshold,3);

end