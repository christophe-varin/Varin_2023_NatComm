function [occurrence,occurrence_rd_threshold,occurrence_rd_dist] = compute_occurrence_rate(calcium_data,calcium_time,calcium_fps,...
                                                    behav_data,behav_time,behav_fps,...
                                                    nrand,percentile_threshold)

% This function calculates how often neurons are active during all episodes
% for each beahvior class

% INPUTS:   calcium_data, neurons x time matrix
%           calcium_time, 1 x time matrix
%           calcium_fps, acquisition rate of endoscope camera
%           behav_data, behaviors x time matrix
%           behav_time, 1 x time matrix
%           behav_fps, acquisition rate of behavior camera
%           nrand, number of permutation for significance level
%           percentile_threshold, percentile of the shuffled distribution to detect significant neuron
%
% OUTPUTS:  occurrence, neurons x behaviors matrix, observed value of activation occurrence as a fraction of the number of episodes
%           occurrence_rd_threshold, neurons x behaviors matrix, threshold value of activation occurrence defined by quantile percentile_threshold over shuffle distribution
%           occurrence_rd_dist, neurons x behaviors x nrand matrix, shuffled distribution of activation occurrence values


%init output
occurrence = NaN(size(calcium_data,1),size(behav_data,1));
occurrence_rd_dist = NaN(size(calcium_data,1),size(behav_data,1),nrand);

%select commun time bins between calcium recording and behavior recording
[~,ind_cal,ind_loco] = intersect(round(calcium_time*calcium_fps)/calcium_fps,round(behav_time*behav_fps)/behav_fps);
event = calcium_data(:,ind_cal);
behav = behav_data(:,ind_loco);

%compute for each behaviors how often neuronal activation for all bouts
for bb=1:1:12
    duration = sum(behav(bb,:)==1)/calcium_fps;
    onsets = find(diff(behav(bb,:))==1)+1;
    offsets = find(diff(behav(bb,:))==-1);
    if onsets(1)>offsets(1)
        onsets = [1,onsets];
    end
    if onsets(end)>offsets(end)
        offsets = [offsets,length(behav(bb,:))];
    end
    listbouts = [onsets;offsets];
    freqbouts = NaN(size(event,1),size(listbouts,2));
    recbehav = zeros(size(behav(bb,:)));
    for zz=1:1:size(listbouts,2)
        recbehav(listbouts(1,zz):listbouts(2,zz)) = 1;
        freqbouts(:,zz) = nanmean(event(:,(listbouts(1,zz):listbouts(2,zz))),2);
    end
    occurrence(:,bb) = sum(freqbouts>0,2)/size(listbouts,2);
end
                            
%shuffle                
for rr=1:1:nrand
    irand = randperm(size(behav,2));
    rbehav = behav(:,irand);
    for bb=1:1:length(behaviors_list)
        onsets = find(diff(rbehav(bb,:))==1)+1;
        offsets = find(diff(rbehav(bb,:))==-1);
        if onsets(1)>offsets(1)
            onsets = [1,onsets];
        end
        if onsets(end)>offsets(end)
            offsets = [offsets,length(rbehav(bb,:))];
        end
        listbouts = [onsets;offsets];
        freqbouts = NaN(size(event,1),size(listbouts,2));
        recbehav = zeros(size(rbehav(bb,:)));
        for zz=1:1:size(listbouts,2)
            recbehav(listbouts(1,zz):listbouts(2,zz)) = 1;
            freqbouts(:,zz) = nanmean(event(:,(listbouts(1,zz):listbouts(2,zz))),2);
        end
        occurrence_rd_dist(:,bb,rr) = sum(freqbouts==0,2)/size(listbouts,2);
    end
end
%calculate threshold values                                 
occurrence_rd_threshold = quantile(occurrence_rd_dist,percentile_threshold,3);

end