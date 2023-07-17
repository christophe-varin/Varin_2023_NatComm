function [mean_population,mean_cells,mean_behaviors] = neuronal_activity_session(calcium_data,calcium_time,calcium_fps,behav_data,behav_time,behav_fps,min_behav_duration)

% This function calculates the average activity for the whold recording and
% during each behavior.

% INPUTS:   calcium_data, neurons x time matrix
%           calcium_time, 1 x time matrix
%           calcium_fps, acquisition rate of endoscope camera
%           behav_data, behaviors x time matrix
%           behav_time, 1 x time matrix
%           behav_fps, acquisition rate of behavior camera
%           min_behav_duration, minimal duration of behavior in sec
%
% OUTPUTS:  mean_population, 1 x 1 matrix, average population activity for all behaviors for all cells
%           mean_cells, neurons x 1 matrix, average population activity for all behaviors
%           mean_behaviors, behaviors x 1, average population activity for all cells for each behavior


%init output
mean_behaviors = NaN(size(behav_data,1),1);

%select commun time bins between calcium recording and behavior recording
[~,ind_cal,ind_loco] = intersect(round(calcium_time*calcium_fps)/calcium_fps,round(behav_time*behav_fps)/behav_fps);
event = calcium_data(:,ind_cal);
behav = behav_data(:,ind_loco);

%compute average activity for the whole recording
mean_population = mean(mean(event,2,'omitnan'),1,'omitnan')*calcium_fps;
mean_cells = mean(event,2,'omitnan')*calcium_fps;

%compute average activity during each behavior
for bb=1:1:size(behav,1)
    dur = sum(behav(bb,:))/fps;
    if dur>=min_behav_duration
        tp_behav = behav(bb,:);
        mean_behaviors(bb) = mean(mean(event(:,tp_behav'==1),2,'omitnan'),1,'omitnan')*calcium_fps;
    else
        disp([' behav ',num2str(bb),' removed / duration ',num2str(dur)])
    end
end

end