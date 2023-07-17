function [similarity,similarity_shuffle,similarity_neighbor] = neuronal_similarity_longitudinal(calcium_data_1,calcium_time_1,...
                                                    calcium_data_2,calcium_time_2,long_reg,...
                                                    centroid_location_1, centroid_location_2,calcium_fps,...
                                                    behav_data_1,behav_time_1,behav_data_2,behav_time_2,...
                                                    behav_fps,min_behav_duration)

% This function calculates the similarity in average neuronal activation 
% for each behavior between two recording sessions in the same mouse. 
% Uses Eucliden distance 

% INPUTS:   calcium_data_1, neurons x time matrix session 1
%           calcium_time_1, 1 x time matrix
%           calcium_data_2, neurons x time matrix session 2
%           calcium_time_2, 1 x time matrix
%           long_reg, neurons x 2: 1st col for session 1, 2nd col for session 2, contains indices of neurons found in both session of in only one session
%           centroid_location_1, neurons x 2, position of neurons in the FOV in session 1
%           centroid_location_2, neurons x 2, position of neurons in the FOV in session 2
%           calcium_fps, acquisition rate of endoscope camera
%           behav_data_1, behaviors x time matrix session 1
%           behav_time_1, 1 x time matrix
%           behav_data_2, behaviors x time matrix session2
%           behav_time_2, 1 x time matrix
%           behav_fps, acquisition rate of behavior camera
%           min_behav_duration, minimal duration of behavior in sec
%           
%
% OUTPUTS:  similarity, behaviors x 1 matrix, neuronal activation similarity between two recording sessions using long_reg registration
%           similarity_shuffle, behaviors x 1 matrix, control for neuronal activation similarity using random shuffle of cell pairs
%           similarity_neighbor, behaviors x 1 matrix, control for neuronal activation similarity using the closest neighbor of the actually registered neuron


%init output
similarity = NaN(size(behav_data_1,1),1);
similarity_shuffle = NaN(size(behav_data_1,1),1);
similarity_neighbor = NaN(size(behav_data_1,1),1);

%select commun time bins between calcium recording and behavior recording
[~,ind_cal1,ind_loco1] = intersect(round(calcium_time_1*calcium_fps)/calcium_fps,round(behav_time_1*behav_fps)/behav_fps);
event1 = calcium_data_1(:,ind_cal1);
behav1 = behav_data_1(:,ind_loco1);
[~,ind_cal2,ind_loco2] = intersect(round(calcium_time_2*calcium_fps)/calcium_fps,round(behav_time_2*behav_fps)/behav_fps);
event2 = calcium_data_2(:,ind_cal2);
behav2 = behav_data_2(:,ind_loco2);

%get indices of neurons found in both sessions
LR = long_reg(long_reg(:,1)>0&long_reg(:,2)>0);
                        
%compute similarit
for ee=1:1:size(behav_data_1,1)
    behav2use1 = behav1(ee,:); event2use1 = event1; 
    behav2use2 = behav2(ee,:); event2use2 = event2; 
    if sum(behav2use1==1)>=min_behav_duration*calcium_fps && sum(behav2use2==1)>=min_behav_duration*calcium_fps
        A = mean(event2use1(LR(:,1),behav2use1==1),2,'omitnan')*calcium_fps;
        B = mean(event2use2(LR(:,2),behav2use1==1),2,'omitnan')*calcium_fps;
        %compute similarity
        if ~isempty(A) && ~isempty(B)
            similarity(ee) = (-norm(A/norm(A)-B/norm(B)));
            %cell shuffle / permute cells in B
            rand_shuff = NaN(1,10);
            for rr=1:1:10
                Br = mean(event2use1(randperm(size(event2use1,1),'omitnan'),behav2use1==1&time2keep1==0),2)*calcium_fps;
                if ~isempty(A) && ~isempty(Br)
                    rand_shuff(rr) = (-norm(A/norm(A)-Br/norm(Br)));
                end
                similarity_shuffle(ee) = mean(rand_shuff,'omitnan');
            end
        end
        %closest neighbor control
        closepairs = [LR(:,1),NaN(size(LR(:,2)))];
        for kk=1:1:size(closepairs,1)
            %compute distances
            XY = centroid_location_1-repmat(centroid_location_2,size(centroid_location_1,1),1);
            dists = sqrt(XY(:,1).^2+XY(:,2).^2);
            ind_cell = find(dists==min(dists),1,'first');
            if ind_cell==LR(kk,2)
                dists(ind_cell) = Inf;
                ind_cell = find(dists==min(dists),1,'first');
            end
            closepairs(kk,2) = ind_cell;
        end
        Bc = mean(event2use2(closepairs(:,2),behav2use1==1),2,'omitnan')*calcium_fps;        
        similarity_neighbor(ee) = (-norm(A/norm(A)-Bc/norm(Bc)));
    end
end

end