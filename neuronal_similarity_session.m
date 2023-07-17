function [similarity,similarity_shuffle] = neuronal_similarity_session(calcium_data,calcium_time,calcium_fps,behav_data,behav_time,behav_fps,min_behav_duration,type_str)

% This function calculates the similarity in average neuronal activation 
% for each behavior. Uses Eucliden distance or dot product between first
% and second halves of recording or inverse of coefficient of variation for
% the whole recording as specified with typr_str string

% INPUTS:   calcium_data, neurons x time matrix
%           calcium_time, 1 x time matrix
%           calcium_fps, acquisition rate of endoscope camera
%           behav_data, behaviors x time matrix
%           behav_time, 1 x time matrix
%           behav_fps, acquisition rate of behavior camera
%           min_behav_duration, minimal duration of behavior in sec
%           type_str, string to set type of similarity computation
%                       can be 'euclidean', euclidean similarity between first and second halves of recording
%                              'dot', dot product between first and second halves of recording
%                              'CVinv', inverse of coefficient of variation, all recording long
%
% OUTPUTS:  similarity, behaviors x 1 matrix, neuronal activation similarity calculated according with type_str
%           similarity_shuffle, behaviors x 1 matrix, shuffle control for neuronal activation similarity not available for inverse of CV


%init output
similarity = NaN(size(behav_data,1),1);
similarity_shuffle = NaN(size(behav_data,1),1);

%select commun time bins between calcium recording and behavior recording
[t,ind_cal,ind_loco] = intersect(round(calcium_time*calcium_fps)/calcium_fps,round(behav_time*behav_fps)/behav_fps);
event = calcium_data(:,ind_cal);
behav = behav_data(:,ind_loco);
                        
%compute similarity according to type_str
for ee=1:1:size(behav_data,1)
    behav2use = behav(ee,:); event2use = event; 
    switch type_str
        %for euclidean distance or dot product
        case {'euclidean','dot'}
            %split duration in half
            t_min = min(t); t_max = max(t); wind = t_min+(t_max-t_min)/2;
            time2keep1 = double(t<=wind);
            if sum(behav2use==1&time2keep1==1)>=min_behav_duration*calcium_fps && sum(behav2use==1&time2keep1==0)>=min_behav_duration*calcium_fps
                A = mean(event2use(:,behav2use==1&time2keep1==1),2,'omitnan')*calcium_fps;
                B = mean(event2use(:,behav2use==1&time2keep1==0),2,'omitnan')*calcium_fps;
                %compute similarity
                if ~isempty(A) && ~isempty(B)
                    if strcmp(type_str,'euclidean')
                        similarity(ee) = (-norm(A/norm(A)-B/norm(B)));
                    elseif strcmp(type_str,'dot')
                        similarity(ee) = dot(A,B);
                    end
                    %cell shuffle / permute cells in B
                    rand_shuff = NaN(1,10);
                    for rr=1:1:10
                        Br = mean(event2use(randperm(size(event2use,1),'omitnan'),behav2use==1&time2keep1==0),2)*calcium_fps;
                        if ~isempty(A) && ~isempty(Br)
                            if strcmp(type_str,'euclidean')
                                rand_shuff(rr) = (-norm(A/norm(A)-Br/norm(Br)));
                            elseif strcmp(type_str,'dot')
                                rand_shuff(rr) = dot(A,Br);  
                            end
                        end
                        similarity_shuffle(ee) = mean(rand_shuff,'omitnan');
                    end
                end
            end
        %for inverse coefficient of variation
        case 'CV'
            if  sum(behav2use==1)>=min_behav_duration
                CVinv = mean(event2use(:,behav2use==1==1),2,'omitnan')./std(event2use(:,behav2use==1==1),[],2,'omitnan');
                CVinv(mean(event2use(:,behav2use==1==1),2,'omitnan')==0) = 0; 
                similarity(ee) = nanmean(CVinv);
            end
    end
end

end