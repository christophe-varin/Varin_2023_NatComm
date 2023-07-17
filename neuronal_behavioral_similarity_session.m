function [similarity_map_neuron,similarity_map_behav,corr_coeff,corr_pval] = neuronal_behavioral_similarity_session(calcium_data,calcium_time,calcium_fps,...
                                                    behav_data,behav_time,behav_fps,...
                                                    behav_distrib_feat,min_behav_duration)

% This function calculates the similarity in average neuronal activation
% between each pair of behaviors and the behavioral distance between each
% pair of behaviors using the Wassenstein distance (aka earth moving
% distance) and calculate the correlation coefficient between behavioral
% similairty and neuronal activation distance

% INPUTS:   calcium_data, neurons x time matrix
%           calcium_time, 1 x time matrix
%           calcium_fps, acquisition rate of endoscope camera
%           behav_data, behaviors x time matrix
%           behav_time, 1 x time matrix
%           behav_distrib_feat, features x time matrix
%           behav_fps, acquisition rate of behavior camera
%           min_behav_duration, minimal duration of behavior in sec
%
% OUTPUTS:  similarity_map_neuron, behavior x behavior matrix, neuronal activation similarity for each pair of behaviors
%           similarity_map_behav, behavior x behavior matrix, behavior similarity for each pair of behaviors
%           corr_coeff, 1 x 1 double, correlation coefficient (Spearman) between pairwise neuronal similarity and pairwise behavioral similarity
%           corr_pval, 1 x 1 double, associated p-value


%init output
similarity_map_neuron = NaN(size(behav_data,1),size(behav_data,1));
similarity_map_behav = NaN(size(behav_data,1),size(behav_data,1));

%select commun time bins between calcium recording and behavior recording
[~,ind_cal,ind_loco] = intersect(round(calcium_time*calcium_fps)/calcium_fps,round(behav_time*behav_fps)/behav_fps);
event = calcium_data(:,ind_cal);
behav = behav_data(:,ind_loco);
feature = behav_distrib_feat(:,ind_loco);
                        
%compute similarity for each pair of behaviors
for bb1=1:1:size(behav_data,1)
    behav2use1 = behav(bb1,:); event2use = event; 
    for bb2=1:1:size(behav_data,1)
        behav2use2 = behav(bb2,:);
        if sum(behav2use1==1)>=min_behav_duration*calcium_fps && sum(behav2use2==1)>=min_behav_duration*calcium_fps
            A = mean(event2use(:,behav2use1==1),2,'omitnan')*calcium_fps;
            B = mean(event2use(:,behav2use2==1),2,'omitnan')*calcium_fps;
            %compute neuronal similarity
            if ~isempty(A) && ~isempty(B)
                similarity_map_neuron(bb1,bb2) = (-norm(A/norm(A)-B/norm(B)));
            end
            %compute behavior similarity, Wasserstein distance
            WD = 0;
            if bb1~=bb2
                for ff=1:1:size(feature,1)
                    D1 = feature(ff,behav2use1==1);
                    D2 = feature(ff,behav2use2==1);
                    WD = WD+(Wasserstein_Dist(D1',D2'));
                end
                similarity_map_behav(bb1,bb2) = -WD; 
            end
        end
    end
end
%compute correlation coefficient & pvalue
B = similarity_map_behav(:); C = similarity_map_neuron(:);
%remove diagonal elements 
linrem = isnan(B)|isnan(C); B(linrem) = []; C(linrem) = [];
[corr_coeff,corr_pval] = corr(B,C,'Type','Spearman');
end


function WS_Dist = Wasserstein_Dist(XX,YY)
YY(~any(~isnan(YY), 2),:)=[];
XX(~any(~isnan(XX), 2),:)=[];
WS_Dist_temp = NaN(size(XX,2),1);
for jj = 1:size(XX,2)
    X = XX(:,jj);
    Y = YY(:,jj);
    nx = length(X);
    ny = length(Y);
    n = nx + ny;
    XY = [X;Y];
    X2 = [(1/nx).*ones(nx,1);zeros(ny,1)];
    Y2 = [zeros(nx,1);(1/ny).*ones(ny,1)];
    [SortedXY ,I] = sort(XY);
    X2_Sorted = X2(I);
    Y2_Sorted = Y2(I);
    Res = 0;
    E_CDF = 0;
    F_CDF = 0;
    power = 1;
    for ii = 1:n-1
      E_CDF = E_CDF + X2_Sorted(ii);
      F_CDF = F_CDF + Y2_Sorted(ii);
      height = abs(F_CDF-E_CDF);
      width = SortedXY(ii+1) - SortedXY(ii);
      Res = Res + (height ^ power) * width;  
    end
    WS_Dist_temp(jj) = Res;  
end
WS_Dist = max(WS_Dist_temp);
end
