function [behav_distrib_feat,gmm_grp,cluster_seeds] = generate_behavior_clusters(DLC_data,subsampling,...
                                                    N_run_tsne,N_run_gmm,min_nb_frames_per_cluster)

% This function performs multiple iterations of tSNE dim reduction followed
% by GMM clustering to behavioral clusters

% INPUTS:   DLC_data, structure generated from DeepLabCut output
%           subsampling, integer, time downsampling
%           N_run_tsne, number of iterations of t-SNE reduction algorithm
%           N_run_gmm, number of iterations of GMM cluster
%                       final number of iterations: N_run_tsne x N_run_gmm
%           min_nb_frames_per_cluster, minimum number of observations for final cluster 
%
% OUTPUTS:  behav_distrib_feat, features x time matrix
%           gmm_grp, all clustering results for all iterations after tSNE + GMM
%           cluster_seeds, final clustering of observations in behav_distrib_feat(1:round(subsampling):end,:)


%get data from DLC_data structure
    %time
time = DLC_data.time;
    %center of mass
center = DLC_data.center; center(center(:,3)<0.9,1:2) = NaN; %filter out points with low likelihood
center(:,1:2) = center(:,1:2)-repmat(DLC_data.centerOF,length(center),1); 
    %neck
neck = DLC_data.neck; neck(neck(:,3)<0.9,1:2) = NaN;
neck(:,1:2) = neck(:,1:2)-repmat(DLC_data.centerOF,length(neck),1); 
    %tail start
tail = DLC_data.tail_start; tail(tail(:,3)<0.8,1:2) = NaN;
tail(:,1:2) = tail(:,1:2)-repmat(DLC_data.centerOF,length(tail),1); 
tail(:,1:2) = fillmissing(tail(:,1:2),'linear','SamplePoints',time);
    %ears
earL = DLC_data.ear_left; earL(earL(:,3)<0.9,1:2) = NaN;
earL(:,1:2) = earL(:,1:2)-repmat(DLC_data.centerOF,length(earL),1); 
earL(:,1:2) = fillmissing(earL(:,1:2),'linear','SamplePoints',time);
earR = DLC_data.ear_right; earR(earR(:,3)<0.9,1:2) = NaN;
earR(:,1:2) = earR(:,1:2)-repmat(DLC_data.centerOF,length(earR),1); 
earR(:,1:2) = fillmissing(earR(:,1:2),'linear','SamplePoints',time);
    %nose
nose = DLC_data.nose; nose(nose(:,3)<0.6,1:2) = NaN;
nose(:,1:2) = nose(:,1:2)-repmat(DLC_data.centerOF,length(nose),1); 
    %camera
camera = DLC_data.camera; camera(camera(:,3)<0.9,1:2) = NaN;
camera(:,1:2) = camera(:,1:2)-repmat(DLC_data.centerOF,length(camera),1);
            
%compute behaviors features
    %average speed 
speed_vec = [[NaN,NaN];[center(3:end,1)-center(1:end-2,1) center(3:end,2)-center(1:end-2,2)]./(time(3:end)-time(1:end-2));[NaN,NaN]];
cos_angle = (speed_vec(:,1).*(neck(:,1)-center(:,1))+speed_vec(:,2).*(neck(:,2)-center(:,2)))./(sqrt(speed_vec(:,1).^2+speed_vec(:,2).^2).*sqrt((neck(:,1)-center(:,1)).^2+(neck(:,2)-center(:,2)).^2));
speed_proj_norm = sqrt(speed_vec(:,1).^2+speed_vec(:,2).^2).*cos_angle;
speed_proj_norm_smooth = smoothdata(speed_proj_norm,'movmean',20);
    %camera movement speed 
speed_proj_vec_smooth = (neck(:,1:2)-center(:,1:2)).*repmat(speed_proj_norm_smooth,1,size(neck(:,1:2),2))./repmat(sqrt((neck(:,1)-center(:,1)).^2+(neck(:,2)-center(:,2)).^2),1,size(neck(:,1:2),2));
speed_cam_vec = [[NaN,NaN];[camera(3:end,1)-camera(1:end-2,1) camera(3:end,2)-camera(1:end-2,2)]./(time(3:end)-time(1:end-2));[NaN,NaN]];
speed_cam_vec2 = (speed_cam_vec-speed_vec);
speed_cam2 = smoothdata(sqrt(speed_cam_vec2(:,1).^2+speed_cam_vec2(:,2).^2),'movmean',20);
    %movement angle 
cos_angle_speed = (speed_proj_vec_smooth(1:end-1,1).*speed_proj_vec_smooth(2:end,1)+speed_proj_vec_smooth(1:end-1,2).*speed_proj_vec_smooth(2:end,2))./speed_proj_norm_smooth(1:end-1)./speed_proj_norm_smooth(2:end); % use cos a = prodiut scalaire / norm vector1 /norm vector 2
sin_angle_speed = (speed_proj_vec_smooth(1:end-1,1).*speed_proj_vec_smooth(2:end,2)-speed_proj_vec_smooth(1:end-1,2).*speed_proj_vec_smooth(2:end,1))./speed_proj_norm_smooth(1:end-1)./speed_proj_norm_smooth(2:end); % use sin a = norme prodiut vectoriel / norm vector 1 / norm vector 2 
angle_speed = [sign(sin_angle_speed).*real(acosd(cos_angle_speed));0];
int_angle_speed = movsum(angle_speed,[0 0.2*fps],'omitnan');
    %neck-center-tail total distance
tot_distance = sqrt((center(:,1)-neck(:,1)).^2+(center(:,2)-neck(:,2)).^2)+sqrt((center(:,1)-tail(:,1)).^2+(center(:,2)-tail(:,2)).^2);
tot_distance(tot_distance>6) = 6;
    %neck elongation
dist_neck_center = sqrt((center(:,1)-neck(:,1)).^2+(center(:,2)-neck(:,2)).^2);
    %distance between neck and orthogonal projection of camera position along the vector orthogonal to ears positions (pretty much collinear to neck-nose vector but always defined)
ortho_ears = [-(earL(:,2)-earR(:,2)),earL(:,1)-earR(:,1)]; ortho_ears = ortho_ears./repmat(sqrt(ortho_ears(:,1).^2+ortho_ears(:,2).^2),1,2);
dproj2 = ((camera(:,1)-neck(:,1)).*(ortho_ears(:,1))+(camera(:,2)-neck(:,2)).*(ortho_ears(:,2)));
dproj2(dproj2<-0.5) = -0.5;

%build features data space
datspace = [speed_proj_norm_smooth,speed_cam2,int_angle_speed,smoothdata(tot_distance,'movmean',20),smoothdata(dist_neck_center,'movmean',20),dproj2];
datspace = fillmissing(datspace,'linear');
behav_distrib_feat = datspace;
datspace = datspace(1:round(subsampling):end,:);
            
%tsne dim reduction
rng(54872643) %for repro
tsne_feats_save = cell(N_run_tsne,1); loss_save = cell(N_run_tsne,1);
p = 100; exag = 20; lr = 2000; %set hyperparameters
[tsne_feats,loss] = tsne(datspace,'Standardize',true,'Exaggeration',exag,'LearnRate',lr,'Perplexity',p,'NumDimensions',3); % Refer to MATLAB tsne function for arguments
tsne_feats_save{1} = tsne_feats; loss_save(1) = loss;
disp(['run ',num2str(1),' / loss ',num2str(loss)])
for runs=2:1:N_run_tsne
    [tsne_feats_tp,loss_tp] = tsne(datspace,'Standardize',true,'Exaggeration',exag,'LearnRate',lr,'Perplexity',p,'NumDimensions',3); % Refer to MATLAB tsne function for arguments
    disp(['run ',num2str(runs),' / loss ',num2str(loss_tp)])
    tsne_feats_save{runs} = tsne_feats; loss_save(runs) = loss;
    if loss_tp<loss
        loss = loss_tp;
        tsne_feats = tsne_feats_tp;
    end
end

%run a Gaussian Mixture Model Expectation Maximization to group the t-SNE clusters
kclass = 50;
grp_all = cell(size(tsne_feats_save,2),N_run_gmm);
for runs=1:1:N_run_tsne
    X = tsne_feats_save{runs}';
    disp(['gmm run ',num2str(runs)])
    [~, ~, ~,~,grp_iters] = em_gmm(X, kclass, N_run_gmm);
    grp_all(runs,:) = grp_iters;
end
gmm_grp = grp_all;

%from grp_iter compute how many times frames clustered together
grp_iters2 = cell2mat(reshape(grp_all,[],1));
count_clust = zeros(size(grp_iters2,2));
for fr1=1:1:size(grp_iters2,2)
    count_clust(fr1,:) = sum(grp_iters2(:,fr1)==grp_iters2,1);
end
count_clust = 1-count_clust/size(grp_iters2,1);
%figure; imagesc(double(count_clust/size(grp_iters2,1)>0.8))

Z = linkage(squareform(count_clust));
cutoff = 0.25;
c = cluster(Z,'cutoff',cutoff);
[C,~,ic] = unique(c);
a_counts = [C,accumarray(ic,1),double(accumarray(ic,1)>=min_nb_frames_per_cluster)];
a_counts(a_counts(:,3)==0,:) = [];
c_new = zeros(size(c));
for cc=1:1:size(a_counts,1)
    c_new(c==a_counts(cc,1))=cc;
end
cluster_seeds = c_new;            
end
    

%Functions to perform GMM

function [grp_max,gmm_max,ll_max,maxll,grp] = em_gmm(datspace, nb_class, n)
% Perform expectation maximization algorithm fitting a mixture of Gaussian distribs

%   datspace, matrix with Rows variables and Columns observations
%   nb_class,maximum number of clusters
%   n, number of random iterations
%
%   grp_max, groups across obersvations
%   gmm_max, Gaussian mixture models with best fitted parameters
%   ll_max, log(likelihood) for the optimal n
%   ll_n, log(likelihood) for the n iterations
%   grp, all realizations of clustering

optloss = 1e-6;
itmax = 500;
grp = cell(n,1); gmm = cell(n,1); lln = cell(n,1); maxll = NaN(n,1);
for n_rand = 1:n
    ll = -inf(1,itmax);
    R = initl(datspace,nb_class);
    for iter = 2:itmax
        [~,grp{n_rand}(1,:)] = max(R,[],2);
        R = R(:,unique(grp{n_rand}));
        gmm{n_rand} = maxim(datspace,R);
        [R, ll(iter)] = expect(datspace,gmm{n_rand});
        if abs(ll(iter)-ll(iter-1)) < optloss*abs(ll(iter)) 
            break; 
        end
    end
    lln{n_rand} = ll(2:iter);
    maxll(n_rand) = max(lln{n_rand});
end
n_max = find(maxll == max(maxll));
grp_max = grp{n_max};
gmm_max = gmm{n_max};
ll_max = lln{n_max};
end

function R = initl(datspace, nb_class)
n = size(datspace,2); 
if numel(nb_class) == 1 
    grp = ceil(nb_class*rand(1,n));
    R = full(sparse(1:n,grp,1,n,nb_class,n));
else
    error('ERROR: init is not valid.');
end
end

function [R, ll] = expect(datspace, gmm)
mu = gmm.mu;
sig = gmm.sig;
w = gmm.w;
n = size(datspace,2);
k = size(mu,2);
R = zeros(n,k);
for i = 1:k
    R(:,i) = lgauspdf(datspace,mu(:,i),sig(:,:,i));
end
R = bsxfun(@plus,R,log(w));
T = logsumexp(R,2);
ll = sum(T)/n; % loglikelihood
R = exp(bsxfun(@minus,R,T));
end

function gmm = maxim(datspace, R)
[d,n] = size(datspace);
k = size(R,2);
nk = sum(R,1);
w = nk/n;
mu = bsxfun(@times, datspace*R, 1./nk);
sig = zeros(d,d,k);
r = sqrt(R);
for i = 1:k
    Xo = bsxfun(@minus,datspace,mu(:,i));
    Xo = bsxfun(@times,Xo,r(:,i)');
    sig(:,:,i) = Xo*Xo'/nk(i)+eye(d)*(1e-6);
end
gmm.mu = mu;
gmm.sig = sig;
gmm.w = w;
end

function y = lgauspdf(datspace, mu, sig)
d = size(datspace,1);
datspace = bsxfun(@minus,datspace,mu); % difference from each data point to the cluster mean
[U,p]= chol(sig); % Cholesky factorization of covariance matrix
if p ~= 0
    error('ERROR: covariance matrix is not positive definite');
end
Q = U'\datspace;
M = dot(Q,Q,1);  % M distance after taking dot product
c = d*log(2*pi)+2*sum(log(diag(U)));   % Gaussian formula for normalization constant
y = -(c+M)/2;
end

function s = logsumexp(X, dim)
y = max(X,[],dim);
s = y+log(sum(exp(bsxfun(@minus,X,y)),dim));
i = isinf(y);
if any(i(:))
    s(i) = y(i);
end
end

