
N_iteration=1000;

X=load('../result/allsbjt_ratio.mat');
X=X.X;


%%%%%%%%%% zero mean and whiten the data %%%%%%%%%%%%%%% 
X=X-mean(X,2)*ones(1,size(X,2));
% %svd
% Cx=X*X'/(size(X,2)-1);
% [U,S,V]=svd(Cx);  
% P=inv(V');
% P1=inv(sqrt(S(:,:)))*P(:,1:5);
% XX=P1'*X;
% %svd

%pca
X=X';
[eigenVector,score,eigenvalue,tsquare] = princomp(X);  %eigenvalue is sorted
trans = eigenVector(:,1:3);     % choose 3

XX = X * trans;
XX=XX';
%pca

%%%%%%%%% Hierarchical clustering %%%%%%%%%%%%%%%%%%
N_clusters=2;

Y = pdist(XX', 'euclidean');    %Generate the dissimilarity matrix, comment out lines 120~133
Z = linkage(Y, 'ward');    %Create hierarchical cluster tree.
cluster_label = cluster(Z, 'maxclust', N_clusters);  %Construct clusters from a hierarchical cluster tree

%%%%%%%%%%%%%%%%%%%% initialization based on the hierarchical clustering %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dim=size(XX,1);
N_points=size(XX,2);
for i=1:N_clusters
    ind=find(cluster_label==i);
    ROI=XX(:,ind);
	mu_init(:,i)=mean(ROI,2); 
	sigma_init(:,:,i)=cov(ROI');
    pi_init(i)=length(ind)/N_points;
end

[R_old,mu_old,sigma_old,pi_old,log_lik_list,count]=EM_MoG(XX,mu_init,sigma_init,pi_init,dim,N_clusters,N_points,N_iteration);






