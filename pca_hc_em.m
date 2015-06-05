clear
close all
%colordef black;
N_iteration=1000;



load guo_roi.mat
X=load('allsbjt_ratio.mat');


%%%%%%%%%% zero mean and whiten the data %%%%%%%%%%%%%%% 
X=X-mean(X,2)*ones(1,size(X,2));
Cx=X*X'/(size(X,2)-1);
[U,S,V]=svd(Cx);    %
P=inv(V');
P1=inv(sqrt(S(:,:)))*P(:,1:10);
XX=P1'*X;

%%%%%%%%% Hierarchical clustering %%%%%%%%%%%%%%%%%%
N_clusters=2;

Y = pdist(XX', 'euclidean');    %Generate the dissimilarity matrix, comment out lines 120~133
Z = linkage(Y, 'ward');    %Create hierarchical cluster tree.
cluster_label = cluster(Z, 'maxclust', N_clusters);  %Construct clusters from a hierarchical cluster tree



plot([1:70],curve(1,:)','r',[1:70],curve(2,:)','g',[1:70],curve(3,:)','b',[1:70],curve(4,:)','k',[1:70],curve(5,:)','m',[1:70],curve(6,:)','y',[1:70],curve(7,:)','c');

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






