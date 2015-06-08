clear
close all
%colordef black;
N_iteration=1000;
color = ['g' 'r' 'm' 'y' 'c' 'b' 'w' ];

load yhkao_roi.mat
X=mix_matrix;
load yhkao_roi_ind.mat
t=im_id;

%%%%%%%%%% zero mean and whiten the data %%%%%%%%%%%%%%% 
X=X-mean(X,2)*ones(1,size(X,2));
Cx=X*X'/(size(X,2)-1);
[U,S,V]=svd(Cx);
P=inv(V');
P1=inv(sqrt(S(:,:)))*P(:,1:5);
XX=P1'*X;

%%%%%%%%% Hierarchical clustering %%%%%%%%%%%%%%%%%%
N_clusters=7;

Y = pdist(XX', 'euclidean');    %Generate the dissimilarity matrix, comment out lines 120~133
Z = linkage(Y, 'ward');    %Create hierarchical cluster tree.
cluster_label = cluster(Z, 'maxclust', N_clusters);  %Construct clusters from a hierarchical cluster tree

figure(1)

rr=zeros(128,128);rr(im_id(find(cluster_label==1)))=1;subplot(2,5,1);imshow(rr);
title('class #1_r')
rr=zeros(128,128);rr(im_id(find(cluster_label==2)))=1;subplot(2,5,2);imshow(rr);
title('class #2_g')
rr=zeros(128,128);rr(im_id(find(cluster_label==3)))=1;subplot(2,5,3);imshow(rr);
title('class #3_b')
rr=zeros(128,128);rr(im_id(find(cluster_label==4)))=1;subplot(2,5,4);imshow(rr);
title('class #4_k')
rr=zeros(128,128);rr(im_id(find(cluster_label==5)))=1;subplot(2,5,5);imshow(rr);
title('class #5_m')
rr=zeros(128,128);rr(im_id(find(cluster_label==6)))=1;subplot(2,5,6);imshow(rr);
title('class #6_y')
rr=zeros(128,128);rr(im_id(find(cluster_label==7)))=1;subplot(2,5,7);imshow(rr);
title('class #7_c')

figure(2)
for i=1:N_clusters
    ind=find(cluster_label==i);
    curve(i,:)=mean(mix_matrix(:,ind),2);
end

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

if 0
% rnndomly assume the number of elements in each cluser 
pis=rand(N_clusters,1);
sum_pis=sum(pis);
pis_n=pis/sum_pis;
for i=1:N_clusters-1
  fr(i) = round(pis_n(i)*N_points);
end
tmp_sum=sum(fr);
fr(N_clusters)=N_points-tmp_sum;

list_cl = [0 cumsum(fr)];

%%% init mu, cov, pi
mu_init=zeros(dim,N_clusters);
sigma_init=zeros(dim,dim,N_clusters);

rand_idx=randperm(N_points); %Random permutation.

for i=1:N_clusters,
    rand_X=XX(:,rand_idx(list_cl(i)+1:list_cl(i+1)));
	mu_init(:,i)=mean(rand_X,2);   	
	sigma_init(:,:,i)=cov(rand_X'); 	
end;
pi_init=ones(1,N_clusters)/N_clusters;

mu_old=mu_init;
sigma_old=sigma_init;
pi_old=pi_init;
mu_new = zeros(dim,N_clusters);
sigma_new = zeros(dim,dim,N_clusters);
pi_new=ones(1,N_clusters);

count=0;
R_old=zeros(N_clusters,N_points);
num=zeros(N_clusters,N_points);
tot_pr=zeros(N_clusters,N_points);
while (count<N_iteration),      
      %%%%%%%%%% compute the responsabilty matrix     
      for i=1:N_clusters
          mu=mu_old(:,i)*ones(1,size(XX,2));
          tmp1=sigma_old(:,:,i);
          tmp2=XX-mu;
          tmp3=sum((tmp2'*inv(tmp1)).*tmp2',2);
          p(i,:)=tmp3';
          num(i,:)=1/((2*pi)^(dim/2)*sqrt(det(tmp1)))*exp(-0.5*p(i,:));
          num(i,:)=num(i,:)*pi_old(i);
      end
      total_pr=sum(num); % sum along row direction and become an 1 x N_points vector
      total_pr=ones(N_clusters,1)*total_pr; % now it becomes N_clusters x N_points
%      R_old= num./total_pr;                 % it is N_clusters x N_points
      R_old=zeros(N_clusters,N_points);
      [Y,IND]=max(num./total_pr,[],1);
      for k=1:N_points
        R_old(IND(k),k)=1;
      end      
      count=count+1;
      log_lik_list(count)=  sum(log(sum(num)));

      if 0
       fprintf('LogLik( %d ) = %20.8f\n',count, log_lik_list(count)); 
      end

      %%%%%%%% update mean and covariance for each cluster %%%%%%%%%%%%%%    
      
      for (i=1:N_clusters),
	    R_old2=ones(dim,1)*R_old(i,:); % dim * N_points
        XX_R = R_old2.*XX;
	    num1=sum(XX_R,2);
	    den=sum(R_old(i,:));
	    mu_new(:,i)=num1/den;        
        num2=(XX_R*XX');      
	    sigma_new(:,:,i) = num2/den - mu_new(:,i)*mu_new(:,i)'+eye(dim)*10E-6;

	    %%%%% compute pi
	    pi_new(i) = den/N_points;

	    %%%%%show figures
	    
      end;
      
      %%%%%%%%%%% update parameters
      mu_old = mu_new;
      sigma_old = sigma_new;
      pi_old = pi_new;
      
      %%%%%%%%%%%%% compute and show loglikelihhod
            
      if (count>1)&(log_lik_list(count)<log_lik_list(count-1)),
	 %fprintf('ERROR!\n');
      end;
      
      if (count>1) & (log_lik_list(count)-log_lik_list(count-1)<10E-8),
    	 fprintf('Done. \nBye!\n');
        break;
      end;
      

end; %% end while loop

end % if 0

figure(3);
title('Loglikelihood');
plot(log_lik_list,'-');
drawnow;


%AIC=-2*log_lik_list(count)+2*(dim*N_clusters+dim*(dim+1)/2*N_clusters+N_clusters-1)
%MDL=-log_lik_list(count)+2.5*(dim*N_clusters+dim*(dim+1)/2*N_clusters+N_clusters-1)*log(N_points)
-log_lik_list(count)
(dim*N_clusters+dim*(dim+1)/2*N_clusters+N_clusters-1)
log(N_points)

MDL=-log_lik_list(count)+0.75*(dim*N_clusters+dim*(dim+1)/2*N_clusters+N_clusters-1)*log(N_points)


%%%%%%%%% compare with Matlab K-means routine

t=im_id;
figure(4)
rr=zeros(128,128);rr(t)=R_old(1,:);subplot(2,4,1);imshow(rr);
title('r')
rr=zeros(128,128);rr(t)=R_old(2,:);subplot(2,4,2);imshow(rr);
title('g')
rr=zeros(128,128);rr(t)=R_old(3,:);subplot(2,4,3);imshow(rr);
title('b')
rr=zeros(128,128);rr(t)=R_old(4,:);subplot(2,4,4);imshow(rr);
title('k')
rr=zeros(128,128);rr(t)=R_old(5,:);subplot(2,4,5);imshow(rr);
title('m')
rr=zeros(128,128);rr(t)=R_old(6,:);subplot(2,4,6);imshow(rr);
title('y')
rr=zeros(128,128);rr(t)=R_old(7,:);subplot(2,4,7);imshow(rr);
title('c')
figure(5)
H=mix_matrix*R_old'*inv(R_old*R_old');
 plot([1:70],H(:,1),'r',[1:70],H(:,2),'g',[1:70],H(:,3),'b',[1:70],H(:,4),'k',[1:70],H(:,5),'m',[1:70],H(:,6),'y',[1:70],H(:,7),'c');

