function [R_old,mu_old,sigma_old,pi_old,log_lik_list,count]=EM_k_mean(Data,mu_init,sigma_init,pi_init,dim,N_clusters,N_points,N_iteration)

% INPUT:
% Data represents reductive sensor signals and its size is dim x N_points
% For example : 10 x 5439
% mu_init represents initial mu of MOG and its size is dim x N_clusters
% For example : 10 x 7
% sigma_init represents initial sigma of MOG and its size is dim x dim x N_clusters
% For example : 10 x 10 x 7
% pi_init represents initial pi of MOG and its size is 1 x N_clusters
% For example : 1 x 7
% N_iteration represents repeated numbers
% For example : 1000
%
% OUTPUT:
% R_old represents estimated responsabilty matrix from EM algorithm and its size is N_clusters x N_points
% For example : 7 x 5439
% mu_old represents estimated mu of MOG and its size is dim x N_clusters
% For example : 10 x 7
% sigma_old represents estimated sigma of MOG and its size is dim x dim x N_clusters
% For example : 10 x 10 x 7
% pi_old represents estimated pi of MOG and its size is 1 x N_clusters
% For example : 1 x 7
% log_lik_list represents log likelihood of p(Xn|mu,sigma,pi) obtained at
% every iteration
% For example : 1 x count
%
% This is a demonstratable verison. Copyright Yen-Chun Chou 04/17/2007

mu_old = mu_init;
sigma_old = sigma_init;
pi_old = pi_init;
mu_new = zeros(dim,N_clusters);
sigma_new = zeros(dim,dim,N_clusters);
pi_new = ones(1,N_clusters);

count = 0;
R_old = zeros(N_clusters,N_points);
num = zeros(N_clusters,N_points);

while count < N_iteration
    
    %%%%%%%%%% compute the responsibility matrix
    
    for i = 1:N_clusters
        mu = mu_old(:,i)*ones(1,N_points);
        tmp1 = sigma_old(:,:,i);
        tmp2 = Data-mu;                        % Importance: if Data is bad, it will cause total_pr=0
        tmp3 = (sum((tmp2'*inv(tmp1)).*tmp2',2))';
        num(i,:) = 1/((2*pi)^(dim/2)*sqrt(det(tmp1)))*exp(-0.5*tmp3);
        num(i,:) = num(i,:)*pi_old(i);
        % num(i,:) = mvnpdf(Data',mu_old(:,i)',sigma_old(:,:,i))';
        % num(i,:) = num(i,:)*pi_old(i);
    end
    clear tmp1 tmp2 tmp3
    
    R_old = zeros(N_clusters,N_points);
    [C,I] = max(num);
    for i = 1:N_clusters
        ind = find(I == i);
        R_old(i,ind) = 1;
    end
    count = count+1;
    log_lik_list(count) = sum(log(C));
    
    %%%%%%%% update mean, covariance, and pi for each cluster %%%%%%%%%%%%%%
    
    for i = 1:N_clusters
        R_old2 = ones(dim,1)*R_old(i,:); % dim * N_points
        Data_R = R_old2.*Data;
	    num1 = sum(Data_R,2);
	    den = sum(R_old(i,:));
	    mu_new(:,i) = num1/den;
        num2 = Data_R*Data';
        % sigma_new(:,:,i) = num2/den-mu_new(:,i)*mu_new(:,i)';
	    sigma_new(:,:,i) = num2/den-mu_new(:,i)*mu_new(:,i)'+eye(dim)*10E-6;
        pi_new(i) = den/N_points;
    end
    
    %%%%%%%%%%% update parameters
    
    mu_old = mu_new;
    sigma_old = sigma_new;
    pi_old = pi_new;
    
    %%%%%%%%%%%%% compute and show log-likelihood
    
    if (count>1) && (log_lik_list(count)-log_lik_list(count-1)<1e-8)
        fprintf('Iteration = %g\nBye!\n',count);
        break;
    end
end