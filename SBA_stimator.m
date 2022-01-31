clc;
clear;
tic
rng(123456);
%format long

%DGP
N = 5000;
k=100;
K = 2;
mu=0;
sigma = 1;
sigma2 =1;
p=3;
kk = 1;
beta_in = zeros(kk,8);
beta_out = zeros(kk,8);% record of each beta
www=100000;
mumu = zeros(1,9);
sigmasigma = eye(9);
sigmasigma(2,3) = 0;
sigmasigma(3,2) = 0;
for i =1:kk
%x =  -1 + 2*rand(N,2);
%x = random('normal',mu,sigma,[N,9]) ; 
x = mvnrnd(mumu,sigmasigma,N);
%x(:,1) = pearsrnd(0,1,-12,150,N,1);
x(:,1) = 1;
x_continuous = x(:,5:9); % continuous
%x = 1 * (x > 0); % discrete
x(:,5:9) = x_continuous;
%e =1*(trnd(1,N,1)>0.2).*( rand(N,1)+0.5)*2 + (1-1*(trnd(1,N,1)>0.8)).*rand(N,1); 
%e = random('normal',mu,sigma,[N,1]) ;
e = random('logistic',mu,sigma2,[N,1]) ;
%e =1*(trnd(1,N,1)>0.2).* random('normal',-2,sigma2,[N,1]) + (1-1*(trnd(1,N,1)>0.8)).*random('normal',5,sigma2,[N,1]); 
%e = trnd(1,N,1);
beta = [0,1,2,4,5,-1,-2,-4,-5];

z = sum(x * beta',2);
y = 1 *(z > e);
xx = x'*x/N;
xx_use = xx^-1;
xx_use = xx_use(2:9,2:9);
L = 1/(2*pi*sigma2^2)^0.5;
w = 1/N* ones(N,1);

%estimation
% First, use logistic to estimate
%beta_0 = [1,0,0,0,0,0,0,0,0];
%beta_0 = [1,1,2,4,5,-1,-2,-4,-5];
beta_0 = [1,1,-100,100,100,1,1,-100,1];
%beta_0 = [100,100,-100,100,100,100,100,-100,-100];
beta_t =beta_0(2:9);
%y_hat = exp(sum(x * beta_0',2))./(1+exp(sum(x * beta_0',2)));
%beta_t = beta_t - (1/L * xx_use /N * x(:,2:9)'*(y_hat -y))';
%iteration
 options1 = optimset('MaxFunEvals',1000000,'MaxIter',10000);
 beta_ave = 0;
for j =1:k
    x_trans = x(:,2:9) * beta_t' ;
    x_trans = x_trans/www;
    x_reg =[x_trans.^0,x_trans];%1111
    %x_reg =x_trans;
    for ii =2:p
    %for ii = 2:p
        x_reg = [x_reg,x_trans.^ii];
    end
    fun =@(S)sum(log(1 +  exp(x_reg *S'))) - y'* x_reg*S';
    S0 = zeros(1,p+1);%2222
    %S0 = zeros(1,p);%2222
    S0(1) =1;
    St= fminsearch(fun,S0,options1);
    g_trans= x_reg * St';
    g_hat = exp(g_trans)./(1+exp(g_trans));
    beta_t1 = beta_t - (200/L * xx_use/j^.6 /N * x(:,2:9)'*(g_hat -y))';
    %if max(abs(beta_t1 - beta_t)) <0.000001 % origial #5 0's
    %    break
    %end
    beta_t = beta_t1;
    beta_t/beta_t(1)
    sum(log(1 +  exp(x_reg *St'))) - y'* x_reg*St'
    beta_ave = (beta_ave*(j-1) + beta_t1)/j;
    beta_ave/beta_ave(1)
    beta_show(j,:) = beta_t;
end
beta_in(i,:) =beta_t/beta_t(1);
beta_out(i,:) = beta_ave/beta_ave(1);

end

xx_vuse = x(:,2:9)'* x(:,2:9)/N;


 x_trans = x(:,2:9) * beta_ave' ;
    x_trans = x_trans/www;
    x_reg =[x_trans.^0,x_trans];%1111
    %x_reg =x_trans;
    for ii =2:p
    %for ii = 2:p
        x_reg = [x_reg,x_trans.^ii];
    end
     fun =@(S)sum(log(1 +  exp(x_reg *S'))) - y'* x_reg*S';
    S0 = zeros(1,p+1);%2222
    %S0 = zeros(1,p);%2222
    S0(1) =1;
    St= fminsearch(fun,S0,options1);
g_trans= x_reg * St';
    g_hat = exp(g_trans)./(1+exp(g_trans));

sigma11 = zeros(8,8);
sigma22 = zeros(8,8);
sigma33 = zeros(8,8);
E_trans = zeros(p+1,8);
%for i =1:N
%sigma11 = sigma11 + g_hat(i)*(1-g_hat(i))*(xx_use^0.5*x(i,2:9)' +x(i,2:9)*beta_ave'* lbeta)*(xx_use^0.5*x(i,2:9)'+x(i,2:9)*beta_ave'* lbeta)';
%end

%for i =1:N
%sigma11 = sigma11 + g_hat(i)*(1-g_hat(i))*(xx_use^-0.5*x(i,2:9)')*(xx_use^-0.5*x(i,2:9)')';
%end

for i =1:N
sigma11 = sigma11 + g_hat(i)*(1-g_hat(i))*x(i,2:9)'*x(i,2:9);
end
sigma11 = sigma11/N



g_prime_trans = x_reg(:,1:3)*[St(2)/www,2*St(3)/www,3*St(4)/www]';
g_prime_hat = exp(g_trans)./(1+exp(g_trans)).^2 .*g_prime_trans;
 g_l_hat = exp(g_trans)./(1+exp(g_trans)).^2;
 gl_hat = g_l_hat *ones(1,p+1);



for i =1:N
    E_trans = E_trans +g_prime_hat(i)* x_reg(i,:)'*x(i,2:9);
end
E_trans = E_trans/N;

for i=1:N
    sigma22 = sigma22 + g_prime_hat(i)*x(i,2:9)'*x(i,2:9) -x(i,2:9)'*x_reg(i,:)* (x_reg'*x_reg/N)^-1 * E_trans;
end
sigma22 = sigma22/N 

%for i=1:N
%   sigma22 = sigma22 + g_prime_hat(i)* xx_vuse^-0.5*x(i,2:9)'*x(i,2:9);
%end
for i=1:N
    sigma33 = sigma33 + g_prime_hat(i)*x(i,2:9)'*x(i,2:9);
end
sigma33 = sigma33/N                                                                                                                                                                                                                             

sigma_beta_ave = sigma22(2:8,2:8)'^-1 *sigma11(2:8,2:8)*sigma22(2:8,2:8)^-1 
st_beta_ave = (diag(sigma_beta_ave).^.5)'/sqrt(N)
z_beta_ave = beta_ave(2:8)./st_beta_ave

sigma_beta_ave = sigma33(2:8,2:8)'^-1 *sigma11(2:8,2:8)*sigma33(2:8,2:8)^-1 
st_beta_ave = (diag(sigma_beta_ave).^.5)'/sqrt(N)
z_beta_ave = beta_ave(2:8)./st_beta_ave

%aver_itime = mean(itime);
%mean_bias_1_10000 = mean(beta_in - ones(kk,1) * beta(2:9));
%rmse_1_10000 = (mean((beta_in - ones(kk,1) * beta(2:9)).^2)).^0.5;
%save('stochasticgradient_discrete_5000.mat','mean_bias_1_5000','rmse_1_5000','aver_itime');
%mean_bias_2_10000 = mean(beta_out - ones(kk,1) * beta(2:9));
%rmse_2_10000 = (mean((beta_out - ones(kk,1) * beta(2:9)).^2)).^0.5;
toc
scatter( x_reg(:,2)*www/beta_t(1),g_hat);
hold on;
scatter(z,exp(z)./(1+exp(z)))
hold off

