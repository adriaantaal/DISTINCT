

clear all;close all;clc;        
M = 64; % observations
N = 256; % the length of the signal x
K = 10;  %the sparsity of the signal x
Index_K = randperm(N);        
x = zeros(N,1);        
x(Index_K(1:K)) = 5*randn(K,1); %x is sparsely K and the position is random
Psi = eye(N); %defining a sparse matrix as a unit matrix x=Psi*theta
Phi = randn(M,N);% measurement matrix is ​​Gaussian matrix
Phi = orth(Phi')';      
A = Phi * Psi; % sensing matrix
sigma = 0.005;      
e = sigma*randn(M,1);    
y = Phi * x + e; %noisy observation vector y
% y = Phi * x; %observation vector y

%% restores the reconstructed signal x
tic        
theta = IHT_Basic(y,A,K);   
% theta = cs_iht(y,A,size(A,2));  
theta = hard_l0_Mterm(y,A,size(A,2),round(1.5*K),'verbose',true);  
xhat = Psi * theta;% x=Psi * theta        
toc        

norm(xhat-x)

%% drawing
figure; hold on;    
stem(xhat,'k.-');% plots the recovery signal of x
plot(x,'r');% plots the original signal x
legend('Recovery','Original')        





