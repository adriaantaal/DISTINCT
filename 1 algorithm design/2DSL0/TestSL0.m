% Using function SL0:

% Goal: Finding sparse solution of the system x=As, where x is n by 1, s is
% m by 1, and m>n.

% Parameters
m = 1000;   % No of sources
n = 400;   % No of sensors
q = 3; % Maximum number of iterations for the internal steepest descent loop
sigma_off = 0.001;

% Important Algorithm's Parameters
sigma_decrease_factor = 0.5;
q = 3;
if sigma_off>0
    sigma_min = sigma_off*4;
else
    sigma_min = 0.00001;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choose Real or Complex case by uncommenting the appripriate lines: %%

% Complex:
% true_s = sparseComplexSigGen4plusNoise(m,floor(n/4),sigma_off); 
% A = randn(n,m) + sqrt(-1)*randn(n,m);
 
% or Real:
true_s = sparseRealSigGen4plusNoise(m,floor(n/4),sigma_off);  
A = randn(n,m);


for i=1:m
    A(:,i) = A(:,i) / norm(A(:,i));
end

%x = A*s + 0.01*randn(n,1);
x = A*true_s;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Global Init:
A_pinv = pinv(A); % Time for computing A_pinv not is counted, because it is done only once.
mu_0 = 2;

fprintf('computing...\n')
tic 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Different forms of calling SL0 (uncomment the desired action): %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% To provide minimum number of input arguments:
%s=SL0(A, x, sigma_min);

% To see the progress:
s=SL0(A, x, sigma_min, sigma_decrease_factor, mu_0, q, A_pinv, true_s);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of calling SL0                                             %  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ElapsedTime=toc;
fprintf('Done. Elapsed Time=%f,  SNR=%f dB\n', ElapsedTime, estimate_SNR(s,true_s))



