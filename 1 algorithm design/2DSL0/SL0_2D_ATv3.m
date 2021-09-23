function x=SL0_2D_ATv3(y, sigma_min, sigma_decrease_factor, mu_0, q, Dpinv, Epinv, projD, projEk, sl0gradconstdiag, max_contrast_ratio, xtrue)
%
% SL0_2D_AT(D, E, y, sigma_min, sigma_decrease_factor, mu_0, q, Apinv, Epinv, xtrue)
%
%   Returns the sparsest vector s which satisfies underdetermined system of
%   linear equations  D*x*E=y, using  Smoothed L0  (SL0) algorithm.
%
%   Sequence of Sigma (sigma_min and sigma_decrease_factor):
%     This is a decreasing geometric sequence of positive numbers:
%       - The  first  element   of  the  sequence  of  sigma is  calculated 
%     automatically. The last  element  is  given  by  'sigma_min', and the 
%     change factor for decreasing sigma is given by 'sigma_decrease_factor'. 
%       - The default value of 'sigma_decrease_factor' is 0.5. Larger value 
%     gives better results  for less sparse sources, but it uses more steps 
%     on   sigma   to  reach  sigma_min,  and  hence  it  requires   higher 
%     computational cost.
%       - There is no default  value for  'sigma_min',  and  it  should  be 
%     provided  by  the  user (depending  on his/her estimated source noise 
%     level,  or  his/her  desired  accuracy).  By `noise' we mean here the
%     noise in the sources, that is, the energy of the inactive elements of
%     'x'.   For example,  by  the  noiseless  case,  we  mean the inactive
%     elements of 'x' are exactly equal to zero. As a rule of tumb, for the
%     noisy case,  sigma_min should be about 2 to 4  times  of the standard
%     deviation of this noise.  For the noiseless case, smaller 'sigma_min'
%     results in  better estimation of the sparsest solution, and hence its
%     value is determined by the desired accuracy.
% 
%   mu_0: 
%        The  value  of  mu_0  scales  the sequence of mu. For each vlue of 
%     sigma, the value of  mu is chosen via mu=mu_0*sigma^2. Note that this 
%     value effects Convergence.
%        The default value is mu_0=2 (see the paper).
%
%   q: 
%        number  of  iterations of the internal (steepest ascent) loop. The
%     default value is q=3.
%
%   Apinv & Epinv: 
%        is the  pseudo-inverse of matrix A defined by A_pinv=A'*inv(A*A'). 
%
%   xtrue: 
%        is the  true value of the  sparse  solution.  This argument is for
%     simulation purposes. If it is provided by the user, then the function
%     will  calculate the SNR of the estimation for each value of sigma and
%     it provides a progress report.
%


% Web-page:
% ------------------
%    http://ee.sharif.ir/~SLzero
%
% Code History:
%--------------
% Version 3.1: January 2021
%        For our case, X = diagonal
%        Therefore, we calculate x instead of calculating full matrix X
%       
%        x -> x - projD*X*col(projEk) + sl0gradconstdiag

% Version 3.0: January 2021
%        Adapted by Adriaan Taal to support 2D
%        Following the paper SPARSE DECOMPOSITION OF TWO DIMENSIONAL SIGNALS
%        By Aboozar Ghaffari, Massoud Babaie-Zadeh, Christian Jutten
%        
%        Added input of E and Epinv
%        X is a diagonal matrix NxN
%        x is the vector version of the diagonal of X
%        Removed input of E and D, they are only used in pinv and projection forms
%       
%       The projection onto the feasible set is now:
%       
%       X -> X - Dpinv*(D*X*E' - y)*Epinv' = X - projD*X*projE' + Dpinv*y*Epinv'
%       
%       Constant Dpinv*y*Epinv' is calculated ahead as "sl0gradconst"



% Version 2.0: 4 April 2010
%        Doing a few small modifications that enable the code to work also
%        for complex numbers (not only for real numbers).
%
% Version 1.3: 13 Sep 2008
%        Just a few modifications in the comments
%
% Version 1.2: Adding some more comments in the help section
%
% Version 1.1: 4 August 2008
%    - Using MATLAB's pseudo inverse function to generalize for the case
%      the matrix A is not full-rank.
%
% Version 1.0 (first official version): 4 July 2008.
%
% First non-official version and algorithm development: Summer 2006

if nargin < 12
    ShowProgress = false;
elseif nargin == 12
    ShowProgress = true;
else
    error('Error in calling SL0 function');
end

% Initialization
x = diag(Dpinv*y*Epinv');
sigma = 2*max(abs(x));

% Main Loop
while sigma>sigma_min
    disp('Decreasing sigma')
    for i=1:q
        %this inner loop is incredibly slow (60s for y = 1024x1024)
        delta = OurDelta(x,sigma); %very fast
        x = x - mu_0*delta;              %SL0 gradient descent
        
        %toss out elements smaller than max_constrast_ratio
        x = wthresh(x,'h',abs(max(x)/max_contrast_ratio));
        activeind = find(abs(x)>0);        
        
        %accelerated 2~3x by calculating only diagonal entries, verified for correct dimensions
        %accelerated another 20x by calculating only the active indices per max_contrast_ratio
        N = length(activeind);
        pDX = projD(activeind,activeind)*diag(x(activeind));
        u = zeros(N,1); for nn = 1:N; u(nn,1) = pDX(nn,:)*projEk(activeind,activeind(nn)); end
        %update x = x - u + sl0gradconstdiag
        x(activeind) = x(activeind) - u + sl0gradconstdiag(activeind);
    end
    
    if ShowProgress
        fprintf('     sigma=%f, SNR=%f\n',sigma,estimate_SNR(x,xtrue))
    end
    
    sigma = sigma * sigma_decrease_factor;
end
