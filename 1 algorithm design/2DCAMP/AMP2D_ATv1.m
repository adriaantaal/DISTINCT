function [x, error]=AMP2D_ATv1(y,D,E,lambda,alpha)
% AMP2D_ATv1: 2D Approximate Message passage algorithms
%
% Approximate Message Passing Algorithm for Complex Separable Compressed Imaging
% Akira Hirabayashi, Jumpei Sugimoto, and Kazushi Mimura
% 
% Also based on:
% Message Passing Algorithms for Compressed Sensing
% David L. Donoho, Arian Maleki, Andrea Montanari
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Input
%
%   Mandatory:
%               y   Data vector to be decomposed
%               D   An nxm matrix (m must be dimension of x)
%               E   An mxk matrix (m must be dimension of x)
%           alpha   Step size
%          lambda   Soft thresholding operator
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Outputs
%
%    x              Solution vector 
%    error          MSE of approximation for each iteration
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Description
% 
% 
% 

[n,m]       = size(D);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Default parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
STOPTOL     = 1e-16;
MAXITER     = 1000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

error       = zeros(MAXITER,0);
X           = zeros(m,m);
oldX        = X;
Z           = y;
oldZ        = Z;
DZE         = D'*Z*E;
init_error  = norm(y);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        Main algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
done = 0;
iter=1;

while ~done
    disp(['Iteration ' num2str(iter)]);

    %oldold means n-1.
    %old means n
    %Updated values are n+1

    oldoldX             =   oldX;
    oldoldZ             =   oldZ;
    oldX                =   X;
        
    %step 1 = update X: soft threshold D'ZE + x
    X = wthresh(D'*oldZ*E+oldX,'s',lambda);
    
    %step 2 = update Z: use X and derivative of oldX
    gradX = mean(mean(gradient(wthresh(D'*oldoldZ*E+oldoldX,'s',lambda))));
    oldZ = y-D*oldX*E' - 1/alpha*oldoldZ*gradX;
    
    %step 3: Calculate error of solution  
    Residual = y-D*X*E';
    error(iter) = norm(Residual); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        Are we done yet?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
     if iter>1 
         if (abs(error(iter-1)-error(iter))/init_error<STOPTOL) || (error(iter)<STOPTOL)
             done=1;
              disp('Stopping. Exact signal representation found!')
         end
     end

     if iter >= MAXITER
         disp('Stopping. Maximum number of iterations reached!')
         done = 1; 
     end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    If not done, take another round
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     if ~done
        iter=iter+1;  
     end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  Only return as many elements as iterations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargout >=2
    error = error(1:iter);
end





function y = wthresh(x,sorh,t)
%WTHRESH Perform soft or hard thresholding. 
%   Y = WTHRESH(X,SORH,T) returns soft (if SORH = 's')
%   or hard (if SORH = 'h') T-thresholding  of the input 
%   vector or matrix X. T is the threshold value.
%
%   Y = WTHRESH(X,'s',T) returns Y = SIGN(X).(|X|-T)+, soft 
%   thresholding is shrinkage.
%
%   Y = WTHRESH(X,'h',T) returns Y = X.1_(|X|>T), hard
%   thresholding is cruder.
%
%   See also WDEN, WDENCMP, WPDENCMP.

%   M. Misiti, Y. Misiti, G. Oppenheim, J.M. Poggi 12-Mar-96.
%   Last Revision: 24-Jul-2007.
%   Copyright 1995-2007 The MathWorks, Inc.

switch sorh
  case 's'
    tmp = (abs(x)-t);
    tmp = (tmp+abs(tmp))/2;
    y   = sign(x).*tmp;
 
  case 'h'
    y   = x.*(abs(x)>t);
    
  otherwise
    y   = x.*(abs(x)>t);
end
