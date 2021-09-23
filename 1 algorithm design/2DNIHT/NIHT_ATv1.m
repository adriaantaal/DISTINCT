function [x, error]=NIHT_ATv1(y,D,L)
% NIHT_ATv1: Hard thresholding algorithm that keeps exactly L elements 
% in each iteration. 
%
% This algorithm has certain performance guarantees as described in [1],
% [2] and [3].
%
% For 1D case
%     gradient  = D'Dx-D'y
%     error     = ||Dx-y||22 = 2-norm of residue 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Usage
%
%   [x, error]=NIHT_ATv1(y,D,L)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Input
%
%   Mandatory:
%               x   Observation vector to be decomposed
%               D   An NxM matrix (m must be dimension of x)
%               E   An MxK matrix (m must be dimension of x)
%               L   non-zero elements to keep in each iteration
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Outputs
%
%    s              Solution vector 
%    err_mse        Vector containing mse of approximation error for each 
%                   iteration
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Description
%
%   Implements the M-sparse algorithm described in [1], [2] and [3].
%   This algorithm takes a gradient step and then thresholds to only retain
%   M non-zero elements. It allows the step-size to be calculated
%   automatically as described in [3] and is therefore now independent from 
%   a rescaling of P.
%   
%   
% References
%   [1]  T. Blumensath and M.E. Davies, "Iterative Thresholding for Sparse 
%        Approximations", submitted, 2007
%   [2]  T. Blumensath and M. Davies; "Iterative Hard Thresholding for 
%        Compressed Sensing" to appear Applied and Computational Harmonic 
%        Analysis 
%   [3] T. Blumensath and M. Davies; "A modified Iterative Hard 
%        Thresholding algorithm with guaranteed performance and stability" 
%        in preparation (title may change) 
% See Also
%   hard_l0_reg
%
% Copyright (c) 2007 Thomas Blumensath
%
% The University of Edinburgh
% Email: thomas.blumensath@ed.ac.uk
% Comments and bug reports welcome
%
% This file is part of sparsity Version 0.4
% Created: April 2007
% Modified January 2009
%
% Part of this toolbox was developed with the support of EPSRC Grant
% D000246/1
%
% Please read COPYRIGHT.m for terms and conditions.

[n,m]       = size(D);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Default parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
STOPTOL     = 1e-16;
MAXITER     = n^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

error       = zeros(MAXITER,0);
x           = zeros(m,1);
Residual    = y;
init_error  = Residual'*Residual/n;
Dx          = D*x;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        Main algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
done = 0;
iter=1;

while ~done
    %Calculate optimal step size and do line search
    oldx                =   x;
    oldDx               =   Dx;
    IND                 =   x~=0;
    grad                =   D'*Residual;
    % If the current vector is zero, we take the largest elements in grad
    if sum(IND)==0
        [dsort sortdind]    =   sort(abs(grad),'descend');
        IND(sortdind(1:L))  =   1;    
     end  

    igrad               =   (IND.*grad);
    Dgrad               =   D*igrad;
    mu                  =   igrad'*igrad/(Dgrad'*Dgrad);
    x                   =   oldx + mu * grad;
    [~, sortind]     =   sort(abs(x),'descend');
    x(sortind(L+1:end)) =   0;
    Dx                  =   D*x;

    % Calculate step-size requirement 
    omega               =   (norm(x-oldx)/norm(Dx-oldDx))^2;

    % As long as the support changes and mu > omega, we decrease mu
    while mu > (0.99)*omega && sum(xor(IND,x~=0))~=0 && sum(IND)~=0
%             display(['decreasing mu'])

                % We use a simple line search, halving mu in each step
                mu                  =   mu/2;
                x                   =   oldx + mu * grad;
                [~, sortind]     =   sort(abs(x),'descend');
                x(sortind(L+1:end)) =   0;
                Dx                  =   D*x;
                % Calculate step-size requirement 
                omega               =   (norm(x-oldx)/norm(Dx-oldDx))^2;
    end
            
     Residual            =   y-Dx;
     error(iter) = Residual'*Residual/n;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        Are we done yet?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
     if iter>1 
         if ((error(iter-1)-error(iter))/init_error<STOPTOL) || (error(iter)<STOPTOL)
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
