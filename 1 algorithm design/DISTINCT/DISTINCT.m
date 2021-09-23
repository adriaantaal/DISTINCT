function [xhatFull,f1,f2,f1f2] = DISTINCT(y,D,Ek,lambda,ETE,DTD,gradconstdiag,intSolver,maxIter,maxIterInt,addSize,mu)
%DISTINCT Solves system lambda ||X||_1 + ||A*X*D' - y||^2_2 with:
    %Diagonal structure of X
    %Nonnegative entries in X 
    
    %The algorithm solves on an active set via the in-crowd 
    %and solves exactly on that in-crowd via exploiting diagonal structure
    
    %https://math.stackexchange.com/questions/2421545/matrix-linear-least-squares-problem-with-diagonal-matrix-constraint
    %derivative without Estar applied, derivative  = D'D*X*B'B-D'yB 
    % gradconst = D'*y*E;
    % ETE = E'*E;
    % DTD = D'*D;
    %
    
%parameters
if ~exist('mu','var')
    mu = 5e-1; end

if ~exist('lambda','var')
    lambda = max(diag(D'*y*Ek))/5; end

if ~exist('intSolver','var')
    intSolver = 0; end

if ~exist('addSize','var')
    addSize = 25; end

if ~exist('maxIter','var')
    maxIter = 1e2; end

if ~exist('maxIterInt','var')
    maxIterInt = 1e4; end

%prepare the matrices. 
if ~exist('gradconstdiag','var')
    gradconstdiag = diag(D'*y*Ek); end

if ~exist('ETE','var')
    ETE = Ek'*Ek; end

if ~exist('DTD','var')
    DTD = D'*D; end

%initialization
N = size(D,2);
xhatFull = zeros(N,1);
old_in_crowd = [];
in_crowd = [];
f1 = zeros(maxIter,1); 
f2 = zeros(maxIter,1); 
R = y;
yvec = y(:);


for ii = 1:maxIter

    %calculate usefullness, full matrix calculation (D'*R*E) is slow
        %accelerated 4x by calculating only diagonal entries
    u = zeros(N,1); DTR = D'*R; for nn = 1:N; u(nn,1) = DTR(nn,:)*Ek(:,nn); end
    
    %find the highest usefullness, all fast
    [us,indus] = sort(abs(u),'descend');
    us = wthresh(us,'s',lambda);
    us(addSize+1:end) = 0;  
    toAdd = indus(us>0);

    %add to in-crowd
    in_crowd = unique([in_crowd; toAdd],'stable');
    Dsmall = D(:,in_crowd);
    Esmall = Ek(:,in_crowd);
    

    %this part is the slowest of all
    %giving the precalculated matrices is 10% faster 
    %all-diagonalized method is 25% faster
    %checking for convergence via derivative is much faster too
    if intSolver == 0
        xhat = solvediagIntProd(yvec,Dsmall,Esmall);        
        Xhat = diag(xhat);
    elseif intSolver == 1 %the direct solution can be achieved by:       
        x0 = xhatFull(in_crowd);
        [xhat,~] = solvediagDXEv2(y,Dsmall,Esmall,x0,maxIterInt,mu,gradconstdiag(in_crowd),ETE(in_crowd,in_crowd),DTD(in_crowd,in_crowd));    
        Xhat = diag(xhat);
    elseif intSolver == 2
        %adding the BB is faster, but less precise, so it needs more in-crowd iterations
        x0 = xhatFull(in_crowd);
        [Xhat,~] = solvediagDXEBB(y,Dsmall,Esmall,X0,maxIterInt,mu);
        xhat = diag(Xhat);
    end

    %these operations are all fast
    %remove zero entries from in-crowd 
    in_crowd = in_crowd(xhat~=0);

    %re-create the full version such that the in-crowd can be added
    xhatFull(old_in_crowd) = 0;
    xhatFull(in_crowd) = xhat(xhat~=0);    

    %calculate the residual, very fast
    R = y-Dsmall*Xhat*Esmall';

    %report values
    f1(ii) = lambda*sum(xhat);
    f2(ii) = norm(R);
    f1f2(ii) = f1(ii)+f2(ii);
    
    %loop back
    old_in_crowd = in_crowd;
    
    if (ii > 1) 
        if (abs(f1f2(ii)-f1f2(ii-1))/f1f2(ii) <= 1e-10) && (f1(ii)~=0)
%             disp('Cost function reached tolerance')
%             disp(['DISSTINCT solution found after ' num2str(ii) ' iterations'])
            break
        end
    end
end

f1(ii+1:end) = [];
f2(ii+1:end) = [];
f1f2(ii+1:end) = [];

end

