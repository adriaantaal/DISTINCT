function [xhatFull,f1,f2,f1f2] = DISTINCTv2(y,D,E,lambda,H,gradconstdiag,intSolver,maxIter,maxIterInt,addSize)
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
if ~exist('lambda','var')
    lambda = max(diag(D'*y*E))/10; end

if ~exist('intSolver','var')
    intSolver = 0; end

if ~exist('addSize','var')
    addSize = 5; end

if ~exist('maxIter','var')
    maxIter = 1e2; end

if ~exist('maxIterInt','var')
    maxIterInt = 1e4; end

if ~exist('gradconstdiag','var')
    gradconstdiag = diag(D'*y*E); end

if ~exist('H','var')
    H = (D'*D).*(E'*E); end

%initialization
xhatFull = zeros(size(D,2),1);
old_in_crowd = [];
in_crowd = [];
f1 = zeros(maxIter,1); 
f2 = zeros(maxIter,1); 
f1f2 = zeros(maxIter,1);
R = y;
yvec = y(:);

N = size(D,2);
options = optimset('MaxIter',maxIterInt,'Display','off','TolFun',1e-6,'TolX',1e-6);
opts.SYM = true;
opts.POSDEF = true;

for ii = 1:maxIter

    u = zeros(N,1);
    %calculate usefullness, full matrix calculation is slow
        %todo: accelerated 4x by calc only diagonal entries
        %kept the smaller matrix multiplication intact
    DTR = D'*R;
    for nn = 1:N
        u(nn,1) = DTR(nn,:)*E(:,nn);
    end

    %find the highest usefullness, all fast
    [us,indus] = sort(abs(u),'descend');
    us = wthresh(us,'s',lambda);
    us(addSize+1:end) = 0;  
    toAdd = indus(us>0);

    %add to in-crowd
    in_crowd = unique([in_crowd; toAdd],'stable');
    Dsmall = D(:,in_crowd);
    Esmall = E(:,in_crowd);

    %this part is the slowest of all
    %giving the precalculated matrices is 10% faster 
    %all-diagonalized method is 25% faster
    %checking for derivative convergence is much faster too
    if intSolver == 0
        xhat = solvediagIntProd(yvec,Dsmall,Esmall);        
    elseif intSolver == 1 %the direct solution can be achieved by:   
        [xhat,~] = lsqr(H(in_crowd,in_crowd),gradconstdiag(in_crowd),[],maxIterInt);
%         [xhat,~] = solvediagDXEv3(y,Dsmall,Esmall,x0,maxIterInt,mu,gradconstdiag(in_crowd),H(in_crowd,in_crowd));
    elseif intSolver == 2
        xhat = H(in_crowd,in_crowd)\gradconstdiag(in_crowd);
    elseif intSolver == 3
        xhat = linsolve(H(in_crowd,in_crowd),gradconstdiag(in_crowd),opts);
%         [xhat,~,~,~] = quadprog(H(in_crowd,in_crowd),gradconstdiag(in_crowd),[],[],[],[],[],[],[],options);
    end
%     xhat = xhat(~isnan(xhat));
    
    %these operations are all fast
    %remove zero entries from in-crowd 
    in_crowd = in_crowd(xhat~=0); 

    %re-create the full version such that the in-crowd can be added
    xhatFull(old_in_crowd) = 0;
    xhatFull(in_crowd) = xhat(xhat~=0);    

    %calculate the residual, very fast
    R = y-(Dsmall.*xhat')*Esmall';
    
    %report values, very fast using Frob norm
    f1(ii) = lambda*norm(xhat,1);
    f2(ii) = norm(R,'fro').^2;
    f1f2(ii) = f1(ii)+f2(ii);

    if (ii > 1) 
        if (abs(f1f2(ii)-f1f2(ii-1))/f1f2(ii) <= 1e-6)
%             disp('Cost function reached tolerance')
%             disp(['DISSTINCT solution found after ' num2str(ii) ' iterations'])
            break
        end
    end

    %loop back
    old_in_crowd = in_crowd;
end

f1(ii+1:end) = [];
f2(ii+1:end) = [];
f1f2(ii+1:end) = [];

end