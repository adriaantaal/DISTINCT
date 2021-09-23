function [xhatFull,DF,EF,f1,f2,f1f2] = DISTINCT_NNMTF(y,lambda,D0,E0,x0,intSolver,maxIter,maxIterInt,addSize)
%DISTINCT_NNMTF Non-negative matrix tri-factorization based on DISTINCT
    %Factorizes Y = DXE 
    %By minimizing cost function lambda ||X||_1 + ||D*X*E' - y||^2_2 with:
        %Diagonal constraint of X
        %Nonnegative entries in D,X,E 
    
    %The NNMTF is solved on an active set via the in-crowd 
    %and solves exactly on that in-crowd via exploiting diagonal structure

%prepare the matrices.
if ~exist('x0','var')
    x0 = zeros(size(D0,2),1); end

DF = D0;
EF = E0;
    
%parameters
if ~exist('lambda','var')
    lambda = max(diag(D0'*y*E0))/5; end

if ~exist('intSolver','var')
    intSolver = 0; end

if ~exist('addSize','var')
    addSize = 25; end

if ~exist('maxIter','var')
    maxIter = 1e2; end

if ~exist('maxIterInt','var')
    maxIterInt = 1e4; end

%initialization
xhatFull = x0;
old_in_crowd = [];
in_crowd = [];
f1 = zeros(maxIter,1); 
f2 = zeros(maxIter,1); 
R = y;

N = size(DF,2);

for ii = 1:maxIter
    disp(['ii= ' num2str(ii) ', Ic = ' num2str(numel(in_crowd))]);
    
    %calculate usefullness, full matrix calculation (D'*R*E) is slow
        %accelerated 4x by calculating only diagonal entries
    u = zeros(N,1); DTR = DF'*R; for nn = 1:N; u(nn,1) = DTR(nn,:)*EF(:,nn); end

    %find the highest usefullness, all fast
    [us,indus] = sort(abs(u),'descend');
    us = wthresh(us,'s',lambda);
    us(addSize+1:end) = 0;  
    toAdd = indus(us>0);

    %add to in-crowd
    in_crowd = unique([in_crowd; toAdd],'stable');
    DFsmall = DF(:,in_crowd);
    EFsmall = EF(:,in_crowd);
    
    %Initializing to zeros fails the calculation
%     x0 = xhatFull(in_crowd);
%     x0 = max(0,diag(pinv(DFsmall)*y*pinv(EFsmall'))); %change this!!!
    x0 = calc_diag_tri(pinv(DFsmall),y,pinv(EFsmall'),1:length(in_crowd)); 
    
    %call the internal solver, solve exactly on the smaller in-crowd
    if intSolver == 0 %this converges reliably
        [DFsmall,xhat,EFsmall,~] = solvediagNNMTFmult(y,DFsmall,x0,EFsmall,maxIterInt);
        Xhat = diag(xhat);
    elseif intSolver == 1 %the pseudoinverse method doesnt always converge
        [DFsmall,xhat,EFsmall,~] = solvediagNNMTFals(y,DFsmall,x0,EFsmall,maxIterInt);
        Xhat = diag(xhat);
    end

    %these operations are all fast
    %remove zero entries from in-crowd 
    in_crowd = in_crowd(xhat~=0);

    %re-create the full version such that the in-crowd can be added
    xhatFull(old_in_crowd) = 0;
    xhatFull(in_crowd) = xhat(xhat~=0);    

    %Same for EF and DF
    DF(:,old_in_crowd) = 0; EF(:,old_in_crowd) = 0;
    DF(:,in_crowd) = DFsmall(:,(xhat~=0));
    EF(:,in_crowd) = EFsmall(:,(xhat~=0));

    %calculate the residual, very fast on small in-crowd
    R = y-DFsmall*Xhat*EFsmall';

    %report values
    f1(ii) = lambda*sum(xhat);
    f2(ii) = norm(R);
    f1f2(ii) = f1(ii)+f2(ii);

    if (ii > 1) 
        if ((f1f2(ii-1)-f1f2(ii))/f1f2(ii) <= 1e-10)
            disp('Cost function reached tolerance')
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

