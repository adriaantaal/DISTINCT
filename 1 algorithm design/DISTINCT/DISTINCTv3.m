function [xhatFull,f1,f2,f1f2] = DISTINCTv3(y,D,E,lambda,H,gradconstdiag,U,s,V,intSolver,maxIter,maxIterInt,addSize,max_contrast_ratio)
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

%Changelog
    %v3: calculated usefullness through the Hessian, speedup of
    %v2: added linear solving methods of the normal equations
    
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

if ~exist('max_contrast_ratio','var')
    max_contrast_ratio = Inf; end

if ~exist('U','var')
    [U,S,V] = svd(y,'econ'); 
    U = U(:,1:L); V = V(:,1:L); s = diag(S(1:L,1:L));
end

%initialization
xhatFull = zeros(size(D,2),1);
xhat = [];
old_in_crowd = [];
in_crowd = [];
f1 = zeros(maxIter,1); 
f2 = zeros(maxIter,1); 
f1f2 = zeros(maxIter,1);
yvec = y(:);

% N = size(D,2);
% options = optimset('MaxIter',maxIterInt,'Display','off','TolFun',1e-6,'TolX',1e-6);
opts.SYM = true;
opts.POSDEF = true;

for ii = 1:maxIter

    %calculate usefullness, full matrix calculation is slow
        %accelerated further through using the Hessian
    if isempty(in_crowd)
        u = gradconstdiag;
    else
        u = gradconstdiag - H(:,in_crowd)*xhat;
    end
        
    %find usefullness bigger than lambda
    [us,indus] = sort(abs(u),'descend');
%     us = wthresh(us,'s',lambda);
    us(addSize+1:end) = 0;  
    toAdd = indus(us>lambda);
    
    %add to in-crowd
    in_crowd = unique([in_crowd; toAdd],'stable');
    Dsmall = D(:,in_crowd);
    Esmall = E(:,in_crowd);

    %Solve F-norm term exactly on the subset
    if intSolver == 0
        xhat = solvediagIntProdv2(U,s,V,Dsmall,Esmall); 
    elseif intSolver == 1    
        x0 = zeros(length(in_crowd),1); x0(1:length(old_in_crowd)) = xhat;
        [xhat,~] = lsqr(H(in_crowd,in_crowd),gradconstdiag(in_crowd),[],maxIterInt,[],[],x0);
    elseif intSolver == 2 %the direct solution through long inversion
        xhat = H(in_crowd,in_crowd)\gradconstdiag(in_crowd);
    elseif intSolver == 3
        xhat = linsolve(H(in_crowd,in_crowd),gradconstdiag(in_crowd),opts);
%         [xhat,~,~,~] = quadprog(H(in_crowd,in_crowd),gradconstdiag(in_crowd),[],[],[],[],[],[],[],options);
    end
    
    %remove indices smaller than max_contrast_ratio
    xhat(abs(xhat)<(max(xhat)./max_contrast_ratio)) = 0;
    
    %these operations are all fast
    %remove zero entries from in-crowd 
    in_crowd = in_crowd(xhat~=0); 

    %re-create the full version such that the in-crowd can be added
    xhat = xhat(xhat~=0);
    xhatFull(old_in_crowd) = 0;
    xhatFull(in_crowd) = xhat;    
    
    %report values, very fast using Frob norm
    f1(ii) = lambda*norm(xhat,1);
    f2(ii) = xhat'*H(in_crowd,in_crowd)*xhat;
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