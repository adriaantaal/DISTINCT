function [xhat,f2] = solvediagDXEv2(y,Dsmall,Esmall,xhat,maxIter,mu,gradconstdiag,ETE,DTD)
%solvediagDXE Solves system ||A*X*D' - y||^2_2 with:
    %Diagonal structure of X
    %Nonnegative entries in X 
    
    %The algorithm using proximal gradient descent 
    %and applies proximal projection operators 
    
    %https://math.stackexchange.com/questions/2421545/matrix-linear-least-squares-problem-with-diagonal-matrix-constraint
    %derivative without Estar applied, derivative  = D'D*X*B'B-D'yB 
    % gradconst = D'*y*E;
    % ETE = E'*E;
    % DTD = D'*D;

    % %derivative with Estar applied
    % gradconst = Estar*E*D'*y*Estar';
    % EsEDTD = Estar*E*DTD;
    % EEs = E'*Estar';

%parameters
if ~exist('maxIter','var')
    maxIter = 1e3; end

if ~exist('mu','var')
    mu = 5e-1; end

%prepare the matrices. 
if ~exist('gradconstdiag','var')
    gradconstdiag = diag(Dsmall'*y*Esmall); end

if ~exist('ETE','var')
    ETE = Esmall'*Esmall; end

if ~exist('DTD','var')
    DTD = Dsmall'*Dsmall; end

f2 = zeros(maxIter,1); 
f2(1) = Inf;
gradXold = 0;

for ii = 1:maxIter
    Xhat = diag(xhat);
    
%     %step 0: reporting
%     if (mod(ii,1000)==0)
%         disp(['Iteration ' num2str(ii)])
%     end
      
    %step 1: proximal gradient descent, find gradient. 
        %This is the slowest part of the algorithm!
        %be careful, DTDETEdiag.*xhat !=! diag(DTD*Xhat*ETE)

    %as this problem is dense, simply the direct calculation is the fastest method:
    gradXdiag = diag(DTD*Xhat*ETE)-gradconstdiag;
    DeltagradX = gradXdiag - gradXold;
        
    %step 2: apply gradient descent (negative), 
    %we are only interested in the diagonal entries anyways
%     Xhat = Xhat - mu*gradX;
    xhat = xhat- mu*gradXdiag; %diagonal version
    
    %step 3: soft threshold in population in xhat
    %omitted as its included in the in-crowd
%     Xhat = wthresh(Xhat,'s',lambda);
    
    %step 4:project onto the set of diagonal matrices
%     Xhat = diag(diag(Xhat));
    
    %step5: project onto nonnegative set
%     Xhat(Xhat<0) = 0;
    xhat(xhat<0) = 0;
    
    
    %step 6: check for convergence
%     f2(ii) = norm(y-Dsmall*Xhat*Esmall')^2; %this is very slow
    %check convergence in the gradient way because the gradient is much smaller
    if (ii > 1) 
        f2(ii) = (DeltagradX'*DeltagradX);
        if (abs(f2(ii)-f2(ii-1))/f2(ii) <= 1e-10)
%             disp(['Internal DXE solver converged after ' num2str(ii) ' iterations'])
            break
        end
    end
    gradXold = gradXdiag;
end
f2(ii+1:end) = [];
end