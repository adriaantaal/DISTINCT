function [xhat,f2] = solvediagDXEv3(y,Dsmall,Esmall,x0,maxIter,mu,gradconstsmall,DTDdotETE)
%solvediagDXE Solves system ||A*X*D' - y||^2_2 with:
    %Diagonal structure of X
    %Nonnegative entries in X 
    
    %The algorithm using proximal gradient descent 
    %and applies proximal projection operators 
    
    %https://math.stackexchange.com/questions/2421545/matrix-linear-least-squares-problem-with-diagonal-matrix-constraint
    %derivative without Estar applied, derivative  = diag(D'D*X*B'B-D'yB)
    % gradconst = D'*y*E;


%parameters
if ~exist('maxIter','var')
    maxIter = 1e3; end

if ~exist('mu','var')
    mu = 1e-1; end

%prepare the matrices. 
if ~exist('gradconstdiag','var')
    gradconstsmall = diag(Dsmall'*y*Esmall); end

if ~exist('DTDdotETE','var')
    DTDdotETE = (Dsmall'*Dsmall).*(Esmall'*Esmall); end

f2 = zeros(maxIter,1); 
f2(1) = Inf;
% gradXold = 0;
xhat = x0;

for ii = 1:maxIter
    
%     %step 0: reporting
%     if (mod(ii,1000)==0)
%         disp(['Iteration ' num2str(ii)])
%     end
      
    %step 1: proximal gradient descent, find gradient. 
    %This is the slowest part of the algorithm!
    %diag(DTD*Xhat*ETE) = DTDdotETE'*xtrue; -> precalc DTDdotETE
    gradX = DTDdotETE*xhat-gradconstsmall;
%     deltagradX = gradX - gradXold;
        
    %apply gradient descent (negative), 
    %we are only interested in the diagonal entries anyways
    xhat = xhat - mu*gradX; %diagonal version
    
    %step 2: soft threshold in population in xhat
    %omitted as its included in the in-crowd
%     Xhat = wthresh(Xhat,'s',lambda);
    
    %step 3:project onto the set of diagonal matrices
%     Xhat = diag(diag(Xhat));
    
    %step4: project onto nonnegative set
%     xhat(xhat<0) = 0;
    
    %this is very time consuming and takes most of the algorithm time
%     f2(ii) = norm(y-Dsmall*diag(xhat)*Esmall','fro')^2; %this is still faster than the svds method!

    %step 5: check for convergence in the gradient way because the gradient is much smaller
    if (ii > 1) 
%         f2(ii) = (deltagradX'*deltagradX);
%         f2(ii) = norm(y-Dsmall*diag(xhat)*Esmall','fro');
        f2(ii) = gradX'*gradX;
%         if (abs(f2(ii)-f2(ii-1))/f2(ii) <= 1e-10)
        if (f2(ii) <= 1e-8)
%             disp(['Internal DXE solver converged after ' num2str(ii) ' iterations'])
            break
        end
    end
%     gradXold = gradX;
end
f2(ii+1:end) = [];
end