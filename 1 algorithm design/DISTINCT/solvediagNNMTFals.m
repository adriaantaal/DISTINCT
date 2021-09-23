function [DFsmall,xhat,EFsmall,f2] = solvediagNNMTFals(y,DFsmall,x0,EFsmall,maxIterInt)
%NNMTF internal solver
%   using alternate least squares (bunch of pseudoinverses)
%   Y = DF * X * EF
    
%parameters
if ~exist('maxIter','var')
    maxIterInt = 1e3; end

f2 = zeros(maxIterInt,1); 
f2(1) = Inf;

Xhat = diag(x0);

for ii = 1:maxIterInt
    
%     %step 0: reporting
%     if (mod(ii,1000)==0)
%         disp(['Iteration ' num2str(ii)])
%     end
      
    %step 1: Update D
    DFsmall = max(0,y*pinv(EFsmall')*pinv(Xhat));
    
    %step 2: Update E
    EFsmall = max(0,pinv(Xhat)*pinv(DFsmall)*y)';
    
    %step 3: Update X
    Xhat = max(0,pinv(DFsmall)*y*pinv(EFsmall'));

    
    %step 4: soft threshold in population in xhat
    %omitted as its included in the in-crowd
%     Xhat = wthresh(Xhat,'s',lambda);
    
    %step 5:project onto the set of diagonal matrices
    %omitted as the update automatically preserves diagonal
%     Xhat = diag(diag(Xhat));
    xhat = diag(Xhat);
    
    %step6: project onto nonnegative set
    %omitted as automcatic for NNMTF
   
    %step 7: check for convergence in norm of DeltaX
    if (ii > 1) 
        DeltaX = xhat-xhatold;
        f2(ii) = (DeltaX'*DeltaX);
        if (abs(f2(ii)-f2(ii-1))/f2(ii) <= 1e-8)
            disp(['Internal NNMTF solver converged after ' num2str(ii) ' iterations'])
            break
        end
    end
    xhatold = xhat;
end
f2(ii+1:end) = [];

end

