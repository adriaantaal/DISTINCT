function [DFsmall,xhat,EFsmall,f2] = solvediagNNMTFmult(y,DFsmall,x0,EFsmall,maxIterInt)
%NNMTF internal solver
%   https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5746986/pdf/13040_2017_Article_160.pdf
%
    
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
    numer = (y*EFsmall*Xhat);
    DFsmall = max(0,DFsmall.*numer./(DFsmall*Xhat*(EFsmall'*EFsmall)*Xhat' + eps(numer)));
    
    %step 2: Update E
    numer = (y'*DFsmall*Xhat);
    EFsmall = max(0,EFsmall.*numer./(EFsmall*Xhat*(DFsmall'*DFsmall)*Xhat + eps(numer)));
    
    %step 3: Update X
    numer = (DFsmall'*y*EFsmall);
    Xhat = max(0,Xhat.*numer./((DFsmall'*DFsmall)*Xhat*(EFsmall'*EFsmall) + eps(numer)));
    
    %step 4: soft threshold in population in xhat
    %omitted as its included in the in-crowd
%     Xhat = wthresh(Xhat,'s',lambda);
    
    %step 5:project onto the set of diagonal matrices
    %omitted as the update automatically preserves diagonal
%     Xhat = diag(diag(Xhat));
     xhat = diag(Xhat);
    
    %step6: project onto nonnegative set
    %omitted as automcatic for NNMTF
   

    %step 7: check for convergence in norm of DeltaX. Takes ~25% of time
    if (ii > 1) 
        DeltaX = xhat-xhatold;
        if ~isempty(DeltaX) 
            f2(ii) = (DeltaX'*DeltaX);            
        else
            f2(ii) = 0;
        end
        if (abs(f2(ii)-f2(ii-1))/f2(ii) <= 1e-8)
            disp(['Internal NNMTF solver converged after ' num2str(ii) ' iterations'])
            break
        end
    end

    xhatold = xhat;
end
f2(ii+1:end) = [];

end

