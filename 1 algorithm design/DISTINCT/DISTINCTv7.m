function [xhatFull,J] = DISTINCTv7(yf2,D,E,lambda,H,gradconstdiag,U,s,V,intSolver,maxIter,maxIterInt,addSize)
%DISTINCT Solves system lambda ||X||_1 + ||A*X*D' - y||^2_2 with:
    %Diagonal structure of X
    %Nonnegative entries in X 
    
    %The algorithm solves on an active set
    %and solves exactly on that active set exploiting diagonal structure

%Changelog
    %v7: Sped up convergence checking. Removes zero entries from H and g
    %v6: Required subtraction of lambda outside the function, us > 0
    %v5: Improved speedup of incrowd and newcrowd calculation, 
    %    removed xhatFull, removed max_contrast_ratio
    %v4: added LSQR solver with prior information 
    %v3: calculated usefullness through the Hessian, speedup of
    %v2: added linear solving methods of the normal equations

%initialization
xhatFull = zeros(size(D,2),1);
xhat = [];
old_in_crowd = [];
in_crowd = [];
J = zeros(maxIter,1); 
u = gradconstdiag;

for ii = 1:maxIter       
        %find usefullness bigger than lambda
        [us,indus] = sort(u,'descend');
        us(addSize+1:end) = 0;  
        toAdd = indus(us>0);

        %detect new crowd
        newcrowd = toAdd(~ismember(toAdd,in_crowd));

        %add new crowd to in-crowd
        in_crowd = [in_crowd;newcrowd];

        %Solve F-norm term exactly on the subset
        if intSolver == 0
            xhat = solvediagIntProdv2(U,s,V,D(:,in_crowd),E(:,in_crowd)); 

        elseif intSolver == 1    
            x0 = zeros(size(in_crowd)); 
            %old in-crowd are previous values
            x0(1:length(old_in_crowd)) = xhat; 
            %new in-crowd are values derived from the normal equations
            [xhat,~] = lsqr(H(in_crowd,in_crowd),gradconstdiag(in_crowd),[],maxIterInt,[],[],x0); 

        elseif intSolver == 2 %the direct solution through long inversion
            xhat = H(in_crowd,in_crowd)\gradconstdiag(in_crowd);

        elseif intSolver == 3
            xhat = linsolve(H(in_crowd,in_crowd),gradconstdiag(in_crowd),opts);
        end
        
        %remove zero entries from g and H to prevent reconsideration
%         if ~isempty(in_crowd(xhat==0))
%             remcrowd = in_crowd(xhat==0);
%             gradconstdiag(remcrowd) = [];
%             H(remcrowd,remcrowd) = [];
% 
%             %remove zero entries from in-crowd 
%             in_crowd = in_crowd(xhat~=0); 
% 
%             %remove from xhat
%             xhat = xhat(xhat~=0);
%         end
        
        in_crowd = in_crowd(xhat~=0); 
        xhat = xhat(xhat~=0);
        
        %new usefulness, also needed for cost function calculation
        u = gradconstdiag - H(:,in_crowd)*xhat;
        
        %check if cost function converged
        J(ii) = yf2 - xhat'*(u(in_crowd) + gradconstdiag(in_crowd) + lambda); %equivalent and much faster
        
        if ii>1
            if (abs(J(ii)-J(ii-1))/J(ii-1) < 1e-3) || (abs(J(ii)/yf2) < 1e-4)
                break;
            end
        end
        
        %loop back
        old_in_crowd = in_crowd;
end

xhatFull(in_crowd) = xhat;    
J(ii+1:end) = [];
end