function [S,costheta] = verify_layout_separability(El,method,noise)
%VERIFY_LAYOUT_SEPARABILITY evaluates signal gain and separability

%X-axis: the absolute signal 
S = vecnorm(El)';

%apply signal weighing vector. Noiseless, this doesnt make a large difference
if nargin > 2
    x = ones(size(El,2),1);
%     Q = cov((El*x+noise)')^-0.5;
    Q = diag(El*x+noise);
    El = Q*El;
end

%Y-axis: the signal separability
%cos(thetai) = (||ai|| * ||wi||)^-1
if method == 1
    R = rank(El);
    %method 1: calculate angle between vector Ei and hyperplane Ej for every source
    %if Ej is full row rank this method doesnt work 
    %because then Ej has no direct normal subspace
    
    for ii = 1:size(El,2)
        disp(['ii = ' num2str(ii)]);
        Ej = El;
        Ej(:,ii) = [];

        %method 1a, find single vector normal to the subspace
%         Ej=Ej-mean(Ej);
%         [U,S,~]=svd(Ej,'econ');
%         normal=U(:,end);
%         costheta(ii) = abs( (El(:,ii)'*normal) ./ (norm(El(:,ii))*norm(normal)) );
        
        %method 1b: find whole subspace normal to Ej
        [U,~,~]=svd(Ej,'econ');
        normal=U(:,R+1:end);
        proj = normal - El(:,ii)*(El(:,ii)'*normal);
        %2 norm is largest singular value, which corresponds to largest
        %inner product between El and Ej
        costheta(ii) = cos(asin(min(1,norm(proj)))); 

        %method 1c, use directly built-in method
        %super slow!!
%         costheta(ii) = cos(subspace_short(El(:,ii),normal));
    end

elseif method == 2
    W = pinv(El);
    %method 2a : calculate inverse norm of W*E, this doesnt work
    costheta = 1./(vecnorm(W').*vecnorm(El))';
%     costheta = 1 ./ (diag(W*E)); %almost equivalent to above
    
elseif method == 3
    %this method also works for full row rank E
    W = pinv(El);
    %method 2b : calculate the proper angle between W and E
    costheta = diag(W*El) ./ (vecnorm(W').*vecnorm(El))';

elseif method == 4
    %Maybe this works?
    W = pinv(El);
    %method 2b : calculate the proper angle between W and E
    costheta = diag(W*El) ./ (vecnorm(W').*vecnorm(El))';
end

end