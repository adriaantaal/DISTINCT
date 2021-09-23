function [S,costheta] = verify_layout_separability2D(D,Ek)
%VERIFY_LAYOUT_SEPARABILITY evaluates signal gain and separability

%X-axis: the absolute signal 
S = vecnorm(D)'.*vecnorm(Ek)';

%Maybe this works for 2D?
W = pinv(D);
P = pinv(Ek');

%calculate the proper angle between W and E
% costheta = diag(W*D) ./ (vecnorm(W').*vecnorm(D))';
costheta = diag(W*D*Ek'*P) ./ (vecnorm(W').*vecnorm(D).*vecnorm(Ek).*vecnorm(P))';

end