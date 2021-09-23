function X_est = nesterovGradDesc_AT_v2(y,A,At,lmbd2,maxItr,tol)
%% Author:

% Vivek Boominathan
% Rice University
% vivekb@rice.edu

%% Size
szX = A.szInp; %(1), A.szInp(2)];

%% Step size
[~,Em] = eigMax_AtA(A,At,100,lmbd2);
L = 2*Em;
invL = 1/L;

%% Initialization
% disp('Intializing..');

t_prv = 1;
X_prv = zeros(szX,'like',y);
V_prv = X_prv;

k = 1;

fidFn = @(X) sum(reshape(A*X - y,[],1).^2) + lmbd2*sum(X(:).^2);
gradFn = @(X) - 2*(At*y) + 2*(At*(A*X)) + 2*lmbd2*X;

%% Reconstruction
% disp('Reconstructing..');
fidFn_nxt = fidFn(X_prv);

tic,
while 1
    fidFn_prv = fidFn_nxt;

    X_nxt = V_prv - invL * gradFn(V_prv);

    % Neterov's
    t_nxt = 0.5*(1+sqrt(1+4*t_prv^2));
    V_nxt = X_nxt + (t_prv-1)*(X_nxt-X_prv)/t_nxt;

    fidFn_nxt = fidFn(X_nxt);
    if fidFn_nxt > fidFn_prv
        t_prv = 1;
    else
        t_prv = t_nxt;
    end

    X_prv = X_nxt;
    V_prv = V_nxt;

    k = k+1;

    if k>maxItr, break; end
    
    relErr = abs(fidFn_nxt-fidFn_prv)/fidFn_prv;
    if relErr < tol, break; end 
end

X_est = X_prv;
        
end



