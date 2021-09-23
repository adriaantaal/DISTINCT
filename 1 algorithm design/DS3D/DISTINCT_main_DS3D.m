%% Welcome to the code that simulates DISTINCT 
%First create a model for Detector and Emitter array
%Then apply the equations to verify reconstruction


%% STEP 0: import all functions and folders
Folder = fullfile(cd, '..');
addpath(genpath(fullfile(Folder)));

%add the whole folder on MAC here:
addpath(genpath('D:\Dropbox\96 Paper DISTINCT\2 algorithm verification'))

addpath(genpath('/proj2/adriaantaal/96 Paper DISTINCT'))

%% STEP 1 dick around with tensors
N = 10;
M = 8;
K = 3;
B = 5;

% X = tensor(zeros(N,N,N));
% X(2,2,2) = 1; X(7,7,7) = 1;
% Y = tensor(zeros(M,K,B));
% D = tensor(randn(M,N));
% E = tensor(randn(N,K));
% T = tensor(randn(N,B));

X = zeros(N,N,N);
X(2,2,2) = 1; X(7,7,7) = 1;
% X = zeros(N,1); X([2,7]) = 1;
D = (randn(M,N)); DTD = D'*D;
E = (randn(N,K)); ETE = E*E';
Q = (randn(N,B)); QTQ = Q*Q';

DX = pagemtimes(D,X);
DXE = pagemtimes(DX,E);
DXE = permute(DXE,[1 3 2]);
DXQE = pagemtimes(DXE,Q);
Y = DXQE;
Yt = tensor(Y);
T = hosvd(Yt,2*sqrt(eps));
%estimate usefullness

DTDXQE = pagemtimes(D',DXQE);
DTDXQTQE = pagemtimes(DTDXQE,Q');
DTDXQTQE = permute(DTDXQTQE,[1 3 2]);
DTDXQTQETE = pagemtimes(DTDXQTQE,E');

yvec = Y(:);
for ii = 1:N
    u(ii) = DTDXQTQETE(ii,ii,ii);
    gradconstdiag(ii) = DTDXQTQETE(ii,ii,ii);
    wii = tovec(tovec(D(:,ii)*E(ii,:))*Q(ii,:));
    x(ii) = yvec'*wii/(wii'*wii);
end


H = DTD.*ETE.*QTQ;
xhat = H\gradconstdiag'

