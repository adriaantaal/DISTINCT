%% Author:
% Adriaan Taal
% Goal of the code
    % 0. The separability is already incorporated in the structure of the model
        %Therefore, there is no further separability to exploit
    % 1. Look at after summing the two matrices -> different
    % 2. Can I use DISTINCT to sparsely solve the T2s?
        % 3a. If the answer is no, then can we use kronecker product on the server?
        % 3b. If the answer is no, how do we incorporate DISTINCT to solve for LED side?
    % 4. If the answer is yes, how do we do LEDs in the 3rd dimension?
        %can we apply DISTINCT for 3d?

%% STEP 0: paths
Folder = fullfile(cd, '..');
addpath(genpath(fullfile(Folder)));

%add the whole folder on MAC here:
addpath(genpath('D:\Dropbox\96 Paper DISTINCT\2 algorithm verification'))
addpath(genpath('/proj2/adriaantaal/96 Paper DISTINCT'))

%% STEP1: Load calibration of matrices
calibDir = 'Device/Calibration_2020_10_13_ECI_WB/calibMat/';
calibFile = 'calibMat_noGap_Cal_ECI_Filter_200um';
calibFile = [calibDir, filesep, calibFile, '.mat'];
calibData = load(calibFile);

[P1, Q1] = deal(calibData.P1, calibData.Q1);
[P2, Q2] = deal(calibData.P2, calibData.Q2);
M_net = calibData.M_net;
gapFlg = calibData.gapFlg;

% sensorSize = calibData.sensorSize;
% sceneSize = calibData.sceneSize;

% %crop the matrix for useful computation. square please for simplicity
% sensorSize = calibData.sensorSize/2; M = prod(sensorSize);
% sceneSize = ceil(calibData.sceneSize/2); 
% sceneCropOffset = (ceil(calibData.sceneSize/2) - sceneSize)/2;
% matSelScene = (1+sceneCropOffset):sceneSize(1)+sceneCropOffset;
% matSelSensor = 1:sensorSize(1);

%or: decimate the matrix of the whole imager
% dec = 2;
% sensorSize = calibData.sensorSize/2; M = prod(sensorSize);
% sceneSize = ceil(calibData.sceneSize/2); N = prod(sceneSize);
% sceneCropOffset = (ceil(calibData.sceneSize/2) - sceneSize)/2;
% matSelScene = (1+sceneCropOffset):dec:sceneSize(1)*dec+sceneCropOffset;
% matSelSensor = 1:dec:sensorSize(1)*dec;

dec = 1;
sensorSize = [64 64]; M = prod(sensorSize);
sceneSize = [128 128]; N = prod(sceneSize);
sceneCropOffset = floor((ceil(calibData.sceneSize) - sceneSize)/2);
matSelScene = (1+sceneCropOffset):dec:sceneSize(1)*dec+sceneCropOffset;
sensorCropOffset = floor((ceil(calibData.sensorSize) - sensorSize)/2);
matSelSensor = (1+sensorCropOffset):dec:sensorSize(1)*dec+sensorCropOffset;


%physical locations of the scene voxels in um
Dpitch = 30;
% X = ((0:sceneSize-1)-(sceneSize(1)-sensorSize(1))/2)*Dpitch;
% Y = ((0:sceneSize-1)-(sceneSize(2)-sensorSize(2))/2)*Dpitch;
X = matSelScene*Dpitch;
Y = matSelScene*Dpitch;
Z = 400;

P1 = P1(matSelSensor,matSelScene);
Q1 = Q1(matSelSensor,matSelScene);
P2 = P2(matSelSensor,matSelScene);
Q2 = Q2(matSelSensor,matSelScene);
M_net = M_net(matSelSensor,matSelSensor);

PTP = P2'*P2; QTQ = Q2'*Q2;

A = Phi_IBIS(P1,Q1,P2,Q2,M_net,'dir');
At = Phi_IBIS(P1,Q1,P2,Q2,M_net,'adj');

% Generate an E-matrix from the LED locations
ledmaskDir = 'Device/Calibration_2019_10_04/LED_mask/';
ledmaskFile = 'LED_mask.png';
M_led = imread([ledmaskDir, ledmaskFile]);
% figure; imshow(M_net);

%lay down the LED locations
Elocs = []; Epitch = 32; Ewidth = 8; 
% Eoff = 0; 
Eoff = 9; 
kk = 0;
for jj = 0:9
    for ee = 0:4
        kk = kk+1;
        Elocs(kk,:) = [Eoff+Ewidth+Epitch*(0.5*mod(jj,2)+ee) Eoff+Ewidth+Epitch*0.5*jj];
    end
end
Elocsxy = Elocs*Dpitch;
B = size(Elocs,1);

figure; scatter3(Elocsxy(:,1),Elocsxy(:,2),Elocsxy(:,2)*0); axis([X(1) X(end) Y(1) Y(end)])

%emitter light field parameters, read them from a far field simulation plot
sigma_main = ones(B,1); mmain = zeros(B,1); mside = zeros(B,1); elevation = ones(B,1); %radians
sigma_elevation = ones(B,1); azimuth = ones(B,1); sigma_azimuth = ones(B,1);

% sigma_main = ones(B,1)*1.3; 
% sigma_azimuth = ones(B,1)*1.55; sigma_elevation = ones(B,1)*1.5; %radians
% mmain = ones(B,1)*100; mside = ones(B,1)*2000; 
% elevation = ones(B,1)*1; %radians
% azimuth = linspace(0,pi/2,B)'; %azimuth = azimuth([1 9 5 13 2 10 6 14 3 11 7 15 4 12 8 16]); %scramble them for maximum effect

deadangle = ones(B,1)*90/180*pi;
Eparams = [sigma_main mmain elevation sigma_elevation azimuth sigma_azimuth mside deadangle];

%emitter matrix
E = Egen(X,Y,Z,Eparams,B,Elocsxy);
E = E./max(max(E)); 
% figure; imagesc(X,Y,reshape(E(1,:),sceneSize))
% figure; imagesc(X,Y,reshape(E(B,:),sceneSize))


%% STEP1b: generate the derivated matrices
D = kron(Q1,P1)+kron(Q2,P2);
DTD = D'*D;
[M,N] = size(D);

% Dpinv = pinv(D);
% projD = Dpinv*D;        %detector projection matrix


%generate Ek from combinations
% comb_haar = comb_haar_6x6();
% comb_haar = eye(B);
% Ek = comb_haar*E;

%generate E from true Haar basis, select K before to guarantee orthogonality
% K = 1024;
% [Ek,x] = wpfun('dct',K,ceil(log2(N)));
% Ek(1,:) = []; Ek(:,[1 2]) = []; 
% Ek = Ek - min(min(Ek)); Ek = Ek./max(max(Ek));

%or make a DCT basis
K = 128;
Ek = dctmtxnonsquare(N,K+1); Ek(1,:) = [];
Ek = Ek - min(min(Ek)); Ek = Ek./max(max(Ek));

ETE = Ek'*Ek;
% Ekpinv = pinv(Ek);
% projEk = Ek'*Ekpinv';   %emitter projection matrix

%uniform illumination
% dotEk = Ek(1,:)'; 
dotEk = ones(N,1); 
DdotEk = D.*dotEk';
% EdotTE = dotEk*dotEk';

%inspect the light fields
% figure; imagesc(X,Y,reshape(dotEk,sceneSize))

%inspect the Hessian qualities
H = DTD.*ETE;
% HdotE = DTD.*EdotTE;

% figure; set(gcf,'Color',[1 1 1]); imagesc(H); title('IBIS detector*emitter Hessian'); 
% figure; imagesc(DTD); title('Detector matrix inner product')
% figure; imagesc(ETE); title('Emitter matrix inner product')

% figure; imagesc(EdotTE); title('Blanket Emitter matrix inner product')
% figure; set(gcf,'Color',[1 1 1]); imagesc(HdotE); title('IBIS blanket illumination Hessian'); 


% %% STEP 1c: analyze all matrices 
% %first derivative: Jacobian
% 
% %second derivative: Hessian
% %For DISTINCT formulation, X is diagonal

%look at eigenvalues
% eigenval(:,1) = eig(DTD); eigenval(:,2) = eig(ETE); eigenval(:,3) = eig(H); 
% eigenval = flipud(eigenval); 
% % eigenval(:,1) = flipud(eigenval(:,1));
% eigenval = eigenval./eigenval(1,:); %normalize
% negeigenval = real(eigenval); negeigenval(negeigenval>0) = 0; 
% poseigenval = real(eigenval); poseigenval(poseigenval<0) = 0; 
% lpos = sum(poseigenval>0)
% lneg = sum(negeigenval<0)
% lpos - lneg
% sum(poseigenval)./sum(-negeigenval)
% sum(abs(real(eigenval)))./sum(abs(imag(eigenval)))
% sumbalance = sum(poseigenval(M:end))./sum(-negeigenval)
% 
% f1 = figure; set(gcf,'Color',[1 1 1]); set(gcf,'Units','Inches'); set(gcf,'Position',[1 1 3.5 3]);
% loglog(abs(eigenval)); 
% legend({['DTD, M = ' num2str(M)],['ETE, K = ' num2str(K)],['H, N = ' num2str(N)],['symmH, N = ' num2str(N)]},'Location','Southwest');
% ylabel('Eigenvalue magnitude'); xlabel('Eigenvalue number'); 
% title(['Kronecker D and ' num2str(K) ' Haar wavelets']);
% ylim([1e-30 1]);
% savefig(f1,['eigenvalues_IBIS_DISTINCT_Haar' num2str(K) '_norm']);
% print(f1,['eigenvalues_IBIS_DISTINCT_Haar' num2str(K) '_norm'],'-dpng','-r600');


% %eigenvalues for the separable imaging
% %for the IBIS formulation, X is dense
% Hs = (P1'*P1)*(Q1'*Q1) + (P1'*P2)*(Q2'*Q1) + ...
%      (P2'*P1)*(Q1'*Q2) + (P2'*P2)*(Q2'*Q2); 
% ss = symmetricity(Hs);
% % eigenvalsep(:,1) = eig(PTP); eigenvalsep(:,2) = eig(QTQ); eigenvalsep(:,3) = eig(Hs);
% % eigenvalsep(:,[1:2]) = flipud(eigenvalsep(:,[1:2]));
% negeigenvalsep = real(eigenvalsep); negeigenvalsep(negeigenvalsep>0) = 0; 
% poseigenvalsep = real(eigenvalsep); poseigenvalsep(poseigenvalsep<0) = 0; 
% lpos = sum(poseigenvalsep>0)
% lneg = sum(negeigenvalsep<0)
% lpos - lneg
% sum(poseigenvalsep)./sum(-negeigenvalsep)
% sum(abs(real(eigenvalsep)))./sum(abs(imag(eigenvalsep)))
% sumbalance = sum(poseigenvalsep(1:end))./sum(-negeigenvalsep)
% f2 = figure; set(gcf,'Color',[1 1 1]); set(gcf,'Units','Inches'); set(gcf,'Position',[1 1 3.5 3]);
% loglog(real(eigenvalsep)); 
% legend({['PTP, M = ' num2str(sqrt(M))],['QTQ, K = ' num2str(sqrt(M))],['Hs, N = ' num2str(sqrt(N))],['symmHs, N = ' num2str(sqrt(N))]},'Location','Southwest');
% ylabel('Eigenvalue magnitude'); xlabel('Eigenvalue number'); 
% title(['Separable imager, no LEDs']);
% savefig(f2,'eigenvalues_IBIS');
% print(f2,'eigenvalues_IBIS','-dpng','-r600');
% 
% 
% %% STEP 1d: PSF analysis
% 
% %location of virtual source
% XYZpsf = [(X(end)-X(1))/2+X(1) (Y(end)-Y(1))/2+Y(1) (Z(end)-Z(1))/2+Z(1)];
% [xpsf, ~, xind] = convert_fluoloc(X,Y,Z,XYZpsf);
% 
% %IBIS+LED
% % psfg = projD(:,xind).*projEk(xind,:)';
% psfg = projH(:,xind);
% PSFg = reshape(psfg,sceneSize);
% xyPSFg = PSFg(:,:,1)';
% xPSFg = PSFg(:,ceil(length(Y)/2),1)';
% fwhm(X,xPSFg)
% 
% %IBIS only
% PSFa = reshape(projD(:,xind),sceneSize);
% xyPSFa = PSFa(:,:,1)';
% xPSFa = PSFa(:,ceil(length(Y)/2),1)';
% fwhm(X,xPSFa)
% 
% figure
% imagesc(X,Y,xyPSFg)
% set(gca,'YDir','normal')
% % set(gca,'colormap',greencolormapHC)
% set(gca,'LineWidth',1.5)
% title('IBIS + LED')
% 
% figure
% imagesc(X,Y,xyPSFa)
% set(gca,'YDir','normal')
% % set(gca,'colormap',greencolormapHC)
% set(gca,'LineWidth',1.5)
% title('IBIS')
% 
% figure; hold on;
% plot(X,xPSFg,'g'); plot(X,xPSFa,'k');
% 
% %% STEP2a: load a simulated measurement image and calculate ideal data
% % res target per px matches sensor px pitch
% img = imread('test/resBars_6px.png'); 
% img = im2double(img(:,:,1));
% img = padarray(img, [20,20]);
% img = imresize(img,sceneSize,'nearest');
% L = numel(find(img));
% 
% %original direct matrix multiplication
% yorig = P1*img*Q1' + P2*img*Q2';
% 
% %kronecker product
% xtrue = tovec(img);
% Xtrue = diag(xtrue);
% y = D*diag(xtrue)*Ek';
% gradconstdiag = diag(D'*y*Ek);

%% STEP2b: load a simulated sparse image

%% STEP2c: load a true measurement

%% STEP3: Reconstruct, compare nesterovGradDesc with DISTINCT
rng(numBatchFile);
%monte carlo params
L=[2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32 64 128 256 512 1024];          %number of sources    
SNR = 100;      %for implementation of noise
nIter = 20;   %number of repetitions in the monte carlo

%nestor params
lmbd2 = 1e-3;
tol = 1e-8;

%DISTINCT params
maxIter = 1e3;
maxIterInt = 1e4;
lambdadiv = 10;

%parameters for SL0
mu_0 = 2;
q = 3;
sigma_min = 0.02; %As a rule of tumb, for the noisy case, sigma_min should be about 2 to 4 times of the standard deviation of this noise
sigma_decrease_factor = 0.5;
max_contrast_ratio = 1e3;

for ll = 1:length(L) 
    for ii = 1:nIter
        [~,~,XYZtrue,~,~,xtrue] = create_randomRGlocs2D(X,Y,Z,L(ll),pi/2,X(1),Y(1));
        %2D detector-emitter formulation
        y = D*diag(xtrue)*Ek';
        noise = 0; y = y+noise;
        gradconstdiag = diag(D'*y*Ek);
%         sl0gradconstdiag = diag(Dpinv*y*Ekpinv');
        
        %1D formulation
%         ydotEk = DdotEk*xtrue+noise; %equivalent

        %original formulation
%         img = reshape(xtrue.*dotEk,sceneSize);
%         yorig = (P1*img*Q1' + P2*img*Q2').*M_net;
        
        %derived parameters
        addSize = max(floor(L(ll)/2),1);
        
        %In the noiseless case, we can use SVD to guess number of sources!
        [U,S,V] = svd(y,'econ'); U = U(:,1:L(ll)); V = V(:,1:L(ll)); S = S(1:L(ll),1:L(ll));
        s = diag(S);
        
    %In the noiseless case, we can use SVD to guess number of sources!
%         s = svd(y,'econ'); sSNR = wthresh(s,'s',s(1)/SNR); Lhat = sum(sSNR>0); 
%         gcdsort = sort(gradconstdiag,'descend'); lambda = gcdsort(Lhat)/lambdadiv;
%         gcdotEksort = sort(ydotEk'*DdotEk,'descend'); lambdadotEk = gcdotEksort(Lhat)/lambdadiv;
        
    lambda = max(gradconstdiag)/lambdadiv;
%     lambdadotEk = max(ydotEk'*DdotEk)/lambdadiv;

    %a. DISTINCT: 
        tic; [xhatDS,~,~,~] = DISTINCTv4(y,D,Ek,lambda,H,gradconstdiag,U,s,V,0,maxIter,maxIterInt,addSize);
        tcompute(1,ll,ii) = toc;
        disp(['Found ' num2str(numel(find(xhatDS))) ' locations'])
        %choose L largest entries -> for exact comparison wth ground truth
        [~,DSind] = sort(xhatDS,'descend'); xhatDS(DSind(L(ll)+1:end)) = 0;
        
%         dxe = calc_tri_opposite(D,diag(xhatDS),Ek',find(xhatDS));
%         J = norm(y-dxe,'fro')^2+lambda*norm(xhatDS,1)
%         Jtr = norm(y,'fro')^2 + xhatDS'*H*xhatDS - 2*gradconstdiag'*xhatDS + lambda*norm(xhatDS,1)
        

%     %b. ICopt, no extra emitter rows, faster but uses less information
%         %so its expected to be far less accurate
%         tic; [xhatIC, ~] = ICOpt(DdotEk,ydotEk,lambdadotEk,0,Inf,0,1,addSize);
%         tcompute(2,ll,ii) = toc;    
%         disp(['Found ' num2str(numel(find(xhatIC))) ' locations'])
%         %choose L largest entries -> for exact comparison wth ground truth
%         [~,ICind] = sort(xhatIC,'descend'); xhatIC(ICind(L(ll)+1:end)) = 0;


%     %c. 2D SL0, adapted for 2D. Already improved speedup ~2-3x 
%         %over the original authors up but still slow as balls
%         %the second speedup round -> speeds up by another 20x
%         tic; xhatSL0=SL0_2D_ATv3(y, sigma_min, sigma_decrease_factor, mu_0, q, Dpinv, Ekpinv, projD, projEk, sl0gradconstdiag, max_contrast_ratio);
%         tcompute(3,ll,ii) = toc;    
%         %choose L largest entries -> for exact comparison wth ground truth
%         [~,SL0ind] = sort(xhatSL0,'descend'); xhatSL0(SL0ind(L(ll)+1:end)) = 0;   


%     %d. competitor #3: 2D NIHT algorithm
%     %Slow because it relies on calculation of many matrix norms
%         tic; [xhatNIHT] = NIHT2D_ATv4(y,D,Ek,L(ll));
%         tcompute(4,ll,ii) = toc;    
%         %choose L largest entries -> for exact comparison wth ground truth
%         [~,NIHTind] = sort(xhatNIHT,'descend'); xhatNIHT(NIHTind(L(ll)+1:end)) = 0;   


    %e. competitor #4 -> our own algorithm using the 
    %gradient descent internal solver, for comparison
        tic; [xhatDSGD,f1,f2,f1f2] = DISTINCTv4(y,D,Ek,lambda,H,gradconstdiag,U,s,V,1,maxIter,maxIterInt,addSize);
        tcompute(2,ll,ii) = toc;
        disp(['Found ' num2str(numel(find(xhatDSGD))) ' locations'])
        %choose L maximum -> for exact comparison wth ground truth
        [~,DSGDind] = sort(xhatDSGD,'descend'); xhatDSGD(DSGDind(L(ll)+1:end)) = 0;

        
%     %f. competitor #5 -> original Nesterov gradient descend method
%         tic; xhatNEST = nesterovGradDesc_AT_v2(yorig,A,At,lmbd2,maxIter,tol);
%         tcompute(6,ll,ii) = toc;
%         disp(['Found ' num2str(numel(find(xhatNEST))) ' locations'])
%         %choose L largest entries -> for exact comparison wth ground truth
%         [~,NESTind] = sort(xhatNEST,'descend'); xhatNEST(NESTind(L(ll)+1:end)) = 0;
        
        
%     %g. DISTINCT using direct inversion: 
%         tic; [xhatDINV,f1,f2,f1f2] = DISTINCTv4(y,D,Ek,lambda,H,gradconstdiag,U,s,V,2,maxIter,maxIterInt,addSize);
%         tcompute(7,ll,ii) = toc;
%         disp(['Found ' num2str(numel(find(xhatDINV))) ' locations'])
%         %choose L largest entries -> for exact comparison wth ground truth
%         [~,DINVind] = sort(xhatDINV,'descend'); xhatDINV(DINVind(L(ll)+1:end)) = 0;
        
%       %h. DISTINCT using linsolve: 
%         tic; [xhatDSQP,f1,f2,f1f2] = DISTINCTv4(y,D,Ek,lambda,H,gradconstdiag,U,s,V,3,maxIter,maxIterInt,addSize);
%         tcompute(8,ll,ii) = toc;
%         disp(['Found ' num2str(numel(find(xhatDSQP))) ' locations'])
%         %choose L largest entries -> for exact comparison wth ground truth
%         [~,DSQPind] = sort(xhatDSQP,'descend'); xhatDSQP(DSQPind(L(ll)+1:end)) = 0;
        
        
    %calculate the accuracy of each solution compared to ground truth.
    %for IBIS, the answer can be negative
        accuracy(1,ll,ii) = some_accuracy_function(abs(xhatDS),xtrue);
%         accuracy(2,ll,ii) = some_accuracy_function(abs(xhatIC),xtrue);
%         accuracy(3,ll,ii) = some_accuracy_function(abs(xhatSL0),xtrue);
%         accuracy(4,ll,ii) = some_accuracy_function(abs(xhatNIHT),xtrue);
        accuracy(2,ll,ii) = some_accuracy_function(abs(xhatDSGD),xtrue);
%         accuracy(6,ll,ii) = some_accuracy_function(abs(xhatNEST(:)),xtrue);
%         accuracy(7,ll,ii) = some_accuracy_function(abs(xhatDINV),xtrue);   
%         accuracy(8,ll,ii) = some_accuracy_function(abs(xhatDSQP),xtrue);   
    end
end


accmean = nanmean(accuracy,3);
accstd = std(accuracy,0,3,'omitnan');

tcomputemean = nanmean(tcompute,3);
tcomputestd = std(tcompute,0,3,'omitnan');


%create directory
t = datetime('now');
mmdd = datestr(t);
mkdir(['myoutput_' mmdd]);

save(['myoutput_' mmdd '/DISTINCTonly_IBIS_DCT128.mat'],'accuracy','tcompute','L','SNR','nIter','addSize','lambdadiv','q','sigma_min','sigma_decrease_factor','max_contrast_ratio','X','Y','Z','numBatchFile','K');


%create figure accuracy
f1 = figure;
plot(L, accmean)
xlim([0 1024])
ylim([0 1])
xlabel("Number of Sources (L)")
ylabel("Accuracy")
legend('DS','ICOpt','SL0','NIHT','DSGD','nesterov','DINV','DSQP')
savefig(f1,['myoutput_' mmdd '/AccuracyGraph']) ;


f2 = figure;
errorbar(L,accmean(1,:),accstd(1,:))
hold on;
errorbar(L,accmean(2,:),accstd(2,:))
errorbar(L,accmean(3,:),accstd(3,:))
errorbar(L,accmean(4,:),accstd(4,:))
errorbar(L,accmean(5,:),accstd(5,:))
errorbar(L,accmean(6,:),accstd(6,:))
errorbar(L,accmean(7,:),accstd(7,:))
hold off;
xlim([0 1024]); set(gca,'xscale','log');
xlabel("Number of Sources (L)"); ylabel("Accuracy");
legend('DS','ICOpt','SL0','NIHT','DSGD','nesterov','DINV','DSQP')
savefig(f2,['myoutput_' mmdd '/AccuracyGraphWithError']) ;


%create figure computing time
f4 = figure;
plot(L, mean(tcompute,3))
xlim([0 1024]); set(gca,'xscale','log');  set(gca,'yscale','log');
xlabel("Number of Sources (L)"); ylabel("Computing Time(s)");
legend('DS','ICOpt','SL0','NIHT','DSGD','nesterov','DINV','DSQP');
savefig(f4,['myoutput_' mmdd '/ComputingTime']) ;



% %% STEP4: show results
% figure; set(gcf,'Color',[1 1 1]);
% subplot(2,3,1); imagesc(X,Y,img); title('Original');
% subplot(2,3,2); imagesc(X,Y,reshape(xhatDS,sceneSize)); title('DISTINCT')
% % subplot(2,3,2); imagesc(X,Y,reshape(xhatDINV,sceneSize)); title('DINV')
% subplot(2,3,3); imagesc(X,Y,reshape(xhatIC,sceneSize)); title('ICOpt')
% % subplot(2,3,4); imagesc(X,Y,reshape(xhatSL0,sceneSize)); title('SL0')
% subplot(2,3,5); imagesc(X,Y,reshape(xhatDSGD,sceneSize)); title('DSGD')
% subplot(2,3,6); imagesc(X,Y,xhatNEST); title('Nesterov');




