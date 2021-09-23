%% Welcome to the code that simulates DISTINCT 
%First create a model for Detector and Emitter array
%Then apply the equations to verify reconstruction

close all; clear all; clc;
ncolor = 50;
greencolormapHC = ([54*linspace(0,1,ncolor)' 255*linspace(0,1,ncolor)' zeros(ncolor,1)]/255).^1; 
redcolormapHC = ([255*linspace(0,1,ncolor)' 0*linspace(0,1,ncolor)' zeros(ncolor,1)]/255).^1;

%% STEP 0: import all functions and folders
Folder = fullfile(cd, '..');
addpath(genpath(fullfile(Folder)));

%add the whole folder on MAC here:
addpath(genpath('D:\Dropbox\96 Paper DISTINCT\2 algorithm verification'))

addpath(genpath('/proj2/adriaantaal/96 Paper DISTINCT'))

%% STEP 1: make a voxel grid and probe layout

%create probe layout
NPixel = 1024;
xpitch = 24.6;
ypitchD = 95.96;
ypitchE = 31;
DProbe = create_Probe(NPixel,xpitch,ypitchD,259.4,1);
EProbe = create_Probe(NPixel,xpitch,ypitchE,259.4,1);
EProbe(:,2) = EProbe(:,2) + (ypitchD-ypitchE)/2; %shift emitters to center of shank

xres = 50;
yres = 50;
zres = 50;
X = (-250+min(DProbe(:,1))):xres:(max(DProbe(:,1))+250);
Y = (-250+min(DProbe(:,2))):yres:(max(DProbe(:,2))+250);
Z = 200:zres:300;
Z = 400;
sizes = [length(X) length(Y) length(Z)]; prod(sizes)


%% STEP 2 make a spatial model for the emission and detection profile

%MIMAS 'emitter' matrix
addpath(genpath('D:\Dropbox\Graduation Research\23 MIM measurements results\7. Pinhole'))
DCR = load('MIMASs5_greenagarpinhole0_DCR.mat'); DCR = DCR.DCR;
DATA = load('MIMASs5_greenagarpinhole_center4.mat'); DATA = DATA.DATA;
DATAr = load('MIMASs5_redagarpinhole_center5.mat'); DATAr = DATAr.DATA;
[hotdeadpix] = Probe_HP_DP(DCR.SPAD_COUNTS, 5, 0); probenhp = 1:1024;
probenhp(unique([hotdeadpix])) = [];
load('D:\Dropbox\Graduation Research\23 MIM measurements results\7. Pinhole\20201213_MIMAS_pinholesagar\favourites\best\MIMASP_h0pad_data_mu1lambdadiv35_L16.mat')
timestart = 8;
YTG = squeeze(DATA.SPAD_COUNTS(probenhp,1,timestart,DATA.mm)).*greentime(timestart:end)/greentime(timestart);
YTGr = squeeze(DATAr.SPAD_COUNTS(probenhp,1,timestart,DATAr.mm)).*redtime(timestart:end)/redtime(timestart);
y = YTG+YTGr;

%MIMAS 'emitter' matrix
Ek = [repmat(mean(YTG)./max(mean(YTG)),prod(sizes),1);repmat(mean(YTGr)./max(mean(YTGr)),prod(sizes),1)]';

%detector light field parameters
clustersize = 16;
probeparams = load('QP2020MIMProbeParams.mat');
MIMparams = probeparams.MIMparams;
alphag = MIMparams(:,1);
betag = MIMparams(:,2);
alphar = MIMparams(:,3);
betar = MIMparams(:,4);
modindex = MIMparams(:,5);
chi = MIMparams(:,6)/180*pi;     %convert to radians
mperp =     zeros(clustersize,1); %dont use
gammaperp = zeros(clustersize,1); %dont use
mpara =     zeros(clustersize,1); %dont use
gammapara = zeros(clustersize,1); %dont use
deadangle = ones(clustersize,1)*pi/2;
Dgparams = [alphag betag modindex chi mperp gammaperp mpara gammapara deadangle];
Drparams = [alphar betar modindex chi mperp gammaperp mpara gammapara deadangle];

%detector matrix
Dg = Agenv7(X,Y,Z,repmat(Dgparams,NPixel/size(Dgparams,1),1),NPixel,DProbe);
Dg = Dg(probenhp,:);
Dg = Dg./max(max(Dg)); 

Dr = Agenv7(X,Y,Z,repmat(Drparams,NPixel/size(Drparams,1),1),NPixel,DProbe);
Dr = Dr(probenhp,:);
Dr = Dr./max(max(Dr));

D =[Dg Dr];


%% generate the derived matrices
Dpinv = pinv(D);
Ekpinv = pinv(Ek);
projD = Dpinv*D;        %detector projection matrix
projDg = pinv(Dg)*Dg; projDr = pinv(Dr)*Dr;
projEk = Ek'*Ekpinv';   %emitter projection matrix
ETE = Ek'*Ek;
DTD = D'*D;
figure; imagesc(ETE); title('')
figure; imagesc(DTD); title('Detector matrix inner product')
figure; imagesc(D*Ek')

%% STEP 3 reconstruct using each method
%data parameters
L=[100];          %number of sources    

%parameters for the solvers
maxIter = 99;
maxIterInt = 1e4;
addSize = 5;
intSolver = 1;
lambdadiv = 1e2;

% figure; semilogy(svd(y),'x');

%a. DISTINCT: 
gradconstdiag = diag(D'*y*Ek);
lambda = max(gradconstdiag)/lambdadiv;
[xhatDS,f1,f2,f1f2] = DISTINCT(y,D,Ek,lambda,ETE,DTD,gradconstdiag,intSolver,maxIter,maxIterInt,addSize);
disp(['Found ' num2str(numel(find(xhatDS))) ' locations'])
%choose L largest entries -> for exact comparison wth ground truth
% [~,DSind] = sort(xhatDS,'descend'); xhatDS(DSind(L+1:end)) = 0;

% xhatDS = diag(Dpinv*y*Ekpinv');

%find the color indices
xhatg = xhatDS(1:prod(sizes));
xhatr = xhatDS(prod(sizes)+1:end);

xhatgCube = reshape(xhatg,sizes);
xhatrCube = reshape(xhatr,sizes);

greenIVI = reshape(projDg*xhatg,sizes);
redIVI = reshape(projDr*xhatr,sizes);

figure; set(gcf,'Color',[1 1 1]);
subplot(2,2,1); 
imagesc(X,Y,real(squeeze(greenIVI(:,:,1))'))
set(gca,'colormap',greencolormapHC); set(gca,'YDir','normal');
title('Green, smoothed')
subplot(2,2,2)
imagesc(X,Y,squeeze(sum(xhatgCube,3))')
set(gca,'colormap',greencolormapHC); set(gca,'YDir','normal');
title('Green, localized')
subplot(2,2,3)
imagesc(X,Y,real(squeeze(redIVI(:,:,1))'))
set(gca,'colormap',redcolormapHC); set(gca,'YDir','normal');
title('Red, smoothed')
subplot(2,2,4)
imagesc(X,Y,squeeze(sum(xhatrCube,3))')
set(gca,'colormap',redcolormapHC); set(gca,'YDir','normal');
title('Red, localized')
set(gcf,'PaperPositionMode','auto')
figname = ['MIMASP_DISTINCT_IVI_lambdadiv' strrep(num2str(lambdadiv),'.','') '_L' num2str(L)];



