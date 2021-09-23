%% Welcome to the code that simulates DISTINCT 
%First create a model for Detector and Emitter array
%Then apply the equations to verify reconstruction

close all; clear all; clc;

%% STEP 0: import all functions and folders
Folder = fullfile(cd, '..');
addpath(genpath(fullfile(Folder)));

%add the whole folder on MAC here:
addpath(genpath('D:\Dropbox\96 Paper DISTINCT\2 algorithm verification'))

addpath(genpath('/proj2/adriaantaal/96 Paper DISTINCT'))

%% STEP 1: make a voxel grid and probe layout


%% STEP 2: do the image compression


%% option 1: blind decomposition of image Y into 3 matrices (DXE)
    %this only works for separate images, as DXE will be different for each image
    %compression is determined by numel(D,x,E)
        %a sparser x gives higher compression because it needs less rows and columns in D & E
    %quality is determined by inspection and by PSNR
    
    %comparison is made to the competitors:
        % SVD
        % 2-matrix NNMF
        % Dense NNMTF (mult + als)

        
%step 1. load the image
% y = double(imread('corn.tif',3));
y = double(imread('rice.png'));
% y = y(1:16,1:16);
[M,K] = size(y);

%step 2. confirm with DCT
Dleft = dctmtx(M);
Dright = dctmtx(K);

yDCT = (Dleft*y*Dright');

figure; imshow(uint8(idct2(yDCT)));

%step 3. confirm with Haar
Eleft=haarmtx(M);
Eright=haarmtx(K);

yHaar = (Eleft*y*Eright');

figure; imshow(uint8(Eleft'*yHaar*Eright));


%% now we want to compress the image using a left DCT and right DCT matrix
% D0 = Dleft;
% E0 = Dright;
[D0,~,E0] = svd(y);
intSolver = 0; % intsolver 0 is reliable
lambda = max(diag(D0'*y*E0))/5000;
maxIter = 50;
maxIterInt = 1e4;
addSize = 10;

%method 1: sparse decomposition of Y with pre known D and E
tic; [xhatDSTF,DF,EF,f1,f2,f1f2] = DISTINCT_NNMTFv2(y,lambda,D0,E0,intSolver,maxIter,maxIterInt,addSize);
tcompute(1) = toc; L = sum(xhatDSTF>0); disp(['Cardinality of x = ' num2str(L)])
storage(1) = numel(DF)+numel(EF)+numel(xhatDSTF);

%method 2: O-NNMF, same rank as L
in_options = [];
in_options.x_init.H = E0(1:L,:); 
% in_options.x_init.W = D0(:,1:L);
tic; [orthnmf, infos] = nmf_orth_mu(y, L, in_options); Ahat = orthnmf.W; Shat = orthnmf.H;
tcompute(2) = toc;
storage(2) = numel(Ahat)+numel(Shat)

%method 3: jpeg
quality = 25;
tic; [ycompJPEG, dct_quantized, dct_qsparse, DCTcols] = jpeg(y,quality); 
    [Dout,~,Eout] = svd(DCTcols); lambda = max(diag(Dout'*DCTcols*Eout))/500;
tcompute(3) = toc;
storage(3) = nnz(dct_qsparse);

%method 4: compressed JPEG. Do a low rank approximation to the DCT components
quality = 50;
tic; [ycompJPEG, dct_quantized, dct_qsparse, DCTcols] = jpeg(y,quality); 
    [Dout,~,Eout] = svd(DCTcols); 
%     Dout = tovec(dct(eye(8)));
%     Eout = Dout'*DCTcols;
    lambda = max(diag(Dout'*DCTcols*Eout))/1e6;
    [xhatDCT,Ddct,Edct,f1,f2,f1f2] = DISTINCT_NNMTFv2(DCTcols,lambda,Dout,Eout,intSolver,maxIter,maxIterInt,addSize);
    L = sum(xhatDCT>0); disp(['Cardinality of x = ' num2str(L)])
    [QX] = jpeg_quality(quality);
    dctcols_DISTINCT = Ddct*diag(xhatDCT)*Edct';
    [ycompJPDS] = jpeg_decodecols(dctcols_DISTINCT,QX);
tcompute(4) = toc;
storage(4) = nnz(Ddct)+nnz(Edct)+nnz(xhatDCT);



%retrieve the images
ycompDSTF = DF*diag(xhatDSTF)*EF';  mse(1) = norm(y-ycompDSTF);
ycompONMF = Ahat*Shat;              mse(2) = norm(y-ycompONMF);
        

%show the images
figure; set(gcf,'Color',[1 1 1]);
subplot(2,4,1); imshow(uint8(ycompDSTF)); subplot(2,4,5); title('DSTF'); imagesc(real(log10(dct2(ycompDSTF))));
subplot(2,4,2); imshow(uint8(ycompONMF)); subplot(2,4,6); title('ONMF'); imagesc(real(log10(dct2(ycompONMF))));
subplot(2,4,3); imshow(uint8(ycompJPEG)); subplot(2,4,7); title('JPEG'); imagesc(real(log10(dct2(ycompJPEG))));
subplot(2,4,4); imshow(uint8(ycompJPDS)); subplot(2,4,8); title('JPEG_DISTINCT'); imagesc(real(log10(dct2(ycompJPDS))));

%calculate PSNR
peaksnr(1) = psnr(ycompDSTF,y);
peaksnr(2) = psnr(ycompONMF,y);
peaksnr(3) = psnr(ycompJPEG,y);
peaksnr(4) = psnr(ycompJPDS,y)

compression = numel(y)./storage



%% option 2: sparse decomposition of Y with pre known D and E
    %this works for a bunch of images, as D&E are constant
    %only X will be different for images, and is no longer diagonal!
    %compression is determined by numel(X)
    %quality is determined by inspection and by PSNR
    
    %the DCT and Haar matrices are adapted to MxN and NxK respectively
    
    %comparison is made to the competitors:
        % SVD
        % 2-matrix Haar compression
        % 2-matrix DCT compression
        % Dense NNMTF (mult + als)        
        



