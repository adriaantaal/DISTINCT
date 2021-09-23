%% Welcome to the code that simulates DISTINCT 
%First create a model for Detector and Emitter array
%Then apply the equations to verify reconstruction

close all; clear all; clc;

%% STEP 0: import all functions and folders
Folder = fullfile(cd, '..');
addpath(genpath(fullfile(Folder)));

%add the whole folder here:
addpath(genpath('D:\Dropbox\96 Paper DISTINCT\2 algorithm verification'))

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
sizes = [length(X) length(Y) length(Z)]; prod(sizes)


%% STEP 2 make a spatial model for the emission and detection profile
clustersize = 16; %16 different types, all 

%emitter light field parameters, read them from a far field simulation plot
sigma_main = ones(clustersize,1)*0.5; %radians
mmain = ones(clustersize,1)*100; 
mside = ones(clustersize,1)*1000;
elevation = ones(clustersize,1)*1; %radians
sigma_elevation = ones(clustersize,1)*0.5; %radians
azimuth = linspace(0,pi/2,clustersize)';
azimuth = azimuth([1 9 5 13 2 10 6 14 3 11 7 15 4 12 8 16]); %scramble them for maximum effect
sigma_azimuth = ones(clustersize,1)*1.5;
deadangle = ones(clustersize,1)*60/180*pi;
Eparams = [sigma_main mmain elevation sigma_elevation azimuth sigma_azimuth mside deadangle];

%emitter matrix
E = Egen(X,Y,Z,repmat(Eparams,NPixel/size(Eparams,1),1),NPixel,EProbe);
E = E./max(max(E)); 

%detector light field parameters
probeparams = load('QP2020MIMProbeParams.mat');
MIMparams = probeparams.MIMparams;
alphag = MIMparams(:,1);
betag = MIMparams(:,2);
modindex = MIMparams(:,5);
chi = MIMparams(:,6)/180*pi;     %convert to radians
mperp =     zeros(clustersize,1); %dont use
gammaperp = zeros(clustersize,1); %dont use
mpara =     zeros(clustersize,1); %dont use
gammapara = zeros(clustersize,1); %dont use
deadangle = ones(clustersize,1)*pi/2;
Dgparams = [alphag betag modindex chi mperp gammaperp mpara gammapara deadangle];

%detector matrix
D = Agenv7(X,Y,Z,repmat(Dgparams,NPixel/size(Dgparams,1),1),NPixel,DProbe);
D = D./max(max(D)); 




% (later TODO): need to design a selection procedure of combinations of LEDs
%               the selected combinations of Emitters in Ek needs to be orthogonal
%               otherwise, DISTINCT will fail
%               We know Ek is orthogonal if ETE=Ek'*Ek = diagonal matrix
%               An orthogonal Ek will give a well-defined gradient which makes it convex
%               And also ensures correctness of the analytical solution in solvediagIntProd.m
%              
%
%               It turns out, Ek doesnt need to be perfectly orthogonal 
%               So ETE doesnt need to be perfectly diagonal
%               But the question is, how orthogonal? 
%               How do we quantify this?

%          a good selection procedure (Ek = combinations*E):
        %minimizes the cross-correlation between rows in Ek
        %maximizes the amount of LEDs on at the same time for higher light input thus SNR
        %minimizes the # of columns in Ek -> each column in Ek is 1 recording
        %the less recordings needed -> higher framerate
        

if 0 %adjustment of emitter matrix into clusters of LEDs on and off
    
    %it turns out, this snippet of code generates an Ek which is not orthogonal enough
    combinations = repmat(emitter_Patterns(14,clustersize),1,NPixel/clustersize);    
   
    
    %so, how do we design an orthogonal combination of LEDs which gives the same reconstruction quality, with less rows in Ek?    
    %     combinations = tristans_combination_procedure(E);
    
    
    Ek = combinations*E; Ek = Ek./max(max(Ek)); 
    
else %else here, each row of E is simply one single pixel enabled
     %that means 1024 measurements for a single image, cutting effective framerate by 1024
     %while each measurement only records from a very small portion of the volume of interest
     %you can imagine there is a lot of room for optimization here
    
    Ek = E; %but... E is orthogonal enough to make DISTINCT converge exactly.
end

K = size(Ek,1);

%generate the derived matrices
Dpinv = pinv(D);
Ekpinv = pinv(Ek);
projD = Dpinv*D;        %detector projection matrix
projEk = Ek'*Ekpinv';   %emitter projection matrix
DdotE = D.*sum(Ek);     %single pattern, all LEDs on, doesnt use any information in turning on different patterns
% DkronE = kron(D,Ek);  %warning: this kronecker matrix will be massive
ETE = Ek'*Ek;
DTD = D'*D;


%because each emitter profile is almost orthogonal to all the others, 
%ETE = E'*E is almost a diagonal matrix as you can see here:
figure; imagesc(E'*E); title('All emmitters separately')
figure; imagesc(Ek'*Ek); title('Adriaans terrible selection procedure')
%same for DTD:
figure; imagesc(DTD); title('Detector matrix inner product')


%% STEP 3 give parameters for monte carlo
%data parameters
L=[4];          %number of sources    
SNR = 100;      %for implementation of noise
nIter = 1000;   %number of repetitions in the monte carlo

%parameters for the solvers
maxIter = 10;
maxIterInt = 1e4;
addSize = 5;
intSolver = 0;
lambdadiv = 10;

%parameters for SL0
mu_0 = 2;
q = 3;
sigma_min = 2/SNR; %As a rule of tumb, for the noisy case, sigma_min should be about 2 to 4 times of the standard deviation of this noise
sigma_decrease_factor = 0.5;
max_contrast_ratio = 1e3;

for ll = 1:length(L) 
    for ii = 1:nIter
        %% STEP 4 put down sources, calculate the data 
        [~,~,XYZtrue,~,~,xtrue] = create_randomRGlocs2D(X,Y,Z,L(ll),pi/2,X(1),Y(1));
        figure; scatter3(XYZtrue(:,1),XYZtrue(:,2),XYZtrue(:,3)),'k';
        axis([X(1) X(end) Y(1) Y(end) Z(1) Z(end)])

        %generate the compressed data as seen by the chip
        y = D*diag(xtrue)*Ek';
%         noise=imnoise(y./SNR,'poisson');
        noise = 0;
        y = y+noise;
        ydotE = DdotE*xtrue+noise;

        %In the noiseless case, we can use SVD to guess number of sources!
    %     [U,S,V] = svd(y,'econ');
    %     figure; semilogy(diag(S),'x')

        %% STEP 5 reconstruct using each method
        %a. DISTINCT: 
            gradconstdiag = diag(D'*y*Ek);
            lambda = max(gradconstdiag)/lambdadiv;
            tic; [xhatDS,f1,f2,f1f2] = DISTINCT(y,D,Ek,lambda,ETE,DTD,gradconstdiag,intSolver,maxIter,maxIterInt,addSize);
            tcompute(1,ll,ii) = toc;
            disp(['Found ' num2str(numel(find(xhatDS))) ' locations'])
            %choose L largest entries -> for exact comparison wth ground truth
            [~,DSind] = sort(xhatDS,'descend'); xhatDS(DSind(L(ll)+1:end)) = 0;


        %b. ICopt, no extra emitter rows, faster but uses less information
            %so its expected to be far less accurate
            lambda = max(ydotE'*DdotE)/lambdadiv;
            tic; [xhatIC, ~] = ICOpt(DdotE,ydotE,lambda,0,1e10,0,1,addSize);
            tcompute(2,ll,ii) = toc;    
            disp(['Found ' num2str(numel(find(xhatIC))) ' locations'])
            %choose L largest entries -> for exact comparison wth ground truth
            [~,ICind] = sort(xhatIC,'descend'); xhatIC(ICind(L(ll)+1:end)) = 0;


        %c. 2D SL0, adapted for 2D. Already improved speedup ~2-3x 
            %over the original authors up but still slow as balls
            %the second speedup round -> speeds up by another 20x
            
            sl0gradconstdiag = diag(Dpinv*y*Ekpinv');
            tic; 
%             xhatSL0=SL0_2D_ATv2(y, sigma_min, sigma_decrease_factor, mu_0, q, Dpinv, Ekpinv, projD, projEk, sl0gradconstdiag); xhatSL0 = wthresh(xhatSL0,'h',max(xhatSL0)/max_contrast_ratio);            
            tic; xhatSL0=SL0_2D_ATv3(y, sigma_min, sigma_decrease_factor, mu_0, q, Dpinv, Ekpinv, projD, projEk, sl0gradconstdiag, max_contrast_ratio);
            tcompute(3,ll,ii) = toc;    
            %choose L largest entries -> for exact comparison wth ground truth
            [~,SL0ind] = sort(xhatSL0,'descend'); xhatSL0(SL0ind(L(ll)+1:end)) = 0;   
            
            
        %d. competitor #3

        
        %e. competitor #4 -> our own algorithm using the 
        %gradient descent internal solver, for comparison
            lambda = max(gradconstdiag)/lambdadiv;
            tic; [xhatDSGD,f1,f2,f1f2] = DISTINCT(y,D,Ek,lambda,ETE,DTD,gradconstdiag,1,maxIter,maxIterInt,addSize);
            tcompute(5,ll,ii) = toc;
            disp(['Found ' num2str(numel(find(xhatDSGD))) ' locations'])
            %choose L maximum -> for exact comparison wth ground truth
            [~,DSGDind] = sort(xhatDSGD,'descend'); xhatDSGD(DSGDind(L(ll)+1:end)) = 0;



        %TODO: 
        %% STEP 6 calculate the accuracy of each solution compared to ground truth
            % -> find a measure of accuracy. 
            %      Simply calculating the least euclidian distance between the 
            %      found locations and the true locations is an NP-hard problem. 
            %      For more than 10 sources it becomes computationally infeasible.
            % -> How do other papers define accuracy?        


            %in terms of the sequence x
        accuracy(1,ll,ii) = some_accuracy_function(xhatDS,xtrue);
        accuracy(2,ll,ii) = some_accuracy_function(xhatIC,xtrue);
        accuracy(3,ll,ii) = some_accuracy_function(xhatSL0,xtrue);
%         accuracy(4,ll,ii) = some_accuracy_function(,xtrue)
        accuracy(5,ll,ii) = some_accuracy_function(xhatDSGD,xtrue);


            %or in terms of the geometric locations
%         [xhatDSind(:,1),xhatDSind(:,2),xhatDSind(:,3)] = ind2sub(sizes,find(xhatDS));
%         [xhatICind(:,1),xhatICind(:,2),xhatICind(:,3)] = ind2sub(sizes,find(xhatIC));
%         accuracy(ii,1) = another_accuracy_function(xhatDSind,XYZloc)
%         accuracy(ii,2) = another_accuracy_function(xhatICind,XYZloc)


        %TODO: STEP 7: measure the RAM usage of the system
            % Consult the sysadmin on slack "computing channel" how to do this
            % on the linux servers 
            % for each algorithm separately while they are running

    end
end

figure
stem(sum(accuracy,3))
xlabel('Algorithm type'); ylabel('Normalized accuracy (??)');
title('Average accuracy');
xlim([0 6])

figure
stem(sum(tcompute,3));
xlabel('Algorithm type'); ylabel('Time (seconds)');
title('Time to compute');
xlim([0 6])

% plot the results
xtrueplot = xtrue; xtrueplot(xtrueplot==0) = NaN;
xhatICplot = xhatIC; xhatICplot(xhatICplot==0) = NaN;
xhatDSplot = xhatDS; xhatDSplot(xhatDSplot==0) = NaN;
xhatSL0plot = xhatSL0; xhatSL0plot(xhatSL0plot==0) = NaN;
xhatDSGDplot = xhatDSGD; xhatDSGDplot(xhatDSGDplot==0) = NaN;

figure; hold on; set(gcf,'Color',[1 1 1]);
stem(xtrueplot,'k')
stem(xhatDSplot./max(xhatDSplot),'m')
stem(xhatICplot./max(xhatICplot),'c')
stem(xhatSL0plot./max(xhatSL0plot),'g')
stem(xhatDSGDplot./max(xhatDSGDplot),'r')
axis([0 prod(sizes) -0.2 1.2])
legend('True locations','DISTINCT','ICopt','SL0','DISTINCT GD')


% scatterplot
clear xhatind xhatICind xhatDSind xhatSL0ind
[xhatDSind(:,1),xhatDSind(:,2),xhatDSind(:,3)] = ind2sub(sizes,find(xhatDS));
[xhatICind(:,1),xhatICind(:,2),xhatICind(:,3)] = ind2sub(sizes,find(xhatIC));
[xhatSL0ind(:,1),xhatSL0ind(:,2),xhatSL0ind(:,3)] = ind2sub(sizes,find(xhatSL0));
[xhatDSGDind(:,1),xhatDSGDind(:,2),xhatDSGDind(:,3)] = ind2sub(sizes,find(xhatDSGD));

%plot the results
figure; hold on; set(gcf,'Color',[1 1 1]);
scatter3(XYZtrue(:,1),XYZtrue(:,2),XYZtrue(:,3)),'k';
scatter3(X(xhatDSind(:,1)),Y(xhatDSind(:,2)),Z(xhatDSind(:,3)),'m*')
scatter3(X(xhatICind(:,1)),Y(xhatICind(:,2)),Z(xhatICind(:,3)),'c^')
scatter3(X(xhatSL0ind(:,1)),Y(xhatSL0ind(:,2)),Z(xhatSL0ind(:,3)),'gd')
scatter3(X(xhatDSGDind(:,1)),Y(xhatDSGDind(:,2)),Z(xhatDSGDind(:,3)),'ro')
axis([X(1) X(end) Y(1) Y(end) Z(1) Z(end)])
xlabel('X [um]'); ylabel('Y [um]'); zlabel('Z [um]');
% set(gca,'PlotBoxAspectRatio',[1 ARy ARz])
legend('True locations','DISTINCT','ICopt','SL0','DISTINCT GD'); view([-45 45]); grid on;



%% make pictures of the matrices for geometrical interpretation
Ecube = squeeze(permute(reshape(E,NPixel,length(X),length(Y),length(Z)),[1 3 2 4]));
Dcube = squeeze(permute(reshape(D,NPixel,length(X),length(Y),length(Z)),[1 3 2 4]));

%normalize each slice so the plotting becomes visible
Ecubenorm = bsxfun(@rdivide,Ecube,max(max(max(Ecube))));
Dcubenorm = bsxfun(@rdivide,Dcube,max(max(max(Dcube))));

ARy = (Y(end)-Y(1))/(X(end)-X(1));
ARz = (Z(end)-Z(1))/(X(end)-X(1));
[Xm,Ym,Zm] = meshgrid(X,Y,Z);

%image the results
pixplot = 178;
for zz = length(Z)
    figure; set(gcf,'Color',[1 1 1]);
    subplot(2,1,1)
    imagesc(X,Y,squeeze(Ecube(pixplot,:,:,zz))); 
    title('Emitter field')
    set(gca,'PlotBoxAspectRatio',[1 ARy ARz])
    xlabel('X [um]'); ylabel('Y [um]');
    set(gca,'YDir','normal')

    subplot(2,1,2)
    imagesc(X,Y,squeeze(Dcube(pixplot,:,:,zz)));
    title('Detector field')
    set(gca,'PlotBoxAspectRatio',[1 ARy ARz])
    xlabel('X [um]'); ylabel('Y [um]');
    set(gca,'YDir','normal')
end
axis([X(1) X(end) Y(1) Y(end)])



%emitter slice plot
xslice = []; 
yslice = [];
zslice = Z;
figure; set(gcf,'Color',[1 1 1])
slice(Xm,Ym,Zm,squeeze(Ecubenorm(pixplot,:,:,:)),xslice,yslice,zslice)
colormap(parula);
axis([X(1) X(end) Y(1) Y(end) 0 Z(end)])
set(gca,'PlotBoxAspectRatio',[1 ARy 1])
shading interp; light;
xlabel('X [um]'); ylabel('Y [um]'); zlabel('Z [um]');
title('Volumetric Emitter field, normalized per slice')
view([-30 30])


%detector slice plot
xslice = []; 
yslice = [];
zslice = Z;
figure; set(gcf,'Color',[1 1 1])
slice(Xm,Ym,Zm,squeeze(Dcubenorm(pixplot,:,:,:)),xslice,yslice,zslice)
colormap(parula);
axis([X(1) X(end) Y(1) Y(end) 0 Z(end)])
set(gca,'PlotBoxAspectRatio',[1 ARy 1])
shading interp; light;
xlabel('X [um]'); ylabel('Y [um]'); zlabel('Z [um]');
title('Volumetric detector field, normalized per slice')
view([-30 30])


