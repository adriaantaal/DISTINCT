%% Welcome to the code that simulates DISTINCT 
%First create a model for Detector and Emitter array
%Then apply the equations to verify reconstruction

% if ~exist('numBatchFile','var')
%     t = datetime('now');
%     mmdd = datestr(t);
%     rng(mmdd(end));
% else
    rng(numBatchFile);
% end


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

xres = 25;
yres = 25;
zres = 25;
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
        
    
%else here, each row of E is simply one single pixel enabled
%that means 1024 measurements for a single image, cutting effective framerate by 1024
%while each measurement only records from a very small portion of the volume of interest
%you can imagine there is a lot of room for optimization here
    
Ek = E; %but... E is orthogonal enough to make DISTINCT converge exactly.

K = size(Ek,1);

%generate the derived matrices
Dpinv = pinv(D);
Ekpinv = pinv(Ek);
projD = Dpinv*D;        %detector projection matrix
projEk = Ek'*Ekpinv';   %emitter projection matrix
DdotE = D.*sum(Ek);     %single pattern, all LEDs on, doesnt use any information in turning on different patterns
ETE = Ek'*Ek;
DTD = D'*D;
H = DTD.*ETE;

%because each emitter profile is almost orthogonal to all the others, 
%ETE = E'*E is almost a diagonal matrix as you can see here:
% figure; imagesc(ETE); title('Emmitter matrix inner product')
%same for DTD:
% figure; imagesc(DTD); title('Detector matrix inner product')


%% STEP 3 give parameters for monte carlo
%data parameters
L=[2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32];          %number of sources    
% L=[2 4 6];          %number of sources    
SNR = 100;      %for implementation of noise
nIter = 10;   %number of repetitions in the monte carlo

%parameters for the solvers
maxIter = 100;
maxIterInt = 1e4;
lambdadiv = 10;

%parameters for SL0
mu_0 = 2;
q = 3;
sigma_min = 0.02;               %As a rule of tumb, for the noisy case, sigma_min should be about 2 to 4 times of the standard deviation of this noise
sigma_decrease_factor = 0.5;
max_contrast_ratio = 1e3;       %higher max_constrast_ratio increases accuracy but scales computing time
% for ss = 1:length(SNR)
for ll = 1:length(L) 
    for ii = 1:nIter
        %% STEP 4 put down sources, calculate the data 
        [~,~,XYZtrue,~,~,xtrue] = create_randomRGlocs2D(X,Y,Z,L(ll),pi/2,X(1),Y(1));
        
        %2D detector-emitter formulation
        y = D*diag(xtrue)*Ek';
        noise = 0; y = y+noise;
        gradconstdiag = diag(D'*y*Ek);
        lambda = max(gradconstdiag)/lambdadiv;
        sl0gradconstdiag = diag(Dpinv*y*Ekpinv');

        %1D formulation
        ydotE = DdotE*xtrue+noise;
        lambdadotEk = max(ydotE'*DdotE)/lambdadiv;
        
        %derived parameters
        addSize = max(floor(L(ll)/2),1);
        
        %In the noiseless case, we can use SVD to guess number of sources!
    %     [U,S,V] = svd(y,'econ');
    %     figure; semilogy(diag(S),'x')

    
        %% STEP 5 reconstruct using each method
        %a. DISTINCT: 
            tic; [xhatDS,~,~,~] = DISTINCTv2(y,D,Ek,lambda,H,gradconstdiag,0,maxIter,maxIterInt,addSize);
            tcompute(1,ll,ii) = toc;
            disp(['Found ' num2str(numel(find(xhatDS))) ' locations'])
            %choose L largest entries -> for exact comparison wth ground truth
            [~,DSind] = sort(xhatDS,'descend'); xhatDS(DSind(L(ll)+1:end)) = 0;


        %b. ICopt, no extra emitter rows, faster but uses less information
            %so its expected to be far less accurate
            tic; [xhatIC, ~] = ICOpt(DdotE,ydotE,lambdadotEk,0,Inf,0,1,addSize);
            tcompute(2,ll,ii) = toc;    
            disp(['Found ' num2str(numel(find(xhatIC))) ' locations'])
            %choose L largest entries -> for exact comparison wth ground truth
            [~,ICind] = sort(xhatIC,'descend'); xhatIC(ICind(L(ll)+1:end)) = 0;


        %c. 2D SL0, adapted for 2D. Already improved speedup ~2-3x 
            %over the original authors up but still slow as balls
            %the second speedup round -> speeds up by another 20x
            tic; xhatSL0=SL0_2D_ATv3(y, sigma_min, sigma_decrease_factor, mu_0, q, Dpinv, Ekpinv, projD, projEk, sl0gradconstdiag, max_contrast_ratio);
            tcompute(3,ll,ii) = toc;    
            %choose L largest entries -> for exact comparison wth ground truth
            [~,SL0ind] = sort(xhatSL0,'descend'); xhatSL0(SL0ind(L(ll)+1:end)) = 0;   
            
            
%         %d. competitor #3: 2D NIHT algorithm
%         %Slow because it relies on calculation of many matrix norms
%             tic; [xhatNIHT] = NIHT2D_ATv4(y,D,Ek,L(ll));
%             tcompute(4,ll,ii) = toc;    
%             %choose L largest entries -> for exact comparison wth ground truth
%             [~,NIHTind] = sort(xhatNIHT,'descend'); xhatNIHT(NIHTind(L(ll)+1:end)) = 0;   
            
                    
        %e. competitor #4 -> our own algorithm using the 
        %gradient descent internal solver, for comparison
            tic; [xhatDSGD,~,~,~] = DISTINCTv2(y,D,Ek,lambda,H,gradconstdiag,1,maxIter,maxIterInt,addSize);
            tcompute(5,ll,ii) = toc;
            disp(['Found ' num2str(numel(find(xhatDSGD))) ' locations'])
            %choose L maximum -> for exact comparison wth ground truth
            [~,DSGDind] = sort(xhatDSGD,'descend'); xhatDSGD(DSGDind(L(ll)+1:end)) = 0;

        %f. competitor #5 -> our own algorithm using the 
        %direct inversion method
            tic; [xhatDINV,~,~,~] = DISTINCTv2(y,D,Ek,lambda,H,gradconstdiag,2,maxIter,maxIterInt,addSize);
            tcompute(6,ll,ii) = toc;
            disp(['Found ' num2str(numel(find(xhatDINV))) ' locations'])
            %choose L maximum -> for exact comparison wth ground truth
            [~,DINVind] = sort(xhatDINV,'descend'); xhatDINV(DINVind(L(ll)+1:end)) = 0;
    
        %g. DISTINCT using linsolve: 
            tic; [xhatDSQP,f1,f2,f1f2] = DISTINCTv2(y,D,Ek,lambda,H,gradconstdiag,3,maxIter,maxIterInt,addSize);
            tcompute(7,ll,ii) = toc;
            disp(['Found ' num2str(numel(find(xhatDSQP))) ' locations'])
            %choose L largest entries -> for exact comparison wth ground truth
            [~,DSQPind] = sort(xhatDSQP,'descend'); xhatDSQP(DSQPind(L(ll)+1:end)) = 0;
        

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
        accuracy(4,ll,ii) = 0;
        accuracy(5,ll,ii) = some_accuracy_function(xhatDSGD,xtrue);
        accuracy(6,ll,ii) = some_accuracy_function(xhatDINV,xtrue);
        accuracy(7,ll,ii) = some_accuracy_function(xhatDSQP,xtrue);
        
        %TODO: STEP 7: measure the RAM usage of the system
            % Consult the sysadmin on slack "computing channel" how to do this
            % on the linux servers 
            % for each algorithm separately while they are running            
    end
end

%create directory
t = datetime('now');
mmdd = datestr(t);
mkdir(['myoutput_' mmdd]);

save(['myoutput_' mmdd '/DISTINCT_sim_large.mat'],'accuracy','tcompute','L','SNR','nIter','addSize','lambdadiv','q','sigma_min','sigma_decrease_factor','max_contrast_ratio','X','Y','Z','numBatchFile');


accmean = nanmean(accuracy,3);
accstd = std(accuracy,0,3,'omitnan');

tcomputemean = nanmean(tcompute,3);
tcomputestd = std(tcompute,0,3,'omitnan');




%create figure accuracy
f1 = figure;
plot(L, accmean)
xlim([0 1024]); ylim([0 1]); set(gca,'xscale','log');
xlabel("Number of Sources (L)"); ylabel("Accuracy");
legend('DS','ICOpt','SL0','NIHT','DSGD','DINV')
savefig(f1,['myoutput_' mmdd '/AccuracyGraph']) ;


f2 = figure;
errorbar(L,accmean(1,:),accstd(1,:))
hold on;
errorbar(L,accmean(2,:),accstd(2,:))
errorbar(L,accmean(3,:),accstd(3,:))
errorbar(L,accmean(4,:),accstd(4,:))
errorbar(L,accmean(5,:),accstd(5,:))
errorbar(L,accmean(6,:),accstd(6,:))
hold off;
xlim([0 1024]); set(gca,'xscale','log');
xlabel("Number of Sources (L)"); ylabel("Accuracy");
legend('DS','ICOpt','SL0','NIHT','DSGD','DINV')
savefig(f2,['myoutput_' mmdd '/AccuracyGraphWithError']) ;


%create figure computing time
f4 = figure;
plot(L, mean(tcompute,3))
xlim([0 1024]); set(gca,'xscale','log');
xlabel("Number of Sources (L)"); ylabel("Computing Time(s)");
legend('DS','ICOpt','SL0','NIHT','DSGD','DINV')
savefig(f4,['myoutput_' mmdd '/ComputingTime']) ;

