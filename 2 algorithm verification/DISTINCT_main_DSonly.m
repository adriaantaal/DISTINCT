%% Welcome to the code that simulates DISTINCT 
%First create a model for Detector and Emitter array
%Then apply the equations to verify reconstruction

if ~exist('numBatchFile','var')
    t = datetime('now');
    mmdd = datestr(t);
    numBatchFile = str2double(mmdd(end));
end

rng(numBatchFile);


%create directory
t = datetime('now');
mmdd = datestr(t);
mkdir(['myoutput_' mmdd]);


%% STEP 0: import all functions and folders
Folder = fullfile(cd, '..');
addpath(genpath(fullfile(Folder)));

%add the whole folder on MAC here:
addpath(genpath('D:\Dropbox\96 Paper DISTINCT\2 algorithm verification'))

addpath(genpath('/proj2/adriaantaal/96 Paper DISTINCT'))


%% STEP 3 give parameters for monte carlo
%data parameters
L=[2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32 64 128 256 512 1024];          %number of sources    
% L=[2 4 6];          %number of sources    
SNR = 100;      %for implementation of noise
nIter = 100;   %number of repetitions in the monte carlo

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
res = [70 50 40 30 25 20];

for rr = 1:length(res)
    %% STEP 1: make a voxel grid and probe layout
    %create probe layout
    NPixel = 1024;
    xpitch = 24.6;
    ypitchD = 95.96;
    ypitchE = 31;
    DProbe = create_Probe(NPixel,xpitch,ypitchD,259.4,1);
    EProbe = create_Probe(NPixel,xpitch,ypitchE,259.4,1);
    EProbe(:,2) = EProbe(:,2) + (ypitchD-ypitchE)/2; %shift emitters to center of shank
    
    xres = res(rr);
    yres = res(rr);
    zres = res(rr);
    X = (-250+min(DProbe(:,1))):xres:(max(DProbe(:,1))+250);
    Y = (-250+min(DProbe(:,2))):yres:(max(DProbe(:,2))+250);
    Z = 200:zres:300;
    sizes = [length(X) length(Y) length(Z)]; N(rr) = prod(sizes);
    
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

    Ek = E; %but... E is orthogonal enough to make DISTINCT converge exactly.
    K = size(Ek,1);

    %generate the derived matrices
%     Dpinv = pinv(D);
%     Ekpinv = pinv(Ek);
%     projD = Dpinv*D;        %detector projection matrix
%     projEk = Ek'*Ekpinv';   %emitter projection matrix
%     DdotEk = D.*sum(Ek);     %single pattern, all LEDs on, doesnt use any information in turning on different patterns
    ETE = Ek'*Ek;
    DTD = D'*D;
    H = DTD.*ETE;
    Dpinv = pinv(D);
    Epinv = pinv(Ek);

    
    for ll = 1:length(L) 
        for ii = 1:nIter
            %% STEP 4 put down sources, calculate the data 
            [~,~,XYZtrue,~,~,xtrue] = create_randomRGlocs2D(X,Y,Z,L(ll),pi/2,X(1),Y(1));
            %2D detector-emitter formulation
            y = D*diag(xtrue)*Ek';
            noise = 0; y = y+noise;
            gradconstdiag = diag(D'*y*Ek);
            yf2 = norm(y,'f')^2;
%             sl0gradconstdiag = diag(Dpinv*y*Ekpinv');
            
            %1D formulation
%             ydotEk = DdotEk*xtrue+noise;

            %derived parameters
            addSize = max(floor(L(ll)/2),1);

            lambda = max(gradconstdiag)/lambdadiv;
%             lambdadotEk = max(ydotEk'*DdotEk)/lambdadiv;
            
            %In the noiseless case, we can use SVD to guess number of sources!
            [U,S,V] = svd(y,'econ'); 
            
            %remove singular values based on value, faster, less accurate
            U = U(:,1:L(ll)); V = V(:,1:L(ll)); S = S(1:L(ll),1:L(ll));
            s = diag(S);
            
            %scale the gradconstdiag and singular values by lambda
            gradconstdiag = gradconstdiag-lambda;
%             s2 = s - lambda*sum((U'*Dpinv').*(Epinv*V)',2); %no good
%             s2 = s - 0.5*lambda*sum(pinv(D'*U).*pinv(V'*E)',2); %maybe good
            
            %need to remove negative entries
%             U2 = U(:,(s2>0)); V2 = V(:,(s2>0)); s2 = s2(s2>0); 
%             U2 = U; V2 = V;

            %% STEP 5 reconstruct using each method
            %a. DISTINCT: 
                tic; [xhatDSV6,~,~,~] = DISTINCTv6(y,D,Ek,lambda,H,gradconstdiag,U,s,V,0,maxIter,maxIterInt,addSize);
                tcompute(1,rr,ll,ii) = toc;
                disp(['Found ' num2str(numel(find(xhatDSV6))) ' locations'])
                %choose L largest entries -> for exact comparison wth ground truth
                [~,DSV6ind] = sort(xhatDSV6,'descend'); xhatDSV6(DSV6ind(L(ll)+1:end)) = 0;

            %d. normal equations using thresholded gradconstdiag
                tic; [xhatDGV6,~,~,~] = DISTINCTv6(y,D,Ek,lambda,H,gradconstdiag,U,s,V,1,maxIter,maxIterInt,addSize);
                tcompute(2,rr,ll,ii) = toc;
                disp(['Found ' num2str(numel(find(xhatDGV6))) ' locations'])
                %choose L maximum -> for exact comparison wth ground truth
                [~,DGV6ind] = sort(xhatDGV6,'descend'); xhatDGV6(DGV6ind(L(ll)+1:end)) = 0;
            
            %a. DISTINCT: 
                tic; [xhatDSV7,j7s] = DISTINCTv7(yf2,D,Ek,lambda,H,gradconstdiag,U,s,V,0,maxIter,maxIterInt,addSize);
                tcompute(3,rr,ll,ii) = toc;
                disp(['Found ' num2str(numel(find(xhatDSV7))) ' locations'])
                %choose L largest entries -> for exact comparison wth ground truth
                [~,DSV7ind] = sort(xhatDSV7,'descend'); xhatDSV7(DSV7ind(L(ll)+1:end)) = 0;
                
            %d. normal equations using thresholded gradconstdiag
                tic; [xhatDGV7,j7g] = DISTINCTv7(yf2,D,Ek,lambda,H,gradconstdiag,U,s,V,1,maxIter,maxIterInt,addSize);
                tcompute(4,rr,ll,ii) = toc;
                disp(['Found ' num2str(numel(find(xhatDGV7))) ' locations'])
                %choose L maximum -> for exact comparison wth ground truth
                [~,DGV7ind] = sort(xhatDGV7,'descend'); xhatDGV7(DGV7ind(L(ll)+1:end)) = 0;
           
                
            %% STEP 6 calculate the accuracy of each solution compared to ground truth
                % -> find a measure of accuracy. 
                %      Simply calculating the least euclidian distance between the 
                %      found locations and the true locations is an NP-hard problem. 
                %      For more than 10 sources it becomes computationally infeasible.
                % -> How do other papers define accuracy?        

                %in terms of the sequence x
            accuracy(1,rr,ll,ii) = some_accuracy_function(xhatDSV6,xtrue);
            accuracy(2,rr,ll,ii) = some_accuracy_function(xhatDGV6,xtrue);
            accuracy(3,rr,ll,ii) = some_accuracy_function(xhatDSV7,xtrue);
            accuracy(4,rr,ll,ii) = some_accuracy_function(xhatDGV7,xtrue);
            
            %TODO: STEP 7: measure the RAM usage of the system
                % Consult the sysadmin on slack "computing channel" how to do this
                % on the linux servers 
                % for each algorithm separately while they are running            
        end
    end
    save(['myoutput_' mmdd '/DISTINCT_sim_DSv7only.mat'],'accuracy','tcompute','L','SNR','nIter','addSize','lambdadiv','q','sigma_min','sigma_decrease_factor','max_contrast_ratio','X','Y','Z','numBatchFile','K','N');
end

accmean = nanmean(accuracy,4);
accstd = std(accuracy,0,3,'omitnan');

tcomputemean = nanmean(tcompute,4);
tcomputestd = std(tcompute,0,3,'omitnan');

%create figure accuracy
f1 = figure;
plot(L, accmean)
xlim([2 1024]); ylim([0 1]); set(gca,'xscale','log');
xlabel('Number of Sources (L)'); ylabel('Accuracy');
legend('DS','DSGD')
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
xlim([2 1024]); set(gca,'xscale','log');
xlabel('Number of Sources (L)'); ylabel('Accuracy');
legend('DS','DSGD')
savefig(f2,['myoutput_' mmdd '/AccuracyGraphWithError']) ;


%create figure computing time
f4 = figure;
plot(L, tcomputemean)
xlim([2 1024]); set(gca,'xscale','log'); set(gca,'yscale','log');
xlabel('Number of Sources (L)'); ylabel('Computing Time(s)');
legend('DS','DSGD')
savefig(f4,['myoutput_' mmdd '/ComputingTime']) ;

figure; set(gcf,'Color',[1 1 1]);
subplot(3,1,1); plot(f1); title('J1 DISTINCT')
subplot(3,1,2); plot(f2); title('J2 DISTINCT') 
subplot(3,1,3); plot((f1+f2)); title('J1+J2 DISTINCT')
