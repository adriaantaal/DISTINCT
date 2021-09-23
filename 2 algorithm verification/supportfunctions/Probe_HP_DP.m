function [hotdeadpix] = Probe_HP_DP(N_DCR, thr, mode)
 %input
    %N_Probe is the 2D full measurement matrix [Npixel x Nframe]
    %N_DC is the dark count measurement matrix
    %mode is 0 for thr*mean, 1 for thr*median
    %thr is this scaler factor 

 %output
    %indices of hot and dead pixels
NPixel = size(N_DCR,1);
deadpix = [];

if nargin < 2
   thr = 10;
end

if nargin < 3
   mode = 0;
end

%use vectors
if size(N_DCR,2) > 1 
    N_DCR = sum(N_DCR,2);
end


%identify hotpixels from N_DC matrix
if (mode == 0)
    Nmean = mean(N_DCR);
    hotpix = find(N_DCR>(thr*Nmean));
    %identify deadpixels from the data matrix
%     deadpix = find(N_DCR==0);
    disp(['Mean DCR is ' num2str(Nmean)])
elseif (mode == 1)
    Nmed = median(N_DCR);
    hotpix = find(N_DCR>(thr*Nmed));
    %identify deadpixels from the data matrix
%     deadpix = find(N_DCR==0);
    disp(['Median DCR is ' num2str(Nmed)])
end

%output
hotdeadpix = unique([hotpix; deadpix]);
disp(['Percentage of hotpix is ' num2str(round(length(hotpix)/NPixel,3)*100) '%' ])
disp(['Percentage of deadpix is ' num2str(round(length(deadpix)/NPixel,3)*100) '%' ])

end