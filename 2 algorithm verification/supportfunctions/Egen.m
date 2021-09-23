function [ E ] = Egen( X,Y,Z,Eparams,NPixel,EProbe)
%Creates an E-matrix for illumination 
%   params = [sigma_main mmain elevation sigma_elevation azimuth sigma_azimuth mside]

for nn = 1:NPixel
    for zz=1:length(Z)
                
        %meshgrids
        [Xm,Ym] = meshgrid(X-EProbe(nn,1),Y-EProbe(nn,2));
        
        [~,elM,rM] = cart2sph(Xm(:),Ym(:),Z(zz));
        elM = reshape(pi/2 - elM,length(Y),length(X));        
        
        %load parameters
        sigma_main = Eparams(nn,1);
        mmain = Eparams(nn,2);
        elevation = Eparams(nn,3);
        sigma_elevation = Eparams(nn,4);        
        azimuth = Eparams(nn,5);
        sigma_azimuth = Eparams(nn,6);
        mside = Eparams(nn,7);
        deadangle = Eparams(nn,8);
        
        %create main lobe        
        sigmain = Z(zz)*tan([sigma_main 0; 0 sigma_main]);
        Emain = mvnpdf([Xm(:) Ym(:)],[0,0],sigmain);
        
        %normalize by sum to equalize the power in each slice
        e = reshape(1+mmain*Emain/sum(sum(Emain)),length(Y),length(X));
        
        %create four side lobes
        for kk = 1:4
            az = azimuth + (kk-1)*pi/2;
            eigs = [cos(az) sin(az);sin(az) -cos(az)];
            sigside = Z(zz)*tan([sigma_elevation 0; 0 sigma_azimuth]);
            sigside = eigs'*sigside*eigs;

            muside =  Z(zz)*tan([elevation*cos(az) elevation*sin(az)]);

            Eside = mvnpdf([Ym(:) Xm(:)],muside,sigside);
            
            %normalize and reshape
            Eside = mside*reshape(Eside/sum(sum(Eside)),length(Y),length(X));
            e = e + Eside;
        end 
        

        %we have 1 lambertian windowing function
%         window = 1;
        window = cos(elM).^2;
        
        %And a scalar for distance
%         scaler = 1;
        scaler = 1./(4*pi*sqrt(rM)); 
        scaler = reshape(scaler,length(Y),length(X));
        
        %And a hard cutoff deadangle window      
        deadangle_gain = ones(length(Y),length(X));                
        deadangle_gain(find(abs(elM)>deadangle)) = 0;

        %add the windowing and scaling
        output(:,:,zz) = (e.*window.*scaler.*deadangle_gain)';
    end
    E(nn,:) = output(:);
end
E(E<0) = 0;
end