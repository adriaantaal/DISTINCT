function [ A ] = Agenv7(X,Y,Z,params,NPixel,DProbe)
%Create a model for ASP
%   params = [alpha beta modindex chi gammapara mpara gammaperp mperp deadangle]

for nn = 1:NPixel
    for zz=1:length(Z)        
        %load params
        my_alpha = params(nn,1);
        my_beta = params(nn,2);
        m = params(nn,3);
        chi = params(nn,4);        

        if size(params,2) > 4
            gammapara = params(nn,5);
            mpara = params(nn,6);
            gammaperp = params(nn,7);
            mperp = params(nn,8);
            deadangle = params(nn,9);
        else
            gammapara = 3/2;
            mpara = 1;
            gammaperp = 0;
            mperp = 0;
            deadangle = 90/180*pi;
        end
        
        %center the X-Y to the pixel
        X0 = (X-DProbe(nn,1));
        Y0 = (Y-DProbe(nn,2));

        %meshgrids
        [Xm,Ym] = meshgrid(X0,Y0);
        
        %rotate the meshgrids with the grating orientation angle chi
        rotmat = [cos(chi) -sin(chi); sin(chi) cos(chi)];
        XYrot = rotmat*[Xm(:)';Ym(:)'];
        
        %reshape
        Xmrot = reshape(XYrot(1,:),length(Y),length(X));
        Ymrot = reshape(XYrot(2,:),length(Y),length(X));
        
        %translate to angle
        Tvr = atan2(Z(zz),Xmrot)-pi/2;
        Pvr = atan2(Z(zz),Ymrot)-pi/2;   
        
        %translate to spherical coordinates        
        [~,elM,rM] = cart2sph(Xm(:),Ym(:),Z(zz));
        elM = reshape(pi/2 - elM,length(Y),length(X));
        rM  = reshape(rM,length(Y),length(X));
                
        %we have 2 windowing functions, 1 parallel and 1 perpendicular
        windowpara = 1+mpara*cos(gammapara*Tvr);

        %perpendicular window 
        windowperp = 1+0.5*mperp*cos(gammaperp*Pvr+pi);
        
        %a global lambertian windowing function
        windowing = cos(elM).^2;
        
        %And a scalar for distance
%         scaler = 1; 
%         scaler = 1./(4*pi*sqrt(rM)); 
%         scaler = 1./(4*pi*rM.^2);
        scaler = 1./(4*pi*rM);
        
        %And a hard cutoff deadangle window      
        deadangle_gain = ones(length(Y),length(X));                
        deadangle_gain(find(abs(elM)>deadangle)) = 0;

        %calculate the response
        a =(1+m*cos(Tvr*my_beta+my_alpha));          

        %add the windowing
        output(:,:,zz) = (a.*windowpara.*windowperp.*windowing.*scaler.*deadangle_gain)';
%         figure; imagesc(output); title(num2str(nn));
    end
    A(nn,:) = output(:);
end
A(A<0) = 0;
end