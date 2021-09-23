function [glocmat,rlocmat,xlocmat,x_int_g,x_int_r,x_int] = create_randomRGlocs2D(X,Y,Z,L,deadangle,xoff,yoff)
    % Create x-vector, m random fluorophore locations with given intensity
    %Everything will be stored in micrometer domain
    %only move to rounded indices for display, otherwise reconstruction fails
    xlocmat =[];         %in indices 0 to sizes
    glocmat = [];
    rlocmat = [];
    xloc = [];

    Xdist = X(end)-X(1);
    Ydist = Y(end)-Y(1);
    Zdist = Z(end)-Z(1);

    Xsize = length(X);
    Ysize = length(Y);
    Zsize = length(Z);

    Xmid = mean(X);
    Ymid = mean(Y);

    sizes = [Xsize Ysize Zsize];
    x_int_g = zeros(prod(sizes),1);
    x_int_r = zeros(prod(sizes),1);
    x_int   = zeros(prod(sizes),1);
    % fprintf('We place %d fluorophores at uniformly random locations with mean middle of volume\n' ,mx);
    ii = 1;
    while size(xlocmat,1) < L
        xyz = rand(size(sizes)).*[Xdist*0.8 Ydist*0.8 Zdist*1]+[Xdist*0.1+X(1) Ydist*0.1+Y(1) Zdist*0.0+Z(1)];
        if ((xyz(1)>(Xmid-xoff-xyz(3)*tan(deadangle))) && (xyz(1)<(Xmid+xoff+xyz(3)*tan(deadangle))) && (xyz(2)>(Ymid-yoff-xyz(3)*tan(deadangle))) && (xyz(2)<(Ymid+yoff+xyz(3)*tan(deadangle))) && (xyz(3)>=Z(1)) && (xyz(3)<=Z(end)))
            xyzround = [X(findnearest(xyz(1),X)) Y(findnearest(xyz(2),Y)) Z(findnearest(xyz(3),Z))]; 
            xlocmat = [xlocmat;xyzround];
            xloc(ii) = sub2ind(sizes,findnearest(xyz(1),X),findnearest(xyz(2),Y),findnearest(xyz(3),Z));
            if (ii <= (L/2))
                glocmat = [glocmat;xyzround];
                x_int_g(xloc(ii)) = 1; 
            else 
                rlocmat = [rlocmat;xyzround];
                x_int_r(xloc(ii)) = 1; 
            end 
            ii = ii+1;
        end
    end    

    x_int(xloc) = 1;
end

