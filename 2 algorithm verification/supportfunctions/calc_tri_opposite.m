function dxe = calc_tri_opposite(D,X,E,activeind)
%when D is a wide matrix
    DX = D*X(:,activeind);
    dxe = DX*E(activeind,:); 
end

