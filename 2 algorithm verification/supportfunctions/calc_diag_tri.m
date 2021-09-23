function dxediag = calc_diag_tri(D,X,E,activeind)
%when D is a skinny matrix
    N = length(activeind);
    DX = D(activeind,:)*X;
    dxediag = zeros(N,1); 
    for nn = 1:N
        dxediag(nn,1) = DX(nn,:)*E(:,activeind(nn)); 
    end
end

