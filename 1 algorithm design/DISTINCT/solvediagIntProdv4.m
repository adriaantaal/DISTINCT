function [xhat] = solvediagIntProdv4(U,s,V,Dsmall,Esmall,sub)
    
%D'DXE'E = D'USV'E
%xhat = D*USV'E* = (D*U)(V'E*)*s because S=diag(s)
%to add subtraction of lambda: D'DXE'E = D'USV'E - lambda*eye

    xhat = ((Dsmall\U).*(V'/Esmall')')*(s-sub);
    
end