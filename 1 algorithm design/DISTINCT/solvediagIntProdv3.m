function [xhat] = solvediagIntProdv3(U,s,V,Dsmall,Esmall,Hsmall,lambda)
    
%D'DXE'E = D'USV'E
%xhat = D*USV'E* = (D*U)(V'E*)*s because S=diag(s)
%to add subtraction of lambda: D'DXE'E = D'USV'E - lambda*eye

    xhat = ((Dsmall\U).*(V'/Esmall')')*s-lambda*sum(inv(Hsmall),2);
    
end