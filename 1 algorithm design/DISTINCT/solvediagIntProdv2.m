function [xhat] = solvediagIntProdv2(U,s,V,Dsmall,Esmall)
%solvediagIntProd Solves system ||A*X*D' - y||^2_2 with:
    %Diagonal structure of X
    %Nonnegative entries in X 
    
    %The algorithm exploits diagonal structure
    %minimized 2 norm directly
    
xhat = ((Dsmall\U).*(V'/Esmall')')*s;

end