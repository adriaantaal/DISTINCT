function [xhat] = solvediagIntProd(yvec,Dsmall,Esmall)
%solvediagIntProd Solves system ||A*X*D' - y||^2_2 with:
    %Diagonal structure of X
    %Nonnegative entries in X 
    
    %The algorithm exploits diagonal structure
    %minimized 2 norm directly
    
    %https://math.stackexchange.com/questions/2421545/matrix-linear-least-squares-problem-with-diagonal-matrix-constraint
    
%parameters
xhat = zeros(size(Dsmall,2),1);
for jj = 1:size(Dsmall,2)
    wjj = tovec(Dsmall(:,jj)*Esmall(:,jj)');
    xhat(jj) = yvec'*wjj/(wjj'*wjj);
end
end