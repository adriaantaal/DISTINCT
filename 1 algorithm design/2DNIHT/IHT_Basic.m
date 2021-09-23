function [ y ] = IHT_Basic( x,Phi,M,mu,epsilon,loopmax )    
%IHT_Basic Summary of this function goes here    
%Version: 1.0 written by jbb0523 @2016-07-30    
%Reference:Blumensath T, Davies M E. Iterative Thresholding for Sparse Approximations[J].   
%Journal of Fourier Analysis & Applications, 2008, 14(5):629-654.   
%(Available at: http://link.springer.com/article/10.1007%2Fs00041-008-9035-z)  
%   Detailed explanation goes here    
    if nargin < 6    
        loopmax = 3000;    
    end    
    if nargin < 5      
        epsilon = 1e-3;      
    end     
    if nargin < 4      
        mu = 1;      
    end     
    [x_rows,x_columns] = size(x);      
    if x_rows<x_columns      
        x = x';%x should be a column vector      
    end    
    n = size(Phi,2);    
    y = zeros(n,1);%Initialize y=0    
    loop = 0;    
    while(norm(x-Phi*y)>epsilon && loop < loopmax)    
        y = y + Phi'*(x-Phi*y)*mu;%update y    
        %the following two lines of code realize functionality of H_M(.)    
        %1st: permute absolute value of y in descending order    
        [ysorted inds] = sort(abs(y), 'descend');    
        %2nd: set all but M largest coordinates to zeros    
        y(inds(M+1:n)) = 0;    
        loop = loop + 1;    
    end    
end   