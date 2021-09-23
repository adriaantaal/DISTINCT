function acc2 = some_accuracy_function(xhatDS,x_int)

%     acc = xhatDS'*x_int/(sum(xhatDS+x_int));

%     correlationmatrix = corrcoef(xhatDS,x_int);
%     acc1 = correlationmatrix(2,1);
      
      %in-crowd uses average of correlation coefficients as accuracy

    locationmatrix = zeros(size(xhatDS));
    for i = 1:length(xhatDS)
        if xhatDS(i) ~= 0
                locationmatrix(i) = 1;
        end
    end
    
    correlationmatrix_method2 = corrcoef(locationmatrix,x_int);
    acc2 = correlationmatrix_method2(2,1);
    
    
%     sizematrix = size(find(x_int));
%     acc3 = xhatDS'*x_int/(sizematrix(1,1)+sum(x_int));
        
        
    
        

 
 % Magnitude of inner product/magnitude of
 
 %   acc = confusionmat(xhatDS,x_int);
 
 % confusion matrix - classify calculated locations as either in or out of
 % a certain range of distance for any of the sources, and use a confusion
 % matrix. Would still be problematic for large number of sources.
 
 %Maybe break up function by number of sources - i.e. low number of
 %sources, evaluate accuracy by doing euclidian distance, while large
 %number of sources, maybe finding mean and SD x,y,z positions of sources and
 %comparing with mean and SD x,y,z of xhadDS
 
 
 
 
        
    
    

end