function s = symmetricity(H)
%metrix for symmetry. -1 fully antisymmetric, +1 fully symmetric
%https://math.stackexchange.com/questions/2048817/metric-for-how-symmetric-a-matrix-is
    symmH = 0.5*(H+H'); 
    asymmH = 0.5*(H-H'); 
    s = (norm(symmH,'fro')-norm(asymmH,'fro'))/(norm(symmH,'fro')+norm(asymmH,'fro'));    
end

