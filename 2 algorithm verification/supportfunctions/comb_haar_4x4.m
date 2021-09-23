function [comb_haar] = comb_haar_4x4()

comb_haar(1,:) = logical([1 1 0 0 1 1 0 0 0 0 1 1 0 0 1 1]); 
comb_haar(2,:) = not(comb_haar(1,:));
% %second order
comb_haar(3,:) = [1 0 1 0 0 0 0 0 1 0 1 0 0 0 0 0];
comb_haar(4,:) = [0 1 0 1 0 0 0 0 0 1 0 1 0 0 0 0];
comb_haar(5,:) = [0 0 0 0 1 0 1 0 0 0 0 0 1 0 1 0];
comb_haar(6,:) = [0 0 0 0 0 1 0 1 0 0 0 0 0 1 0 1];

end

