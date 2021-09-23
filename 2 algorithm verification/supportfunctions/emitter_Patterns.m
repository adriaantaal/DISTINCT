function [combinations] = emitter_Patterns(minEmitter,clustersize)
%EMITTER_PATTERNS
%   minEmitter: minimum of emitters on within cluster 
        %slashes computation cost and guarantees minimum signal
%   clusterSize: number of pixels in repeated cluster

    combinations = dec2bin(2^clustersize-1:-1:1)-'0';
    combinations(sum(combinations,2)<minEmitter,:) = [];
end

