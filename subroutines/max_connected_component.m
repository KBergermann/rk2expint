function [ind] = max_connected_component(G)
% Detects largest connected component of a graph.
% 
% Input:
%       G:      Matlab graph object
% 
% Output:
%       ind:    node indices of largest connected component
% 

    bins = conncomp(G);
    labels = unique(bins);
    maxval = 0; maxindex = 1;
    for i = 1 : length(labels)
        val = sum(bins==labels(i));
        if val > maxval, maxval = val; maxindex = i; end
    end
    ind = (bins==maxindex);
end