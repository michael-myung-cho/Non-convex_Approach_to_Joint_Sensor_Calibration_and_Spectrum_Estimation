function [supp , suppidx] = GenerateFrequency(sparsity,sep1,sep2,stype,RL)
%% function to generate sources with the separation condition given
%
% Input:  sparsity: number of sources
%         sep1: minimal separation
%         sep2: maximal separation
%
% Output: supp: frequency locations \in [0,1)

if strcmp(stype,'on grid') && nargin < 5
    display('sources on grid, then input must be 5'); return;
end
if sparsity*sep2>1
    display('sparsity and separation not compatible: sparsity too large'); return;
end


supp = 0.1;
while (supp(end)<1-2*sep2) && (length(supp)<sparsity)
    move = sep1+rand*(sep2-sep1);
    supp = [supp ; supp(end)+move];
end

switch stype
    case 'on grid'
        suppidx = round(supp/RL);
        supp    =  suppidx*RL;
    case 'off grid'
        suppidx = [];
end
