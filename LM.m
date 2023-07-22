function [tmax , tmaxind] = LM(ell,x,RL,sparsity,band)
% function to locate the local maximum of x
t       = 0 : ell : 1-ell;
if nargin == 4
    Lband   = round(0.5*RL/ell);
else
    Lband   = round(band/ell);
end

N = length(t);
tmax = [];
tmaxind = [];

for k = 1 : N
    Neighbor = k-Lband : k+Lband;
    if Neighbor(1)<0.5
        Ind = Neighbor<0.5;
        Neighbor = Neighbor.*(1-Ind)+(Neighbor+N).*Ind;
    elseif Neighbor(end)>(N+0.5)
        Ind = Neighbor>(N+0.5);
        Neighbor = Neighbor.*(1-Ind)+(Neighbor-N).*Ind;
    end
    if x(k) == max(x(Neighbor))
        tmax = [tmax t(k)];
        tmaxind = [tmaxind k];
    end
end

[temv , temind] = sort(x(tmaxind),'descend');
temind  = temind(1:sparsity);
tmax    = tmax(temind);
tmaxind = tmaxind(temind);


