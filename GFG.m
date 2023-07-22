function gfg = GFG(GF)
% function to compute the matrix GFG' where
%                     G = diag(g) with g = gr + 1i * gi
%                     F = diag(f,1:-1:-N)+diag(conj(f(2:end)),2:N)

N     = size(GF,1);
gr    = GF(:,1);    gi = GF(:,2);     fr = GF(:,3);     fi = GF(:,4);

G     = diag(gr+1i*gi);
f     = fr + 1i * fi;

F     = diag(f(1)*ones(N,1));
for k = 2:N
    F = F + diag(f(k)*ones(N-k+1,1),-(k-1)) + diag(conj(f(k))*ones(N-k+1,1),k-1);
end

gfg   = G*F*G';