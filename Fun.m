function objective = Fun(GF,Ry)
% function to compute the objective function ||GFG'-Ry||_F^2 where 
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

objective   = sum(sum(abs(G*F*G'-Ry).^2));