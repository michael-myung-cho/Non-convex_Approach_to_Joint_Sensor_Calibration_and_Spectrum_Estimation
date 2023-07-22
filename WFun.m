function [obj , obj_loss , obj_reg] = WFun(GF,Ry,rh,n0)
% function to compute the objective function ||GFG'-Ry||_F^2 + rh*G0 where 
%                     G = diag(g) with g = gr + 1i * gi
%                     F = diag(f,1:-1:-N)+diag(conj(f(2:end)),2:N)

N     = size(GF,1);
g     = GF(:,1);    
f     = GF(:,2);

G     = diag(g);
F     = diag(f(1)*ones(N,1));
for k = 2:N
    F = F + diag(f(k)*ones(N-k+1,1),-(k-1)) + diag(conj(f(k))*ones(N-k+1,1),k-1);
end

obj_loss   = sum(sum(abs(G*F*G'-Ry).^2));
obj_reg    = rh*max( norm(f)^2/(2*n0)-1 , 0 )^2 + rh*max( norm(g)^2/sqrt(2*n0)-1 , 0 )^2;
obj        = obj_loss + obj_reg;