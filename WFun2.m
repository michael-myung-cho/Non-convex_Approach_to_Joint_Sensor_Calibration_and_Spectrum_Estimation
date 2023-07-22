function [obj , obj_loss , obj_reg1 , obj_reg2] = WFun2(HU,Ry,rh1,rh2)
% function to compute the objective function ||diag(h) Ry diag(hbar) - Toep(u)||_F^2 

N     = size(HU,1);
h     = HU(:,1);    
u     = HU(:,2);

H     = diag(h);
Toepu =diag(u(1)*ones(N,1));
for k = 2:N
    Toepu = Toepu + diag(u(k)*ones(N-k+1,1),-(k-1)) + diag(conj(u(k))*ones(N-k+1,1),k-1);
end


obj_loss   = sum(sum(abs(Toepu-H*Ry*H').^2));
obj_reg1   = rh1*max( 1 - norm(h)^2/(1/4) , 0 )^2; 
obj_reg2   = rh2*max( norm(h)^2/(4)-1 , 0 )^2; 
obj        = obj_loss + obj_reg1 + obj_reg2;

% regularize on u
% obj_reg1   = rh1*max( 1 - norm(u)^2/(n0/4) , 0 )^2; 
% obj_reg2   = rh2*max( norm(u)^2/(4*n0)-1 , 0 )^2; 