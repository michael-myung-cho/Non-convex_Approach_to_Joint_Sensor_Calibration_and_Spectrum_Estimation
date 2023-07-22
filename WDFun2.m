function D = WDFun2(HU,Ry,rh1,rh2)
% function to compute the gradient of the objective function through
% Wirtinger derivative
% ||diag(h) Ry diag(hbar) - Toep(u)||_F^2  
N   =size(HU,1);                                  % number of sensors 

%% Form matrix G and F                   
h   = HU(:,1);      hbar = conj(h);
u   = HU(:,2);      ubar = conj(u);

%% Compute the Wirtinger Derivative 
Dh  =  zeros(N,1); 
Du  =  zeros(N,1); 

for k = 1 : N
    Dh(k) =((hbar(1:k).').*Ry(k,1:k))*(hbar(k)*(h(1:k).*Ry(k,1:k)')-ubar(k:-1:1))...
           +(hbar(k:N).*conj(Ry(k:N,k))).'*(hbar(k)*(h(k:N).*Ry(k:N,k))-u(1:N-k+1));
end

for k = 1:N
    Du(k)=-sum(conj(h(k:N).*hbar(1:N-k+1).*diag(Ry,-k+1)-u(k)));
end

Dh    = Dh - rh1/(1/4) * 2*max( 1 - norm(h)^2/(1/4) , 0 ) * hbar + rh2/(4*1) * 2*max( norm(h)^2/(4) -1 , 0 ) * hbar;

D = [Dh Du];







