function D = WDFun(GF,Ry,rh,n0)
% function to compute the gradient of the objective function through
% Wirtinger derivative
% ||GfG'-Ry||_F^2 where G=diag(g) with g=gr+i*gi
%                       G = diag(g) with g = gr + 1i * gi
%                       F = diag(f,1:-1:-N)+diag(conj(f(2:end)),2:N)  

N   =size(GF,1);                                  % number of sensors 

%% Form matrix G and F                   
g   = GF(:,1);      gbar = conj(g);
f   = GF(:,2);      fbar = conj(f);
F=diag(f(1)*ones(N,1));
for k = 2:N
    F = F + diag(f(k)*ones(N-k+1,1),-(k-1)) + diag(conj(f(k))*ones(N-k+1,1),k-1);
end

%% Compute the Wirtinger Derivative 
DG  =  zeros(N,1); 
DF  =  zeros(N,1); 

for k = 1 : N
    DG(k) =(gbar(1:k).*f(k:-1:1)).'*(gbar(k)*(g(1:k).*fbar(k:-1:1))-Ry(k,1:k)')+(gbar(k:N).*fbar(1:N-k+1)).'*(gbar(k)*(g(k:N).*f(1:N-k+1))-Ry(k:N,k));
end

DG    = DG + rh/sqrt(2*n0) * 2*max( norm(g)^2/sqrt(2*n0) -1 , 0) * gbar;

DF(1) =((g.*gbar)*f(1)-diag(Ry))'*(g.*gbar);
for k=2:N
    DF(k)=(g(k:N).*gbar(1:N-k+1)*f(k)-diag(Ry,-k+1))'*(g(k:N).*gbar(1:N-k+1));
end

DF    = DF + rh/(2*n0) * 2*max( norm(f)^2/(2*n0) -1 , 0 ) * fbar;

% DF(1) =((g.*gbar)*fbar(1)-diag(Ry))'*(g.*gbar);
% for k=2:N
%     DF(k)=(gbar(k:N).*g(1:N-k+1)*fbar(k)-diag(Ry,-k+1)).'*(g(k:N).*gbar(1:N-k+1));
% end


D = [DG DF];







% %% Compute the Wirtinger Derivative 
% DGr=zeros(N,1); DGi=zeros(N,1); DFr=zeros(N,1); DFi=zeros(N,1);
% 
% for k=1:N
%     DGr(k)=(gbar(1:k).*f(k:-1:1)).'*(gbar(k)*(g(1:k).*fbar(k:-1:1))-Ry(k,1:k)')+(gbar(k:N).*fbar(1:N-k+1)).'*(gbar(k)*(g(k:N).*f(1:N-k+1))-Ry(k:N,k));
%     DGi(k)=-2*imag(DGr(k));
%     DGr(k)=2*real(DGr(k));
%    %DGi(k)=-2*imag(2*DGr(k)-2*f(1)*gbar(k)*(f(1)*g(k)*gbar(k)-Ry(k,k)));
%    %DGr(k)=2*real(2*DGr(k)-2*f(1)*gbar(k)*(f(1)*g(k)*gbar(k)-Ry(k,k)));
% end
% 
% %DGr(N)=2*(gbar(1:N).*f(N:-1:1)).'*(gbar(N)*(g(1:N).*fbar(N:-1:1))-Ry(N,1:N)')-f(1)*gbar(N)*(f(1)*g(N)*gbar(N)-Ry(N,N));
% %DGi(N)=-1*imag(DGr(N));
% %DGr(N)=1*real(DGr(N));
% 
% 
% DFr(1)=((g.*gbar)*fbar(1)-diag(Ry))'*(g.*gbar);
% DFi(1)=2*imag(DFr(1));
% DFr(1)=2*real(DFr(1));
% 
% for k=2:N
%     DFr(k)=(gbar(k:N).*g(1:N-k+1)*fbar(k)-diag(Ry,-k+1)).'*(g(k:N).*gbar(1:N-k+1));
%     DFi(k)=-2*imag(DFr(k));
%     DFr(k)=2*real(DFr(k));
% end
% 
%     
%     
% D = [DGr DGi DFr DFi];
% 
% 
% 
% 
