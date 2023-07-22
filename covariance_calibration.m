clear;


clc;

r = 64;

n = 64;


k = 0:n-1;    

            
t1spike = rand(1,r);%[0.2,0.2+0.5/n,0.315,0.6,0.8];%

amplitude = rand(r,1);

x = exp(-1i*2*pi*(k'*t1spike))*amplitude;

 
Y = toeplitz(x); 

%g = rand(n,1).*exp(-1i*2*pi*rand(n,1)); 
%g = g/g(1);

dynamic_range=20;

g = (1 + 10.^(rand(n,1).*(dynamic_range/20))).*exp(-1i*2*pi*rand(n,1)); 

g = g/g(1);

Ry = diag(g)*Y*diag(g');

h_truth =1./g;



abs(g'*ones(n,1))/sqrt(n)/norm(g)

%%

cvx_solver mosek%sdpt3%
cvx_begin sdp quiet
    variable H(n,n) complex hermitian
    variable T(n,n) hermitian toeplitz
    minimize trace(H)
    subject to
       Ry.*H == T; 
       H == hermitian_semidefinite(n);
       T == hermitian_semidefinite(n);
       H(2,1) == h_truth(2);

cvx_end

[u,h] = eigs(H,1);
%%
after_calibration = sqrt(h)*(u.*g);

abs(after_calibration'*ones(n,1))/sqrt(n)/norm(after_calibration)