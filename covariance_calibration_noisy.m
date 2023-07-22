clear;


clc;

r = 5;

n = 16;


k = 0:n-1;    

            
t1spike = rand(1,r);%[0.2,0.2+0.5/n,0.315,0.6,0.8];%

dynamic_range=10;

for monte = 1:10
g = (1 + 10.^(rand(n,1).*(dynamic_range/20))).*exp(-1i*2*pi*rand(n,1)); 

g = g/g(1);

h_truth =1./g;


abs(g'*ones(n,1))/sqrt(n)/norm(g)

Llist = 100:100:2000;
%%
for iter = 1:length(Llist)
    %iter,
    L = Llist(iter);

amplitude = randn(r,L);

x = exp(1i*2*pi*(k'*t1spike))*amplitude + (randn(n,L)+1i*randn(n,L))/sqrt(2)*0.05;

 
Y = x*x'/L; 

Ry = diag(g)*Y*diag(g');
%

cvx_solver mosek%sdpt3%
cvx_begin sdp quiet
    variable H(n,n) complex hermitian
    variable T(n,n) hermitian toeplitz
    minimize square_pos(norm(Ry.*H - T,'fro'))+ 0.5* trace(H)
    subject to
       %Ry.*H == T;
       H == hermitian_semidefinite(n);
       T == hermitian_semidefinite(n);
       H(2,1) == h_truth(2);

cvx_end

[u,h] = eigs(H,1);
%
after_calibration = sqrt(h)*(u.*g);

after_calibration = after_calibration*exp(-1i*angle(after_calibration(1)));

localization_after=rootmusic(T,r,'corr')/(2*pi);
localization_before=rootmusic(Ry,r,'corr')/(2*pi);
localization_sample=rootmusic(Y,r,'corr')/(2*pi);

err_after(iter,monte) = H_dist(localization_after,t1spike);
err_before(iter,monte) = H_dist(localization_before,t1spike);
err_sample(iter,monte) = H_dist(localization_sample,t1spike);
%err(iter) = abs(after_calibration'*ones(n,1))/sqrt(n)/norm(after_calibration)
end



end
%%
figure;hold on;
plot(Llist,mean(err_after,2),'bo-'); hold on;
plot(Llist,mean(err_before,2),'rx-'); 
plot(Llist,mean(err_sample,2),'ms-'); 
legend('after','before','sample');