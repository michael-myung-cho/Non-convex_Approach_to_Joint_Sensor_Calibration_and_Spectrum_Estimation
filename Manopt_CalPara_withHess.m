function [ x, Tu ] = Manopt_CalPara_noHess( Y_cov,h0 )
%
% n: signal dimension
% r: rank
% Y_cov: covariance matrix
%
% -----------------------

%% initialization
n=size(Y_cov,1);
 
%% Create the problem structure.
manifold = spherecomplexfactory(n, 1);
problem.M = manifold;
 
%% Define the problem cost function and its Euclidean gradient.
problem.cost  = @(h) norm(diag(h)*Y_cov*diag(h') - toep(diag(h)*Y_cov*diag(h')),'fro')^2;
problem.egrad = @(h) 2*diag(conj(Y_cov'*diag(h')*(diag(h)*Y_cov*diag(h') - toep(diag(h)*Y_cov*diag(h'))))); 
problem.ehess = @hess;
function H11 = hess(h)
    for a=0:n-1
        for b=0:n-1
            sum1=0;
            sum2=0;
            sum3=0;
            sum4=0;
            sum5=0;
            if b <= n-a-1
                for j=0:a
                    sum1=sum1+conj(Y_cov(j+1,a+1))*conj(h(j+1))/(n+j-a)*Y_cov(b+1,b+1-j+a)*conj(h(b+1-j+a));
                end
                for j=a+1:a+b
                    sum2=sum2+conj(Y_cov(j+1,a+1))*conj(h(j+1))/(n-j+a)*Y_cov(b+1,b+1-j+a)*conj(h(b+1-j+a));
                end
            else
                for j=b-n+a+1:a
                    sum1=sum1+conj(Y_cov(j+1,a+1))*conj(h(j+1))/(n+j-a)*Y_cov(b+1,b+1-j+a)*conj(h(b+1-j+a));
                end
                for j=a+1:n-1
                    sum2=sum2+conj(Y_cov(j+1,a+1))*conj(h(j+1))/(n-j+a)*Y_cov(b+1,b+1-j+a)*conj(h(b+1-j+a));
                end
            end
            if b <= a
                for j=a-b:a
                    sum3=sum3+conj(Y_cov(j+1,a+1))*conj(h(j+1))/(n+j-a)*Y_cov(b+1+j-a,b+1)*h(b+1+j-a);
                end
                for j=a-b:n-1
                    sum4=sum4+Y_cov(j+1+b-a,j+1)*h(j+1+b-a)*conj(h(j+1));
                end
                sum4=conj(Y_cov(b+1,a+1))/(n+b-a)*sum4;
                for j=a+1:n-1
                    sum5=sum5+conj(Y_cov(j+1,a+1))*conj(h(j+1))/(n-j+a)*Y_cov(b+1+j-a,b+1)*conj(h(b+1+j-a));
                end                   
            else 
                for j=0:a
                    sum3=sum3+conj(Y_cov(j+1,a+1))*conj(h(j+1))/(n+j-a)*Y_cov(b+1+j-a,b+1)*h(b+1+j-a);
                end
                for j=0:n-b+a-1
                    sum4=sum4+Y_cov(j+1+b-a,j+1)*h(j+1+b-a)*conj(h(j+1));
                end
                sum4=conj(Y_cov(b+1,a+1))/(n-b+a)*sum4;
                for j=a+1:n+a-1-b
                    sum5=sum5+conj(Y_cov(j+1,a+1))*conj(h(j+1))/(n-j+a)*Y_cov(b+1+j-a,b+1)*conj(h(b+1+j-a));
                end                     
            end
            H00(a+1,b+1)=2*conj(Y_cov(b+1,a+1))*Y_cov(b+1,a+1)*conj(h(b+1))*cong(h(a+1))-2*sum1-2*sum2;
            H11(a+1,b+1)=conj(H00(a+1,b+1));
            H01(a+1,b+1)=2*conj(Y_cov(b+1,a+1))*Y_cov(b+1,a+1)*h(b+1)*cong(h(a+1))-2*sum3-2*sum4-2*sum5;
            H10(a+1,b+1)=conj(H10(a+1,b+1));
        end
    end
    H=[H00,H01;H10,H11];
end

%% Numerically check gradient consistency (optional).
checkgradient(problem);
 
%% Solve.
% options.maxiter = 200;
% options.batchsize = 10;
% options.stepsize_type = 'decay';
% options.stepsize_init = 1e2;
% options.stepsize_lambda = 1e-3;
% options.verbosity = 2;
%[x, cost, info, options] = conjugategradient(problem,WGFk0(:,1));

%options.rho_prime = 0.249;
%options.useRand = true;
% options.rho_regularization = 1e6;
[x, xcost, info, options] = trustregions(problem,h0);
%g     = WGFk0(:,1);  
%[x, xcost, info, options] = trustregions(problem,1./g);


%% Display some statistics.
% figure;
% semilogy([info.iter], [info.gradnorm], '.-');
% xlabel('Iteration number');
% ylabel('Norm of the gradient of f');

Tu=toep(Y_cov.*(x*x'));
%save('manopt.mat');
end

function [ B ] = toep( A )
    [n,m]=size(A);
    ind=1;
    for ii=0:-1:-n+1
        vA(ind)=sum(diag(A,ii))/size(diag(A,ii),1);
        ind=ind+1;
    end
    for ii=0:n-1
        vA(ind)=sum(diag(A,ii))/size(diag(A,ii),1);
        ind=ind+1;
    end
    B=toeplitz(vA(1:n),vA(n+1:end));
end