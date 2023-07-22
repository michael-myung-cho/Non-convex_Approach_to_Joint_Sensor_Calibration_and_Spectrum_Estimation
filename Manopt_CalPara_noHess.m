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