function [g_rec,f_rec,Outs] = OptimizationMethod(Y_noisy_cov_rec,Opts)
% function to recover the calibration parameters via gradient descent
% Input:      - Y_noisy_cov_rec:   covariance matrix of y
%             - opts:              MaxIter:   maximal number of iterations
%                                  OptsFun:   objective value at true
% Output:     - g_rec:             recovered calibration parameters
%             - f_rec:             recovered f

GFk          = Opts.GFk0;
Outs         = struct();
Outs.Fun     = zeros(Opts.MaxIter,1);
rho          = 1/2;        lsc          = 1/2;  %parameters in backtracking line search
for iter = 1 : Opts.MaxIter
    Funk     = Fun(GFk,Y_noisy_cov_rec);
    DFunk    = DFun(GFk,Y_noisy_cov_rec);
    pk       = -DFunk;
    DFunkpk  = - (norm(DFunk,'fro'))^2;
    alpha    = norm(Funk,'fro')/norm(DFunk,'fro');
    GFkp1    = GFk + alpha * pk;
    while Fun(GFkp1,Y_noisy_cov_rec)> Funk + lsc*alpha*DFunkpk
        alpha    = alpha*rho;
        GFkp1    = GFk + alpha * pk;
    end
    if alpha < Opts.alphatol 
        break;
    end
    GFk                = GFkp1;
    Outs.Fun(iter)     = Fun(GFk,Y_noisy_cov_rec);
%     fprintf('iter = %6.0f  Fun = %e Opt = %e  ',iter,Outs.Fun(iter),Opts.OptFun);
%     fprintf('alpha = %6.5f\n',alpha)
end
fprintf('iter = %6.0f  Fun = %e Opt = %e  \n',iter,Outs.Fun(iter),Opts.OptFun);
Outs.Iter = iter-1;

g_rec        = GFk(:,1) + 1i * GFk(:,2);
f_rec        = GFk(:,3) + 1i * GFk(:,4);
