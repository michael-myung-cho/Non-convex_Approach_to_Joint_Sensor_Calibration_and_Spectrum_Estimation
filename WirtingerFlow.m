function [g_rec,f_rec,Outs] = WirtingerFlow(Y_noisy_cov_rec,Opts)
% function to recover the calibration parameters via wirtinger gradient descent
% Input:      - Y_noisy_cov_rec:   covariance matrix of y
%             - opts:              MaxIter:   maximal number of iterations
%                                  OptsFun:   objective value at true
% Output:     - g_rec:             recovered calibration parameters
%             - f_rec:             recovered f

GFk          = Opts.GFk0;
Outs         = struct();
Outs.Fun     = zeros(Opts.MaxIter,1);
Outs.Fun_loss= zeros(Opts.MaxIter,1);
Outs.Fun_reg = zeros(Opts.MaxIter,1);
rho          = 1/2;        lsc          = 1/2;  %parameters in backtracking line search
for iter = 1 : Opts.MaxIter
    Funk     = WFun(GFk,Y_noisy_cov_rec,Opts.rh,Opts.n0);
    DFunk    = WDFun(GFk,Y_noisy_cov_rec,Opts.rh,Opts.n0);
    pk       = -conj(DFunk);
    DFunkpk  = - (norm(DFunk,'fro'))^2;
    alpha    = norm(Funk,'fro')/norm(DFunk,'fro');
    GFkp1    = GFk + alpha * pk;
    while WFun(GFkp1,Y_noisy_cov_rec,Opts.rh,Opts.n0)> Funk + lsc*alpha*DFunkpk
        alpha    = alpha*rho;
        GFkp1    = GFk + alpha * pk;
    end
    if alpha < Opts.alphatol 
        break;
    end
    GFk                = GFkp1;
    [Outs.Fun(iter) , Outs.Fun_loss(iter) , Outs.Fun_reg(iter)] = WFun(GFk,Y_noisy_cov_rec,Opts.rh,Opts.n0);
%     fprintf('%6.0f Fun = %e Opt = %e Loss = %e OptLoss = %e Reg = %e OptReg = %e  \n',iter,Outs.Fun(iter),Opts.OptFun,Outs.Fun_loss(iter),Opts.OptFun_loss,Outs.Fun_reg(iter),Opts.OptFun_reg);
%     fprintf('alpha = %6.5f\n',alpha)
end
% fprintf('%6.0f Fun = %e Opt = %e \n        Loss = %e OptLoss = %e \n        Reg = %e OptReg = %e  \n',iter,Outs.Fun(iter),Opts.OptFun,Outs.Fun_loss(iter),Opts.OptFun_loss,Outs.Fun_reg(iter),Opts.OptFun_reg);
% Outs.Iter = iter-1;

g_rec        = GFk(:,1);
f_rec        = GFk(:,2);
end