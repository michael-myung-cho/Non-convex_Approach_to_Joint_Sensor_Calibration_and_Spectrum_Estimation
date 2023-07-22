function [h_rec,u_rec,Outs] = WirtingerFlow2(Y_noisy_cov_rec,Opts)
% function to recover the calibration parameters via wirtinger gradient descent
% Input:      - Y_noisy_cov_rec:   covariance matrix of y
%             - opts:              MaxIter:   maximal number of iterations
%                                  OptsFun:   objective value at true
% Output:     - u_rec:             recovered calibration parameters inverse
%             - u_rec:             recovered u

HUk          = Opts.HUk0;
Outs         = struct();
Outs.Fun     = zeros(Opts.MaxIter,1);
Outs.Fun_loss= zeros(Opts.MaxIter,1);
Outs.Fun_reg1= zeros(Opts.MaxIter,1);
Outs.Fun_reg2= zeros(Opts.MaxIter,1);
rho          = 1/2;        lsc          = 1/2;  %parameters in backtracking line search
for iter = 1 : Opts.MaxIter
    Funk     = WFun2(HUk,Y_noisy_cov_rec,Opts.rh1,Opts.rh2);
    DFunk    = WDFun2(HUk,Y_noisy_cov_rec,Opts.rh1,Opts.rh2);
    pk       = -conj(DFunk);
    DFunkpk  = - (norm(DFunk,'fro'))^2;
    alpha    = norm(Funk,'fro')/norm(DFunk,'fro');
    HUkp1    = HUk + alpha * pk;
    while WFun2(HUkp1,Y_noisy_cov_rec,Opts.rh1,Opts.rh2)> Funk + lsc*alpha*DFunkpk
        alpha    = alpha*rho;
        HUkp1    = HUk + alpha * pk;
    end
%     if alpha < Opts.alphatol 
%         break;
%     end
    HUk                = HUkp1;
    [Outs.Fun(iter) , Outs.Fun_loss(iter) , Outs.Fun_reg1(iter) , Outs.Fun_reg2(iter)] = WFun2(HUk,Y_noisy_cov_rec,Opts.rh1,Opts.rh2);
%    fprintf('%6.0f Fun = %6.4f Opt = %6.4f Loss = %6.4f OptLoss = %6.4f Reg1 = %6.4f Reg2 = %6.4f \n',iter,Outs.Fun(iter),Opts.OptFun,Outs.Fun_loss(iter),Opts.OptFun_loss,Outs.Fun_reg1(iter),Outs.Fun_reg2(iter));
%    fprintf('alpha = %6.5f\n',alpha)
end
Outs.Iter = iter-1;
%fprintf('%6.0f Fun = %6.4f Opt = %6.4f Loss = %6.4f OptLoss = %6.4f Reg1 = %6.4f OptReg1 =%6.4f Reg2 = %6.4f OptReg2 =%6.4f \n',Outs.Iter,Outs.Fun(Outs.Iter),Opts.OptFun,Outs.Fun_loss(Outs.Iter),Opts.OptFun_loss,Outs.Fun_reg1(Outs.Iter),Opts.OptFun_reg1,Outs.Fun_reg2(Outs.Iter),Opts.OptFun_reg2);


h_rec        = HUk(:,1);
u_rec        = HUk(:,2);
