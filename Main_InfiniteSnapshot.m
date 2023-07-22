%% Main code: source localization + sensor calibration with multiple snapshot of measurements
clear all
close all

display('============================================================')
%% Set-up
sparsity = 2;                % number of sources
stype    = 'off grid';        % source type: on gird - sources on the RL grid; off grid - sources on the continuum.
N        = 10;                % number of sensors
RL       = 1/N;               % Rayleigh length (resolution limit)
srange   = 2;                 % range of source amplitudes
crange   = 2;                 % range of calibration amplitudes
nlevel   = 100;                 % percentage of noise
fprintf('%6.0f sources %6.0f measurements  ',sparsity,N)
fprintf('noise level = %6.0f \n',nlevel)

%% Generate source location
sepn1               = 2;
sepn2               = 3;
sep1                = sepn1*RL;
sep2                = sepn2*RL;
[Supp , SuppIdx]    = GenerateFrequency(sparsity,sep1,sep2,stype,RL);

%% Generate source amplitudes
Samp     = 1+rand(sparsity,1)*(srange-1);

%% Sensor gain
Alpha    = 1+rand(N,1)*(crange-1);
Beta     = 2*pi*rand(N,1);
Gain     = Alpha.*exp(1i*Beta); % calibration para. (MCho)
AlphaL2  = norm(Alpha);


%% Generate measurement
A        = exp(2*pi*1i*((0:N-1)')*(Supp'));
f        = A*(Samp).^2/N;
F        = diag(f(1)*ones(N,1));
for k = 2:N
    F = F + diag(f(k)*ones(N-k+1,1),-(k-1)) + diag(conj(f(k))*ones(N-k+1,1),k-1);
end


%% Noisey measurement
Y_cov             = diag(Gain)*F*diag(conj(Gain));
nsigma            = nlevel/100*norm(Y_cov,'fro')/(sqrt(2)*N);
E                 = normrnd(0,nsigma,N,N)+1i*normrnd(0,nsigma,N,N);
%E                 = normrnd(0,nlevel/100,N,N)+1i*normrnd(0,nlevel/100,N,N); %MCho
Y_noisy_cov       = Y_cov + E;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y_noisy_cov_rec   = Y_noisy_cov;



ell               = 1/100*RL;
Grid              = 0 : ell : 1;
GridLen           = length(Grid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Method 1: algebraic method
% display('===================Algebraic method=========================')
% 
% % Reconstruction
% [Gain_rec ,  ImagingFun , f_rec] = AlgebraicMethod(Y_noisy_cov_rec,sparsity,ell,AlphaL2);
% n0                   = norm(Gain_rec)^2*norm(f_rec);
% Gain_rec_algebraic   = Gain_rec;
% f_rec_algebraic      = f_rec;
% % Local Maximum extraction in MUSIC
% [Supp_rec , supprecind] = LM(ell,ImagingFun,RL,sparsity);
% % Compute error and adjust trivial ambiguities
% [Supp_rec_c2,Beta_rec_c12,Gain_rec_c,Supp_err,Alpha_reerr,Beta_reerr,Gain_reerr,c0,c1,c2,c] = ComputeReconstructionError(Supp,Supp_rec,Gain,Gain_rec,ell);
% fprintf('Support recovery error = %6.4f RL\n',Supp_err)
% fprintf('Amplitude error mean =%6.4f 0.9 0.5 0.1 quantitle =  %6.4f %6.4f %6.4f \n',mean(Alpha_reerr),quantile(Alpha_reerr,0.9),median(Alpha_reerr),quantile(Alpha_reerr,0.1) );
% fprintf('Phase     error mean =%6.4f 0.9 0.5 0.1 quantitle = %6.4f %6.4f %6.4f \n',mean(Beta_reerr),quantile(Beta_reerr,0.9),median(Beta_reerr),quantile(Beta_reerr,0.1) );
% fprintf('Cali      error mean =%6.4f 0.9 0.5 0.1 quantitle = %6.4f %6.4f %6.4f \n',mean(Gain_reerr),quantile(Gain_reerr,0.9),median(Gain_reerr),quantile(Gain_reerr,0.1) );
% fprintf('cali constant        = %6.2f + i %6.2f  magnitude =%6.2f\n',real(c),imag(c),abs(c))
% % Imaging function
% figure; plot(mod(Grid+c2/(2*pi),1),ImagingFun,'b*-'); hold on;
% stem(Supp,zeros(sparsity,1),'ro','filled');
% title('Imaging function in the algebraic method')
% xlabel(['support error = ' num2str(Supp_err) ' RL'])

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Method 2: full algebraic method
% display('===================Full Algebraic method=========================')
%
% % Reconstruction
% [Gain_rec ,  ImagingFun] = AlgebraicMethodFull(Y_noisy_cov_rec,sparsity,ell,AlphaL2);
% % Local Maximum extraction in MUSIC
% [Supp_rec , supprecind] = LM(ell,ImagingFun,RL,sparsity);
% % Compute error and adjust trivial ambiguities
% [Supp_rec_c2,Beta_rec_c12,Gain_rec_c,Supp_err,Alpha_reerr,Beta_reerr,Gain_reerr,c0,c1,c2,c] = ComputeReconstructionError(Supp,Supp_rec,Gain,Gain_rec,ell);
% fprintf('Support recovery error = %6.4f RL\n',Supp_err)
% fprintf('Amplitude error mean =%6.4f 0.9 0.5 0.1 quantitle =  %6.4f %6.4f %6.4f \n',mean(Alpha_reerr),quantile(Alpha_reerr,0.9),median(Alpha_reerr),quantile(Alpha_reerr,0.1) );
% fprintf('Phase     error mean =%6.4f 0.9 0.5 0.1 quantitle = %6.4f %6.4f %6.4f \n',mean(Beta_reerr),quantile(Beta_reerr,0.9),median(Beta_reerr),quantile(Beta_reerr,0.1) );
% fprintf('Cali      error mean =%6.4f 0.9 0.5 0.1 quantitle = %6.4f %6.4f %6.4f \n',mean(Gain_reerr),quantile(Gain_reerr,0.9),median(Gain_reerr),quantile(Gain_reerr,0.1) );
% fprintf('cali constant        = %6.2f + i %6.2f  magnitude =%6.2f\n',real(c),imag(c),abs(c))
% % Imaging function
% figure; stem(Supp,Samp,'ro'); hold on;
% plot(Grid,ImagingFun,'b*-')
% title('Imaging function for the algebraic method')
% xlabel(['support error = ' num2str(Supp_err) ' RL'])
%
%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Method 3: Optimization
% display('===================Optimization with real variables======================')
%
% Initialization
% GFk0          = rand(N,4);
% GFk0          = zeros(N,4);
% GFk0(:,1)     = real(n0^(1/4)*Gain_rec_algebraic/norm(Gain_rec_algebraic)); %GFk0(:,1:2)*crange/sqrt(2);
% GFk0(:,2)     = imag(n0^(1/4)*Gain_rec_algebraic/norm(Gain_rec_algebraic));
% GFk0(:,3)     = real(sqrt(n0)*f_rec_algebraic/norm(f_rec_algebraic));      %GFk0(:,1:2)*srange^2/sqrt(2)*sparsity/N;
% GFk0(:,4)     = imag(sqrt(n0)*f_rec_algebraic/norm(f_rec_algebraic));
% Opts          = struct('MaxIter',100,'GFk0',GFk0,'alphatol',10^(-4));
% Opts.OptFun   = Fun([real(Gain) imag(Gain) real(f) imag(f)],Y_noisy_cov_rec);
%
% % Gradient descent
% [g_rec,f_rec,Outs] = OptimizationMethod(Y_noisy_cov_rec,Opts);
%
% figure; plot(log(2:length(Outs.Fun)),log(Outs.Fun(2:end)),'b*-');hold on;
% plot(log(2:length(Outs.Fun)),ones(length(Outs.Fun)-1,1)*log(Opts.OptFun),'r')
% title('Objective function versus iteration in log-log scale')
%
% Alpha_rec    = abs(g_rec);
% Alpha_rec    = Alpha_rec/norm(Alpha_rec)*AlphaL2;
% Beta_rec     = angle(g_rec);
% Gain_rec     = Alpha_rec.*exp(1i*Beta_rec);
%
% % MUSIC
% ImagingFun   = MUSIC(f_rec,sparsity,ell);
% % Local Maximum extraction in MUSIC
% [Supp_rec , supprecind] = LM(ell,ImagingFun,RL,sparsity);
% % Compute error and adjust trivial ambiguities
% [Supp_rec_c2,Beta_rec_c12,Gain_rec_c,Supp_err,Alpha_reerr,Beta_reerr,Gain_reerr,c0,c1,c2,c] = ComputeReconstructionError(Supp,Supp_rec,Gain,Gain_rec,ell);
% fprintf('Support recovery error = %6.4f RL\n',Supp_err)
% fprintf('Amplitude error mean =%6.4f 0.9 0.5 0.1 quantitle =  %6.4f %6.4f %6.4f \n',mean(Alpha_reerr),quantile(Alpha_reerr,0.9),median(Alpha_reerr),quantile(Alpha_reerr,0.1) );
% fprintf('Phase     error mean =%6.4f 0.9 0.5 0.1 quantitle = %6.4f %6.4f %6.4f \n',mean(Beta_reerr),quantile(Beta_reerr,0.9),median(Beta_reerr),quantile(Beta_reerr,0.1) );
% fprintf('Cali      error mean =%6.4f 0.9 0.5 0.1 quantitle = %6.4f %6.4f %6.4f \n',mean(Gain_reerr),quantile(Gain_reerr,0.9),median(Gain_reerr),quantile(Gain_reerr,0.1) );
% fprintf('cali constant        = %6.2f + i %6.2f  magnitude =%6.2f\n',real(c),imag(c),abs(c))
%
%
% % Imaging function
% figure; plot(mod(Grid+c2/(2*pi),1),ImagingFun,'b*-'); hold on;
% stem(Supp,zeros(sparsity,1),'ro','filled');
% title('Imaging function in the optimization approach')
% xlabel(['support error = ' num2str(Supp_err) ' RL'])
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Method 2: Optimization with Wirtinger flow
% display('======================Witinger flow=========================')
%
% % Initialization
% rh            = 0;
% WGFk0         = [GFk0(:,1)+1i*GFk0(:,2) GFk0(:,3)+1i*GFk0(:,4)];
% Opts          = struct('MaxIter',500,'GFk0',WGFk0,'alphatol',10^(-4),'rh',rh,'n0',n0);
% [Opts.OptFun , Opts.OptFun_loss , Opts.OptFun_reg]   = WFun([Gain f],Y_noisy_cov_rec,rh,n0);
%
% % Gradient descent
% [g_rec,f_rec,Outs] = WirtingerFlow(Y_noisy_cov_rec,Opts);
%
% figure; plot((2:Outs.Iter),log10(Outs.Fun(2:Outs.Iter)),'b*-');hold on;
% plot((2:Outs.Iter),ones(Outs.Iter-1,1)*log10(Opts.OptFun),'r');
% axis tight
% title('Objective function versus iteration in Wirtinger gradient descent')
% xlabel('Iteration')
% ylabel('Objective function in log_{10} scale')
%
% Gain_rec    = g_rec/norm(g_rec)*AlphaL2;
%
% % MUSIC
% ImagingFun   = MUSIC(f_rec,sparsity,ell);
% % Local Maximum extraction in MUSIC
% [Supp_rec , supprecind] = LM(ell,ImagingFun,RL,sparsity);
% % Compute error and adjust trivial ambiguities
% [Supp_rec_c2,Beta_rec_c12,Gain_rec_c,Supp_err,Alpha_reerr,Beta_reerr,Gain_reerr,c0,c1,c2,c] = ComputeReconstructionError(Supp,Supp_rec,Gain,Gain_rec,ell);
% fprintf('Support recovery error = %6.4f RL\n',Supp_err)
% fprintf('Amplitude error mean =%6.4f 0.9 0.5 0.1 quantitle =  %6.4f %6.4f %6.4f \n',mean(Alpha_reerr),quantile(Alpha_reerr,0.9),median(Alpha_reerr),quantile(Alpha_reerr,0.1) );
% fprintf('Phase     error mean =%6.4f 0.9 0.5 0.1 quantitle = %6.4f %6.4f %6.4f \n',mean(Beta_reerr),quantile(Beta_reerr,0.9),median(Beta_reerr),quantile(Beta_reerr,0.1) );
% fprintf('Cali      error mean =%6.4f 0.9 0.5 0.1 quantitle = %6.4f %6.4f %6.4f \n',mean(Gain_reerr),quantile(Gain_reerr,0.9),median(Gain_reerr),quantile(Gain_reerr,0.1) );
% fprintf('cali constant        = %6.2f + i %6.2f  magnitude =%6.2f\n',real(c),imag(c),abs(c))
%
% % Imaging function
% figure; plot(mod(Grid+c2/(2*pi),1),ImagingFun,'b*-'); hold on;
% stem(Supp,zeros(sparsity,1),'ro','filled');
% title('Imaging function in Wirtinger gradient descent')
% xlabel(['support error = ' num2str(Supp_err) ' RL'])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Method 2: Optimization with Witinger flow on the transformed problem
% display('======================Witinger flow 2=========================')
%
% % Initialization
% rh1           = n0;
% rh2           = n0;
% WHUk0         = [1./(GFk0(:,1)+1i*GFk0(:,2)) GFk0(:,3)+1i*GFk0(:,4)];
% Opts          = struct('MaxIter',500,'HUk0',WHUk0,'alphatol',10^(-4),'rh1',rh1,'rh2',rh2);
% [Opts.OptFun , Opts.OptFun_loss , Opts.OptFun_reg1 , Opts.OptFun_reg2]   = WFun2([1./Gain f],Y_noisy_cov_rec,rh1,rh2);
%
% % Gradient descent
% [h_rec,u_rec,Outs] = WirtingerFlow2(Y_noisy_cov_rec,Opts);
%
%
% figure; plot((2:Outs.Iter),log10(Outs.Fun_loss(2:Outs.Iter)),'b*-');hold on;
% plot((2:Outs.Iter),ones(Outs.Iter-1,1)*log10(Opts.OptFun_loss),'ro-');
% axis tight
% title('Loss function versus iteration in Wirtinger gradient descent 2')
% xlabel('Iteration')
% ylabel('Loss function in log_{10} scale')
% legend('Loss at each iteration','Loss at the truth')
%
% g_rec       = 1./h_rec;
% Gain_rec    = g_rec/norm(g_rec)*AlphaL2;
%
% % MUSIC
% ImagingFun   = MUSIC(u_rec,sparsity,ell);
% % Local Maximum extraction in MUSIC
% [Supp_rec , supprecind] = LM(ell,ImagingFun,RL,sparsity);
% % Compute error and adjust trivial ambiguities
% [Supp_rec_c2,Beta_rec_c12,Gain_rec_c,Supp_err,Alpha_reerr,Beta_reerr,Gain_reerr,c0,c1,c2,c] = ComputeReconstructionError(Supp,Supp_rec,Gain,Gain_rec,ell);
% fprintf('Support recovery error = %6.4f RL\n',Supp_err)
% fprintf('Amplitude error mean =%6.4f 0.9 0.5 0.1 quantitle =  %6.4f %6.4f %6.4f \n',mean(Alpha_reerr),quantile(Alpha_reerr,0.9),median(Alpha_reerr),quantile(Alpha_reerr,0.1) );
% fprintf('Phase     error mean =%6.4f 0.9 0.5 0.1 quantitle = %6.4f %6.4f %6.4f \n',mean(Beta_reerr),quantile(Beta_reerr,0.9),median(Beta_reerr),quantile(Beta_reerr,0.1) );
% fprintf('Cali      error mean =%6.4f 0.9 0.5 0.1 quantitle = %6.4f %6.4f %6.4f \n',mean(Gain_reerr),quantile(Gain_reerr,0.9),median(Gain_reerr),quantile(Gain_reerr,0.1) );
% fprintf('cali constant        = %6.2f + i %6.2f  magnitude =%6.2f\n',real(c),imag(c),abs(c))
%
% % Imaging function
% figure; plot(mod(Grid+c2/(2*pi),1),ImagingFun,'b*-'); hold on;
% stem(Supp,zeros(sparsity,1),'ro','filled');
% title('Imaging function in Wirtinger gradient descent')
% xlabel(['support error = ' num2str(Supp_err) ' RL'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Method 3: Alternating method by Friedlander and Weiss
% display('======================Friedlander and Weiss=========================')
% Opts          = struct('MaxIter',200,'g0',WGFk0(:,1),'epsilon',10^(-8));   % WGFk0(:,1) initialized by algebraic  rand(N,1)+1i*rand(N,1) random initialzation
% Opts.Jtrue    = FriedlanderWeiss_ComputeJtrue(Y_noisy_cov_rec,Gain,Supp);
% [g_rec,Supp_rec,Outs] = FriedlanderWeiss(Y_noisy_cov_rec,sparsity,ell,Opts);
%
% Gain_rec    = g_rec/norm(g_rec)*AlphaL2;
%
% [Supp_rec_c2,Beta_rec_c12,Gain_rec_c,Supp_err,Alpha_reerr,Beta_reerr,Gain_reerr,c0,c1,c2,c] = ComputeReconstructionError(Supp,Supp_rec,Gain,Gain_rec,ell);
% fprintf('Support recovery error = %6.4f RL\n',Supp_err)
% fprintf('Amplitude error mean =%6.4f 0.9 0.5 0.1 quantitle =  %6.4f %6.4f %6.4f \n',mean(Alpha_reerr),quantile(Alpha_reerr,0.9),median(Alpha_reerr),quantile(Alpha_reerr,0.1) );
% fprintf('Phase     error mean =%6.4f 0.9 0.5 0.1 quantitle = %6.4f %6.4f %6.4f \n',mean(Beta_reerr),quantile(Beta_reerr,0.9),median(Beta_reerr),quantile(Beta_reerr,0.1) );
% fprintf('Cali      error mean =%6.4f 0.9 0.5 0.1 quantitle = %6.4f %6.4f %6.4f \n',mean(Gain_reerr),quantile(Gain_reerr,0.9),median(Gain_reerr),quantile(Gain_reerr,0.1) );
% fprintf('cali constant        = %6.2f + i %6.2f  magnitude =%6.2f\n',real(c),imag(c),abs(c))



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Method 4: Trace minimization
% display('======================Trace minimization=========================')
% GainInv    = 1./Gain;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if nlevel == 0 % Call cvx for trace minimization
%     cvx_solver sdpt3%
%     cvx_begin sdp quiet
%     variable HH(N,N) complex hermitian
%     variable TT(N,N) hermitian toeplitz
%     minimize trace(HH)
%     subject to
%     Y_noisy_cov_rec.*HH == TT;
%     HH == hermitian_semidefinite(N);
%     TT == hermitian_semidefinite(N);
%     HH(2,1) == GainInv(2)*conj(GainInv(1));
%     cvx_end
%     [uu , hh]      = eigs(HH,2);
%     hh             = diag(hh);
%     Opts           = struct();
%     Opts.objrec    = trace(HH);
%     Opts.objtruth  = GainInv'*GainInv;
%     fprintf('objective at the truth = %6.4f objective at the recovery = %6.4f\n',Opts.objtruth,Opts.objrec)
%     fprintf('loss at the truth =%6.4f loss at the recovery %6.4f\n',norm(Y_noisy_cov_rec.*(GainInv*GainInv') - toeplitz(f'),'fro'),norm(Y_noisy_cov_rec.*HH - TT,'fro'))
%     fprintf('top eigenvalue = %6.4f   second eigenvalue = %6.4f\n',hh(1),hh(2));
% else % Call cvx for reguarlized trace minimization
%     %     lambda = 0.001;
%     %     cvx_solver sdpt3%
%     %     cvx_begin sdp quiet
%     %         variable HH(N,N) complex hermitian
%     %         variable TT(N,N) hermitian toeplitz
%     %         minimize square_pos(norm(Y_noisy_cov_rec.*HH - TT,'fro')) + lambda * trace(HH)
%     %         subject to
%     %             HH == hermitian_semidefinite(N);
%     %             TT == hermitian_semidefinite(N);
%     %             HH(2,1) == GainInv(2).*conj(GainInv(1));
%     %     cvx_end
%     lambda=0;
%     cvx_solver sdpt3%
%     cvx_begin sdp quiet
%     variable HH(N,N) complex hermitian
%     variable TT(N,N) hermitian toeplitz
%     minimize trace(HH)
%     subject to
%     HH == hermitian_semidefinite(N);
%     TT == hermitian_semidefinite(N);
%     square_pos(norm(Y_noisy_cov_rec.*HH - TT,'fro')) <= square_pos(norm(Y_noisy_cov_rec.*(GainInv*GainInv') - toeplitz(f'),'fro'))
%     HH(2,1) == GainInv(2).*conj(GainInv(1));
%     cvx_end
%     [uu , hh]        = eigs(HH,2);
%     hh               = diag(hh);
%     Opts             = struct();
%     Opts.objrec      = norm(Y_noisy_cov_rec.*HH - TT,'fro')+ lambda* trace(HH);
%     Opts.lossrec     = norm(Y_noisy_cov_rec.*HH - TT,'fro');
%     Opts.tracerec    = trace(HH);
%     Opts.objtruth    = norm(Y_noisy_cov_rec.*(GainInv*GainInv') - toeplitz(f'),'fro')+ lambda* (GainInv'*GainInv);
%     Opts.losstruth   = norm(Y_noisy_cov_rec.*(GainInv*GainInv') - toeplitz(f'),'fro');
%     Opts.tracetruth  = (GainInv'*GainInv);
%     fprintf('objective at the truth = %6.4f objective at the recovery = %6.4f\n',Opts.objtruth,Opts.objrec)
%     fprintf('loss at the truth = %6.4f loss at the recovery = %6.4f\n',Opts.losstruth,Opts.lossrec)
%     fprintf('trace at the truth = %6.4f trace at the recovery = %6.4f\n',Opts.tracetruth,Opts.tracerec)
%     fprintf('top eigenvalue = %6.4f   second eigenvalue = %6.4f\n',hh(1),hh(2));
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GainInv_rec         = uu(:,1);
% g_rec               = 1./GainInv_rec;
% Gain_rec            = g_rec/norm(g_rec)*AlphaL2;
% %% MUSIC and error
% ImagingFun   = MUSIC(TT,sparsity,ell);
% % Local Maximum extraction in MUSIC
% [Supp_rec , supprecind] = LM(ell,ImagingFun,RL,sparsity);                  % recovered support via our algorithm
% %Supp_rec_rootmusic      = rootmusic(TT,sparsity,'corr')/(2*pi);            % recovered support via rootmusic
% % Compute error and adjust trivial ambiguities
% [Supp_rec_c2,Beta_rec_c12,Gain_rec_c,Supp_err,Alpha_reerr,Beta_reerr,Gain_reerr,c0,c1,c2,c] = ComputeReconstructionError(Supp,Supp_rec,Gain,Gain_rec,ell);
% fprintf('Support recovery error = %6.4f RL\n',Supp_err)
% fprintf('Amplitude error mean =%6.4f 0.9 0.5 0.1 quantitle =  %6.4f %6.4f %6.4f \n',mean(Alpha_reerr),quantile(Alpha_reerr,0.9),median(Alpha_reerr),quantile(Alpha_reerr,0.1) );
% fprintf('Phase     error mean =%6.4f 0.9 0.5 0.1 quantitle = %6.4f %6.4f %6.4f \n',mean(Beta_reerr),quantile(Beta_reerr,0.9),median(Beta_reerr),quantile(Beta_reerr,0.1) );
% fprintf('Cali      error mean =%6.4f 0.9 0.5 0.1 quantitle = %6.4f %6.4f %6.4f \n',mean(Gain_reerr),quantile(Gain_reerr,0.9),median(Gain_reerr),quantile(Gain_reerr,0.1) );
% fprintf('cali constant        = %6.2f + i %6.2f  magnitude =%6.2f\n',real(c),imag(c),abs(c))
% 
% % Imaging function
% figure; plot(mod(Grid+c2/(2*pi),1),ImagingFun,'b*-'); hold on;
% stem(Supp,zeros(sparsity,1),'ro','filled');
% % stem(Supp_rec_rootmusic,zeros(sparsity,1),'g+','filled');
% title('Imaging function in trace minimization')
% xlabel(['support error = ' num2str(Supp_err) ' RL'])

%% Method 5: Modified Trace minimization (modified by MCho)
% display('======================modified Trace minimization=========================')
% GainInv    = 1./Gain;
% if nlevel == 0 % Call cvx for trace minimization
%     %    cvx_solver sdpt3
%     cvx_begin sdp quiet
%         variable HH(N,N) complex hermitian
%         variable u(N) complex
%         dual variable W1
%         dual variable W2
%         dual variable W3
%         
%         TT=toeplitz(u);
%         TT=0.5*(TT+TT'); % taking care of non-symmetric warning
%     minimize real(trace(TT))
%     subject to
%         TT - Y_noisy_cov_rec.*HH >= 0: W1;
%         HH >= 0: W2;
%         HH(2,1) == GainInv(2).*conj(GainInv(1)): W3;
%     % HH(2,1) == GainInv(2)*conj(GainInv(1));
%     cvx_end
%     [uu, hh]=eigs(HH,2);
%     hh=diag(hh);
%     Opts           = struct();
%     Opts.objrec    = trace(HH);
%     Opts.objtruth  = GainInv'*GainInv;
%     fprintf('objective at the truth = %6.4f objective at the recovery = %6.4f\n',Opts.objtruth,Opts.objrec)
%     fprintf('loss at the truth =%6.4f loss at the recovery %6.4f\n',norm(Y_noisy_cov_rec.*HH - TT,'fro'),norm(Y_noisy_cov_rec.*(GainInv*GainInv') - toeplitz(f'),'fro'))
%     fprintf('top eigenvalue = %6.4f   second eigenvalue = %6.4f\n',hh(1),hh(2));
% else % Call cvx for reguarlized trace minimization
%     lambda = 1;
%     cvx_solver sdpt3
%     cvx_begin sdp quiet
%     variable HH(N,N) complex hermitian
%     variable TT(N,N) hermitian toeplitz
%     minimize square_pos(norm(Y_noisy_cov_rec.*HH - TT,'fro')) + lambda * trace(TT)
%     subject to
%     TT - Y_noisy_cov_rec.*HH >= 0;
%     HH >= 0;
%     TT >= 0;
%     HH(2,1) == GainInv(2).*conj(GainInv(1));
%     cvx_end
%     [uu , hh]        = eigs(HH,2);
%     hh               = diag(hh);
%     Opts             = struct();
%     Opts.objrec      = norm(Y_noisy_cov_rec.*HH - TT,'fro')+ lambda* trace(HH);
%     Opts.lossrec     = norm(Y_noisy_cov_rec.*HH - TT,'fro');
%     Opts.tracerec    = trace(HH);
%     Opts.objtruth    = norm(Y_noisy_cov_rec.*(GainInv*GainInv') - toeplitz(f'),'fro')+ lambda* (GainInv'*GainInv);
%     Opts.losstruth   = norm(Y_noisy_cov_rec.*(GainInv*GainInv') - toeplitz(f'),'fro');
%     Opts.tracetruth  = (GainInv'*GainInv);
%     fprintf('objective at the truth = %6.4f objective at the recovery = %6.4f\n',Opts.objtruth,Opts.objrec)
%     fprintf('loss at the truth = %6.4f loss at the recovery = %6.4f\n',Opts.losstruth,Opts.lossrec)
%     fprintf('trace at the truth = %6.4f trace at the recovery = %6.4f\n',Opts.tracetruth,Opts.tracerec)
%     fprintf('top eigenvalue = %6.4f   second eigenvalue = %6.4f\n',hh(1),hh(2));
% end
% GainInv_rec         = uu(:,1);
% g_rec               = 1./GainInv_rec;
% Gain_rec            = g_rec/norm(g_rec)*AlphaL2;
% %% MUSIC and error
% ImagingFun   = MUSIC(TT,sparsity,ell);
% % Local Maximum extraction in MUSIC
% [Supp_rec , supprecind] = LM(ell,ImagingFun,RL,sparsity);                  % recovered support via our algorithm
% %Supp_rec_rootmusic      = rootmusic(TT,sparsity,'corr')/(2*pi);            % recovered support via rootmusic
% % Compute error and adjust trivial ambiguities
% [Supp_rec_c2,Beta_rec_c12,Gain_rec_c,Supp_err,Alpha_reerr,Beta_reerr,Gain_reerr,c0,c1,c2,c] = ComputeReconstructionError(Supp,Supp_rec,Gain,Gain_rec,ell);
% fprintf('Support recovery error = %6.4f RL\n',Supp_err)
% fprintf('Amplitude error mean =%6.4f 0.9 0.5 0.1 quantitle =  %6.4f %6.4f %6.4f \n',mean(Alpha_reerr),quantile(Alpha_reerr,0.9),median(Alpha_reerr),quantile(Alpha_reerr,0.1) );
% fprintf('Phase     error mean =%6.4f 0.9 0.5 0.1 quantitle = %6.4f %6.4f %6.4f \n',mean(Beta_reerr),quantile(Beta_reerr,0.9),median(Beta_reerr),quantile(Beta_reerr,0.1) );
% fprintf('Cali      error mean =%6.4f 0.9 0.5 0.1 quantitle = %6.4f %6.4f %6.4f \n',mean(Gain_reerr),quantile(Gain_reerr,0.9),median(Gain_reerr),quantile(Gain_reerr,0.1) );
% fprintf('cali constant        = %6.2f + i %6.2f  magnitude =%6.2f\n',real(c),imag(c),abs(c))
% 
% % Imaging function
% figure; plot(mod(Grid+c2/(2*pi),1),ImagingFun,'b*-'); hold on;
% stem(Supp,zeros(sparsity,1),'ro','filled');
% % stem(Supp_rec_rootmusic,zeros(sparsity,1),'g+','filled');
% title('Imaging function in trace minimization')
% xlabel(['support error = ' num2str(Supp_err) ' RL'])


%% Method 4: manifold minimization
display('======================Manifold minimization========================')

GainInv    = 1./Gain;
[ GainInv_rec,TT ] = Manopt_CalPara( Y_noisy_cov_rec );
g_rec=1./GainInv_rec;
Gain_rec=g_rec/norm(g_rec)*AlphaL2;

% MUSIC
ImagingFun   = MUSIC(TT,sparsity,ell);
% Local Maximum extraction in MUSIC
[Supp_rec , supprecind] = LM(ell,ImagingFun,RL,sparsity);
% Compute error and adjust trivial ambiguities
[Supp_rec_c2,Beta_rec_c12,Gain_rec_c,Supp_err,Alpha_reerr,Beta_reerr,Gain_reerr,c0,c1,c2,c] = ComputeReconstructionError(Supp,Supp_rec,Gain,Gain_rec,ell);
fprintf('Support recovery error = %6.4f RL\n',Supp_err)
fprintf('Amplitude error mean =%6.4f 0.9 0.5 0.1 quantitle =  %6.4f %6.4f %6.4f \n',mean(Alpha_reerr),quantile(Alpha_reerr,0.9),median(Alpha_reerr),quantile(Alpha_reerr,0.1) );
fprintf('Phase     error mean =%6.4f 0.9 0.5 0.1 quantitle = %6.4f %6.4f %6.4f \n',mean(Beta_reerr),quantile(Beta_reerr,0.9),median(Beta_reerr),quantile(Beta_reerr,0.1) );
fprintf('Cali      error mean =%6.4f 0.9 0.5 0.1 quantitle = %6.4f %6.4f %6.4f \n',mean(Gain_reerr),quantile(Gain_reerr,0.9),median(Gain_reerr),quantile(Gain_reerr,0.1) );
fprintf('cali constant        = %6.2f + i %6.2f  magnitude =%6.2f\n',real(c),imag(c),abs(c))

% Imaging function 
% figure; plot(mod(Grid+c2/(2*pi),1),ImagingFun,'b*-'); hold on;
% stem(Supp,zeros(sparsity,1),'ro','filled');
% title('Imaging function in Wirtinger gradient descent')
% x