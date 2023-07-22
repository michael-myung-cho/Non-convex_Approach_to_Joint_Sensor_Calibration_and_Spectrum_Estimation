%% Error versus # snapshot
clear all
close all

display('============================================================')


%% Set-up
sparsity = 20;                 % number of sources
stype    = 'off grid';        % source type: on gird - sources on the RL grid; off grid - sources on the continuum.
N        = 100;               % number of sensors
RL       = 1/N;               % Rayleigh length (resolution limit)
L        = 500;               % number of snapshot
srange   = 2;                 % range of source amplitudes
crange   = 2;                 % range of calibration amplitudes
fprintf('%6.0f sources %6.0f measurements\n',sparsity,N)

%% source location
sepn1               = 2;
sepn2               = 3;
sep1                = sepn1*RL;
sep2                = sepn2*RL;

Log10SIGMA            = -1:0.1:1;
SIGMA                 = 10.^Log10SIGMA;
NOISE                 = 100*SIGMA;
LenN                  = length(NOISE);
Times                 = 100;
competitors           = struct('noise',NOISE,'sigma',SIGMA,'alg',[],'col',[],'mar',[],...
                        'supperr',[],'supperr_mean',[],'alphaerr',[],'alphaerr_mean',[],'betaerr',[],'betaerr_mean',[]);
competitors.alg            = {'algebraic','algebraic full','opt','wirtinger','FriedlanderWeiss'};
competitors.supperr        = {zeros(Times,LenN),zeros(Times,LenN),zeros(Times,LenN),zeros(Times,LenN),zeros(Times,LenN)};
competitors.supperr_mean   = {zeros(1,LenN),zeros(1,LenN),zeros(1,LenN),zeros(1,LenN),zeros(1,LenN)};
competitors.alphaerr       = {zeros(Times,LenN),zeros(Times,LenN),zeros(Times,LenN),zeros(Times,LenN),zeros(Times,LenN)};
competitors.alphaerr_mean  = {zeros(1,LenN),zeros(1,LenN),zeros(1,LenN),zeros(1,LenN),zeros(1,LenN)};
competitors.betaerr        = {zeros(Times,LenN),zeros(Times,LenN),zeros(Times,LenN),zeros(Times,LenN),zeros(Times,LenN)};
competitors.betaerr_mean   = {zeros(1,LenN),zeros(1,LenN),zeros(1,LenN),zeros(1,LenN),zeros(1,LenN)};
competitors.gainerr        = {zeros(Times,LenN),zeros(Times,LenN),zeros(Times,LenN),zeros(Times,LenN),zeros(Times,LenN)};
competitors.gainerr_mean   = {zeros(1,LenN),zeros(1,LenN),zeros(1,LenN),zeros(1,LenN),zeros(1,LenN)};
competitors.col            = {'b','g','m','r','k'};
competitors.mar            = {'o','*','s','d','+'};


for inoise  = 1 : LenN
    nsigma = SIGMA(inoise);
    nlevel = nsigma*100;
    for itime = 1 : Times
        display('=================================================')
        fprintf('noise = %6.0f time = %6.0f\n',nlevel,itime)
        % generate source location
        [Supp , SuppIdx]    = GenerateFrequency(sparsity,sep1,sep2,stype,RL);
        % Generate source amplitudes
        Samp     = 1+rand(sparsity,1)*(srange-1);
        X        = zeros(sparsity,L);
        for j = 1:sparsity
            X(j,:) = Samp(j)*exp(2*pi*1i*rand(1,L));
        end
        % Sensor gain
        Alpha    = 1+rand(N,1)*(crange-1);
        Beta     = 2*pi*rand(N,1);
        Gain     = Alpha.*exp(1i*Beta);
        AlphaL2  = norm(Alpha);
        % Generate measurement
        A        = exp(2*pi*1i*((0:N-1)')*(Supp'));
        Y        = diag(Gain)*A*X/sqrt(N);
        f        = A*(Samp).^2/N;
        % Add noise
        % nsigma   = mean(mean(abs(Y)))*nlevel/100;
        Noise    = nsigma*exp(2*pi*1i*rand(N,L));
        % Noisey measurement
        Y_noisy      = Y + Noise;
        Y_noisy_cov  = Y_noisy*Y_noisy'/L;
        % Step 1: Estimate noise level
        [V , D]      = eig(Y_noisy_cov);
        D            = flipud(diag(D));
        V            = fliplr(V);
        nsigma_rec   = sqrt(mean(D(sparsity+1:end)));
        fprintf('exact noise sigma = %6.8f recovered = %6.8f\n',nsigma,nsigma_rec)
        Y_noisy_cov_rec   = Y_noisy_cov-nsigma_rec^2*eye(N);
        ell               = 1/50*RL;
        Grid              = 0 : ell : 1;
        GridLen           = length(Grid);
        %% Method 1: algebraic method
        display('------------------Algebraic method-------------------')
        % Reconstruction
        [Gain_rec ,  ImagingFun , f_rec] = AlgebraicMethod(Y_noisy_cov_rec,sparsity,ell,AlphaL2);
        n0                   = norm(Gain_rec)^2*norm(f_rec);
        Gain_rec_algebraic   = Gain_rec;
        f_rec_algebraic      = f_rec;
        % Local Maximum extraction in MUSIC
        [Supp_rec , supprecind] = LM(ell,ImagingFun,RL,sparsity);
        % Compute error and adjust trivial ambiguities
        [Supp_rec_c2,Beta_rec_c12,Gain_rec_c,Supp_err,Alpha_reerr,Beta_reerr,Gain_reerr,c0,c1,c2,c] = ComputeReconstructionError(Supp,Supp_rec,Gain,Gain_rec,ell);
        fprintf('Support recovery error = %6.4f RL\n',Supp_err)
        fprintf('Amplitude error mean =%6.4f 0.9 0.5 0.1 quantitle =  %6.4f %6.4f %6.4f \n',mean(Alpha_reerr),quantile(Alpha_reerr,0.9),median(Alpha_reerr),quantile(Alpha_reerr,0.1) );
        fprintf('Phase     error mean =%6.4f 0.9 0.5 0.1 quantitle = %6.4f %6.4f %6.4f \n',mean(Beta_reerr),quantile(Beta_reerr,0.9),median(Beta_reerr),quantile(Beta_reerr,0.1) );
        fprintf('Cali      error mean =%6.4f 0.9 0.5 0.1 quantitle = %6.4f %6.4f %6.4f \n',mean(Gain_reerr),quantile(Gain_reerr,0.9),median(Gain_reerr),quantile(Gain_reerr,0.1) );
        fprintf('cali constant        = %6.2f + i %6.2f  magnitude =%6.2f\n',real(c),imag(c),abs(c))
        competitors.supperr{1}(itime,inoise)  = Supp_err;
        competitors.alphaerr{1}(itime,inoise) = mean(Alpha_reerr);
        competitors.betaerr{1}(itime,inoise)  = mean(Beta_reerr);
        competitors.gainerr{1}(itime,inoise)  = mean(Gain_reerr);
        %% Method 2: full algebraic method
        display('--------------------Full Algebraic method-----------------')
        % Reconstruction
        [Gain_rec ,  ImagingFun] = AlgebraicMethodFull(Y_noisy_cov_rec,sparsity,ell,AlphaL2);
        % Local Maximum extraction in MUSIC
        [Supp_rec , supprecind] = LM(ell,ImagingFun,RL,sparsity);
        % Compute error and adjust trivial ambiguities
        [Supp_rec_c2,Beta_rec_c12,Gain_rec_c,Supp_err,Alpha_reerr,Beta_reerr,Gain_reerr,c0,c1,c2,c] = ComputeReconstructionError(Supp,Supp_rec,Gain,Gain_rec,ell);
        fprintf('Support recovery error = %6.4f RL\n',Supp_err)
        fprintf('Amplitude error mean =%6.4f 0.9 0.5 0.1 quantitle =  %6.4f %6.4f %6.4f \n',mean(Alpha_reerr),quantile(Alpha_reerr,0.9),median(Alpha_reerr),quantile(Alpha_reerr,0.1) );
        fprintf('Phase     error mean =%6.4f 0.9 0.5 0.1 quantitle = %6.4f %6.4f %6.4f \n',mean(Beta_reerr),quantile(Beta_reerr,0.9),median(Beta_reerr),quantile(Beta_reerr,0.1) );
        fprintf('Cali      error mean =%6.4f 0.9 0.5 0.1 quantitle = %6.4f %6.4f %6.4f \n',mean(Gain_reerr),quantile(Gain_reerr,0.9),median(Gain_reerr),quantile(Gain_reerr,0.1) );
        fprintf('cali constant        = %6.2f + i %6.2f  magnitude =%6.2f\n',real(c),imag(c),abs(c))
        competitors.supperr{2}(itime,inoise)  = Supp_err;
        competitors.alphaerr{2}(itime,inoise) = mean(Alpha_reerr);
        competitors.betaerr{2}(itime,inoise)  = mean(Beta_reerr);
        competitors.gainerr{2}(itime,inoise)  = mean(Gain_reerr);
        % Method 3: Optimization
        display('--------------------Optimization method--------------------')
        % Initialization
        GFk0          = rand(N,4);
        GFk0(:,1)     = real(n0^(1/4)*Gain_rec_algebraic/norm(Gain_rec_algebraic)); %GFk0(:,1:2)*crange/sqrt(2);
        GFk0(:,2)     = imag(n0^(1/4)*Gain_rec_algebraic/norm(Gain_rec_algebraic)); 
        GFk0(:,3)     = real(sqrt(n0)*f_rec_algebraic/norm(f_rec_algebraic));      %GFk0(:,1:2)*srange^2/sqrt(2)*sparsity/N;
        GFk0(:,4)     = imag(sqrt(n0)*f_rec_algebraic/norm(f_rec_algebraic));  
        Opts          = struct('MaxIter',500,'GFk0',GFk0,'alphatol',10^(-4));
        Opts.OptFun   = Fun([real(Gain) imag(Gain) real(f) imag(f)],Y_noisy_cov_rec);
        % Gradient descent
        [g_rec,f_rec,Outs] = OptimizationMethod(Y_noisy_cov_rec,Opts);
        Alpha_rec    = abs(g_rec);
        Alpha_rec    = Alpha_rec/norm(Alpha_rec)*AlphaL2;
        Beta_rec     = angle(g_rec);
        Gain_rec     = Alpha_rec.*exp(1i*Beta_rec);
        % MUSIC
        ImagingFun   = MUSIC(f_rec,sparsity,ell);
        % Local Maximum extraction in MUSIC
        [Supp_rec , supprecind] = LM(ell,ImagingFun,RL,sparsity);
        % Compute error and adjust trivial ambiguities
        [Supp_rec_c2,Beta_rec_c12,Gain_rec_c,Supp_err,Alpha_reerr,Beta_reerr,Gain_reerr,c0,c1,c2,c] = ComputeReconstructionError(Supp,Supp_rec,Gain,Gain_rec,ell);
        fprintf('Support recovery error = %6.4f RL\n',Supp_err)
        fprintf('Amplitude error mean =%6.4f 0.9 0.5 0.1 quantitle =  %6.4f %6.4f %6.4f \n',mean(Alpha_reerr),quantile(Alpha_reerr,0.9),median(Alpha_reerr),quantile(Alpha_reerr,0.1) );
        fprintf('Phase     error mean =%6.4f 0.9 0.5 0.1 quantitle = %6.4f %6.4f %6.4f \n',mean(Beta_reerr),quantile(Beta_reerr,0.9),median(Beta_reerr),quantile(Beta_reerr,0.1) );
        fprintf('Cali      error mean =%6.4f 0.9 0.5 0.1 quantitle = %6.4f %6.4f %6.4f \n',mean(Gain_reerr),quantile(Gain_reerr,0.9),median(Gain_reerr),quantile(Gain_reerr,0.1) );
        fprintf('cali constant        = %6.2f + i %6.2f  magnitude =%6.2f\n',real(c),imag(c),abs(c))
        competitors.supperr{3}(itime,inoise)  = Supp_err;
        competitors.alphaerr{3}(itime,inoise) = mean(Alpha_reerr);
        competitors.betaerr{3}(itime,inoise)  = mean(Beta_reerr);
        competitors.gainerr{3}(itime,inoise)  = mean(Gain_reerr);
        %% Method 3: Optimization with Witinger flow
        display('--------------------Witinger flow-----------------------')
        % Initialization
        rh            = n0;
        WGFk0         = [GFk0(:,1)+1i*GFk0(:,2) GFk0(:,3)+1i*GFk0(:,4)];
        Opts          = struct('MaxIter',500,'GFk0',WGFk0,'alphatol',10^(-4),'rh',rh,'n0',n0);
        [Opts.OptFun , Opts.OptFun_loss , Opts.OptFun_reg]   = WFun([Gain f],Y_noisy_cov_rec,rh,n0);
        % Gradient descent
        [g_rec,f_rec,Outs] = WirtingerFlow(Y_noisy_cov_rec,Opts);
        Alpha_rec    = abs(g_rec);
        Alpha_rec    = Alpha_rec/norm(Alpha_rec)*AlphaL2;
        Beta_rec     = angle(g_rec);
        Gain_rec     = Alpha_rec.*exp(1i*Beta_rec);
        % MUSIC
        ImagingFun   = MUSIC(f_rec,sparsity,ell);
        % Local Maximum extraction in MUSIC
        [Supp_rec , supprecind] = LM(ell,ImagingFun,RL,sparsity);
        % Compute error and adjust trivial ambiguities
        [Supp_rec_c2,Beta_rec_c12,Gain_rec_c,Supp_err,Alpha_reerr,Beta_reerr,Gain_reerr,c0,c1,c2,c] = ComputeReconstructionError(Supp,Supp_rec,Gain,Gain_rec,ell);
        fprintf('Support recovery error = %6.4f RL\n',Supp_err)
        fprintf('Amplitude error mean =%6.4f 0.9 0.5 0.1 quantitle =  %6.4f %6.4f %6.4f \n',mean(Alpha_reerr),quantile(Alpha_reerr,0.9),median(Alpha_reerr),quantile(Alpha_reerr,0.1) );
        fprintf('Phase     error mean =%6.4f 0.9 0.5 0.1 quantitle = %6.4f %6.4f %6.4f \n',mean(Beta_reerr),quantile(Beta_reerr,0.9),median(Beta_reerr),quantile(Beta_reerr,0.1) );
        fprintf('Cali      error mean =%6.4f 0.9 0.5 0.1 quantitle = %6.4f %6.4f %6.4f \n',mean(Gain_reerr),quantile(Gain_reerr,0.9),median(Gain_reerr),quantile(Gain_reerr,0.1) );
        fprintf('cali constant        = %6.2f + i %6.2f  magnitude =%6.2f\n',real(c),imag(c),abs(c))
        competitors.supperr{4}(itime,inoise)  = Supp_err;
        competitors.alphaerr{4}(itime,inoise) = mean(Alpha_reerr);
        competitors.betaerr{4}(itime,inoise)  = mean(Beta_reerr);
        competitors.gainerr{4}(itime,inoise)  = mean(Gain_reerr);
        %% Method 5: Alternating method by Friedlander and Weiss
        display('--------------------Friedlander and Weiss-----------------------')        
        % Opts          = struct('MaxIter',500,'g0',WGFk0(:,1),'epsilon',10^(-4));
        Opts          = struct('MaxIter',500,'g0',rand(N,1)+1i*rand(N,1),'epsilon',10^(-4));
        Opts.Jtrue    = FriedlanderWeiss_ComputeJtrue(Y_noisy_cov_rec,Gain,Supp);
        [g_rec,Supp_rec,~] = FriedlanderWeiss(Y_noisy_cov_rec,sparsity,ell,Opts);
        Gain_rec      = g_rec/norm(g_rec)*AlphaL2;
        [Supp_rec_c2,Beta_rec_c12,Gain_rec_c,Supp_err,Alpha_reerr,Beta_reerr,Gain_reerr,c0,c1,c2,c] = ComputeReconstructionError(Supp,Supp_rec,Gain,Gain_rec,ell);
        fprintf('Support recovery error = %6.4f RL\n',Supp_err)
        fprintf('Amplitude error mean =%6.4f 0.9 0.5 0.1 quantitle =  %6.4f %6.4f %6.4f \n',mean(Alpha_reerr),quantile(Alpha_reerr,0.9),median(Alpha_reerr),quantile(Alpha_reerr,0.1) );
        fprintf('Phase     error mean =%6.4f 0.9 0.5 0.1 quantitle = %6.4f %6.4f %6.4f \n',mean(Beta_reerr),quantile(Beta_reerr,0.9),median(Beta_reerr),quantile(Beta_reerr,0.1) );
        fprintf('Cali      error mean =%6.4f 0.9 0.5 0.1 quantitle = %6.4f %6.4f %6.4f \n',mean(Gain_reerr),quantile(Gain_reerr,0.9),median(Gain_reerr),quantile(Gain_reerr,0.1) );
        fprintf('cali constant        = %6.2f + i %6.2f  magnitude =%6.2f\n',real(c),imag(c),abs(c))
        competitors.supperr{5}(itime,inoise)  = Supp_err;
        competitors.alphaerr{5}(itime,inoise) = mean(Alpha_reerr);
        competitors.betaerr{5}(itime,inoise)  = mean(Beta_reerr);
        competitors.gainerr{5}(itime,inoise)  = mean(Gain_reerr);
    end
end

for ialgorithm = 1:5
    competitors.supperr_mean{ialgorithm}  = squeeze(mean(competitors.supperr{ialgorithm},1));
    competitors.alphaerr_mean{ialgorithm} = squeeze(mean(competitors.alphaerr{ialgorithm},1));
    competitors.betaerr_mean{ialgorithm}  = squeeze(mean(competitors.betaerr{ialgorithm},1));
    competitors.gainerr_mean{ialgorithm}  = squeeze(mean(competitors.gainerr{ialgorithm},1));
end

save(['ErrorVersusNoiseL' num2str(L),'.mat'],'NOISE','competitors','sparsity','N','srange','crange','nlevel','sepn1','sepn2','ell');



