% Main code: source localization + sensor calibration with multiple snapshot of measurements

clear all
close all
clc

% save data to xls file
countXLS=1;
clk=clock;
datestmp=strcat([date,'-',num2str(clk(4)),'h',num2str(clk(5))]);
filename = strcat('N20_k5_Prob-Sensor_Cali-',datestmp,'.xlsx');
col_name = {'N','sparsity','RL','nLevel','L','nSuccAlge','nSuccOptWF','nSuccAltFW','nSuccPANM';};
xlswrite(filename,col_name);
resBuf=[];

nTri=1;
for nlevel=200
    for ii=1:7
        L=10^ii;
        
        nSuccPANM=0;
        nSuccAltFW=0;
        nSuccAlge=0;
        nSuccOptWF=0;
        nSuccRPANM=0;
        nSuccNON=0;
        
        for itri=1:nTri
            display('============================================================')
            % Set-up
            sparsity = 5;                 % number of sources
            stype    = 'off grid';        % source type: on gird - sources on the RL grid; off grid - sources on the continuum.
            N        = 20;               % number of sensors
            RL       = 1/N;               % Rayleigh length (resolution limit)
            L        = 10^4;               % number of snapshot
            srange   = 2;                 % range of source amplitudes
            crange   = 2;                 % range of calibration amplitudes
            nlevel   = 0;                 % percentage of noise
            fprintf('%6.0f sources %6.0f measurements and %6.0f snapshot  ',sparsity,N,L)
            fprintf('noise sigma = %6.0f \n',nlevel/100)
            
            % Generate source location
            sepn1               = 2;
            sepn2               = 3;
            sep1                = sepn1*RL;
            sep2                = sepn2*RL;
            [Supp , SuppIdx]    = GenerateFrequency(sparsity,sep1,sep2,stype,RL);
            
            % Generate source amplitudes
            Samp     = 1+rand(sparsity,1)*(srange-1);
            X        = zeros(sparsity,L);
            for j = 1:sparsity
                X(j,:) = Samp(j)*exp(2*pi*1i*rand(1,L));
            end
            
            % Correlated source generation (MCho)
            Samp     = 1+rand(sparsity,1)*(srange-1);
            X        = zeros(sparsity,L);
            a=rand(1,L);
            for j = 1:sparsity
                X(j,:) = Samp(j)*exp(2*pi*1i*a);
            end
            
            
            X_cov = X*X'/L;    % covariance matrix of source
            
            % Sensor gain
            Alpha    = 1+rand(N,1)*(crange-1);
            Beta     = 2*pi*rand(N,1);
            Gain     = Alpha.*exp(1i*Beta);
            AlphaL2  = norm(Alpha);
            
            
            % Generate measurement
            A        = exp(2*pi*1i*((0:N-1)')*(Supp'));
            Y        = diag(Gain)*A*X/sqrt(N);
            f        = A*(Samp).^2/N;
            
            Ttemp=A*X*X'*A'/L;
            % Add noise
            nsigma   = nlevel/100;  %mean(mean(abs(Y)))*nlevel/100;
            Noise    = nsigma*exp(2*pi*1i*rand(N,L));
            
            Noise_cov = Noise*Noise'/L;
            
            
            % Noisey measurement
            Y_noisy      = Y + Noise;
            Y_noisy_cov  = Y_noisy*Y_noisy'/L;
            X_cov        = (A*X)*(A*X)'/(N*L);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Step 1: Estimate noise level
            [V , D]      = eig(Y_noisy_cov);
            D            = flipud(diag(D));
            V            = fliplr(V);
            
            nsigma_rec   = sqrt(mean(D(sparsity+1:end)));
            fprintf('exact noise sigma = %6.8f recovered = %6.8f\n',nsigma,nsigma_rec)
            
            Y_noisy_cov_rec   = Y_noisy_cov-nsigma_rec^2*eye(N);
            
            ell               = 1/100*RL;
            Grid              = 0 : ell : 1;
            GridLen           = length(Grid);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Method 1: algebraic method
            display('===================Algebraic method=========================')
            
            Reconstruction
            [Gain_rec ,  ImagingFun , f_rec] = AlgebraicMethod(Y_noisy_cov_rec,sparsity,ell,AlphaL2);
            
            n0                   = norm(Gain_rec)^2*norm(f_rec);
            Gain_rec_algebraic   = Gain_rec;
            f_rec_algebraic      = f_rec;
            
            Local Maximum extraction in MUSIC
            [Supp_rec , supprecind] = LM(ell,ImagingFun,RL,sparsity);
            Compute error and adjust trivial ambiguities
            [Supp_rec_c2,Beta_rec_c12,Gain_rec_c,Supp_err,Alpha_reerr,Beta_reerr,Gain_reerr,c0,c1,c2,c] = ComputeReconstructionError(Supp,Supp_rec,Gain,Gain_rec,ell);
            fprintf('Support recovery error = %6.4f RL\n',Supp_err)
            fprintf('Amplitude error mean =%6.4f 0.9 0.5 0.1 quantitle =  %6.4f %6.4f %6.4f \n',mean(Alpha_reerr),quantile(Alpha_reerr,0.9),median(Alpha_reerr),quantile(Alpha_reerr,0.1) );
            fprintf('Phase     error mean =%6.4f 0.9 0.5 0.1 quantitle = %6.4f %6.4f %6.4f \n',mean(Beta_reerr),quantile(Beta_reerr,0.9),median(Beta_reerr),quantile(Beta_reerr,0.1) );
            fprintf('Cali      error mean =%6.4f 0.9 0.5 0.1 quantitle = %6.4f %6.4f %6.4f \n',mean(Gain_reerr),quantile(Gain_reerr,0.9),median(Gain_reerr),quantile(Gain_reerr,0.1) );
            fprintf('cali constant        = %6.2f + i %6.2f  magnitude =%6.2f\n',real(c),imag(c),abs(c))
            Imaging function
            figure;  hold on;
            plot(mod(Grid+c2/(2*pi),1),ImagingFun,'b*-'); hold on;
            stem(Supp,zeros(sparsity,1),'ro','filled');
            title('Imaging function in the algebraic method')
            xlabel(['support error = ' num2str(Supp_err) ' RL'])
            
            
            % tracking # of success
            if Supp_err < 0.2/N
                nSuccAlge=nSuccAlge+1;
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Method 2: full algebraic method
            display('===================Full Algebraic method=========================')
            
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
            % Imaging function
            figure; stem(Supp,Samp,'ro'); hold on;
            plot(Grid,ImagingFun,'b*-')
            title('Imaging function for the algebraic method')
            xlabel(['support error = ' num2str(Supp_err) ' RL'])
            
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Method 3: Optimization
            display('===================Optimization method with real variables======================')
            
            Initialization in optimization is given by the algebraic method
            GFk0          = rand(N,4);
            GFk0(:,1)     = real(n0^(1/4)*Gain_rec_algebraic/norm(Gain_rec_algebraic)); %GFk0(:,1:2)*crange/sqrt(2);
            GFk0(:,2)     = imag(n0^(1/4)*Gain_rec_algebraic/norm(Gain_rec_algebraic));
            GFk0(:,3)     = real(sqrt(n0)*f_rec_algebraic/norm(f_rec_algebraic));      %GFk0(:,1:2)*srange^2/sqrt(2)*sparsity/N;
            GFk0(:,4)     = imag(sqrt(n0)*f_rec_algebraic/norm(f_rec_algebraic));
            Opts          = struct('MaxIter',500,'GFk0',GFk0,'alphatol',10^(-8));
            Opts.OptFun   = Fun([real(Gain) imag(Gain) real(f) imag(f)],Y_noisy_cov_rec);
            
            % Gradient descent
            [g_rec,f_rec,Outs] = OptimizationMethod(Y_noisy_cov_rec,Opts);
            
            figure; plot((2:Outs.Iter),log10(Outs.Fun(2:Outs.Iter)),'b*-');hold on;
            plot((2:Outs.Iter),ones(Outs.Iter-1,1)*log10(Opts.OptFun),'r');
            axis tight
            title('Objective function versus iteration in optimization')
            xlabel('Iteration')
            ylabel('Objective function in log_{10} scale')
            
            Gain_rec    = g_rec/norm(g_rec)*AlphaL2;
            
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
            
            Imaging function
            figure; plot(mod(Grid+c2/(2*pi),1),ImagingFun,'b*-'); hold on;
            stem(Supp,zeros(sparsity,1),'ro','filled');
            title('Imaging function in the optimization approach')
            xlabel(['support error = ' num2str(Supp_err) ' RL'])
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Method 2: Optimization with Witinger flow
            display('======================Witinger flow=========================')
            
            Initialization
            rh            = 0;
            WGFk0         = [GFk0(:,1)+1i*GFk0(:,2) GFk0(:,3)+1i*GFk0(:,4)];
            Opts          = struct('MaxIter',500,'GFk0',WGFk0,'alphatol',10^(-4),'rh',rh,'n0',n0);
            [Opts.OptFun , Opts.OptFun_loss , Opts.OptFun_reg]   = WFun([Gain f],Y_noisy_cov_rec,rh,n0);
            
            Gradient descent
            [g_rec,f_rec,Outs] = WirtingerFlow(Y_noisy_cov_rec,Opts);
            
            figure; plot((2:Outs.Iter),log10(Outs.Fun(2:Outs.Iter)),'b*-');hold on;
            plot((2:Outs.Iter),ones(Outs.Iter-1,1)*log10(Opts.OptFun),'r');
            axis tight
            title('Objective function versus iteration in Wirtinger gradient descent')
            xlabel('Iteration')
            ylabel('Objective function in log_{10} scale')
            
            Gain_rec    = g_rec/norm(g_rec)*AlphaL2;
            
            MUSIC
            ImagingFun   = MUSIC(f_rec,sparsity,ell);
            Local Maximum extraction in MUSIC
            [Supp_rec , supprecind] = LM(ell,ImagingFun,RL,sparsity);
            Compute error and adjust trivial ambiguities
            [Supp_rec_c2,Beta_rec_c12,Gain_rec_c,Supp_err,Alpha_reerr,Beta_reerr,Gain_reerr,c0,c1,c2,c] = ComputeReconstructionError(Supp,Supp_rec,Gain,Gain_rec,ell);
            fprintf('Support recovery error = %6.4f RL\n',Supp_err)
            fprintf('Amplitude error mean =%6.4f 0.9 0.5 0.1 quantitle =  %6.4f %6.4f %6.4f \n',mean(Alpha_reerr),quantile(Alpha_reerr,0.9),median(Alpha_reerr),quantile(Alpha_reerr,0.1) );
            fprintf('Phase     error mean =%6.4f 0.9 0.5 0.1 quantitle = %6.4f %6.4f %6.4f \n',mean(Beta_reerr),quantile(Beta_reerr,0.9),median(Beta_reerr),quantile(Beta_reerr,0.1) );
            fprintf('Cali      error mean =%6.4f 0.9 0.5 0.1 quantitle = %6.4f %6.4f %6.4f \n',mean(Gain_reerr),quantile(Gain_reerr,0.9),median(Gain_reerr),quantile(Gain_reerr,0.1) );
            fprintf('cali constant        = %6.2f + i %6.2f  magnitude =%6.2f\n',real(c),imag(c),abs(c))
            
            Imaging function
            figure; plot(mod(Grid+c2/(2*pi),1),ImagingFun,'b*-'); hold on;
            stem(Supp,zeros(sparsity,1),'ro','filled');
            title('Imaging function in Wirtinger gradient descent')
            xlabel(['support error = ' num2str(Supp_err) ' RL'])
            
            % tracking # of success
            if Supp_err < 0.2/N
                nSuccOptWF=nSuccOptWF+1;
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Method 2: Optimization with Witinger flow on the transformed problem
            display('======================Witinger flow 2=========================')
            
            % Initialization
            rh1           = n0;
            rh2           = n0;
            WHUk0         = [1./(GFk0(:,1)+1i*GFk0(:,2)) GFk0(:,3)+1i*GFk0(:,4)];
            Opts          = struct('MaxIter',500,'HUk0',WHUk0,'alphatol',10^(-4),'rh1',rh1,'rh2',rh2);
            [Opts.OptFun , Opts.OptFun_loss , Opts.OptFun_reg1 , Opts.OptFun_reg2]   = WFun2([1./Gain f],Y_noisy_cov_rec,rh1,rh2);
            
            % Gradient descent
            [h_rec,u_rec,Outs] = WirtingerFlow2(Y_noisy_cov_rec,Opts);
            
            
            % figure; plot((2:Outs.Iter),log10(Outs.Fun_loss(2:Outs.Iter)),'b*-');hold on;
            % plot((2:Outs.Iter),ones(Outs.Iter-1,1)*log10(Opts.OptFun_loss),'ro-');
            % axis tight
            % title('Loss function versus iteration in Wirtinger gradient descent 2')
            % xlabel('Iteration')
            % ylabel('Loss function in log_{10} scale')
            % legend('Loss at each iteration','Loss at the truth')
            
            g_rec       = 1./h_rec;
            Gain_rec    = g_rec/norm(g_rec)*AlphaL2;
            
            % MUSIC
            ImagingFun   = MUSIC(u_rec,sparsity,ell);
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
            % xlabel(['support error = ' num2str(Supp_err) ' RL'])
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Method 3: Alternating method by Friedlander and Weiss
            display('======================Friedlander and Weiss=========================')
            Opts          = struct('MaxIter',200,'g0',rand(N,1)+1i*rand(N,1),'epsilon',10^(-4));
            Opts.Jtrue    = FriedlanderWeiss_ComputeJtrue(Y_noisy_cov_rec,Gain,Supp);
            [g_rec,Supp_rec,Outs] = FriedlanderWeiss(Y_noisy_cov_rec,sparsity,ell,Opts);
            
            Gain_rec    = g_rec/norm(g_rec)*AlphaL2;
            
            [Supp_rec_c2,Beta_rec_c12,Gain_rec_c,Supp_err,Alpha_reerr,Beta_reerr,Gain_reerr,c0,c1,c2,c] = ComputeReconstructionError(Supp,Supp_rec,Gain,Gain_rec,ell);
            fprintf('Support recovery error = %6.4f RL\n',Supp_err)
            fprintf('Amplitude error mean =%6.4f 0.9 0.5 0.1 quantitle =  %6.4f %6.4f %6.4f \n',mean(Alpha_reerr),quantile(Alpha_reerr,0.9),median(Alpha_reerr),quantile(Alpha_reerr,0.1) );
            fprintf('Phase     error mean =%6.4f 0.9 0.5 0.1 quantitle = %6.4f %6.4f %6.4f \n',mean(Beta_reerr),quantile(Beta_reerr,0.9),median(Beta_reerr),quantile(Beta_reerr,0.1) );
            fprintf('Cali      error mean =%6.4f 0.9 0.5 0.1 quantitle = %6.4f %6.4f %6.4f \n',mean(Gain_reerr),quantile(Gain_reerr,0.9),median(Gain_reerr),quantile(Gain_reerr,0.1) );
            fprintf('cali constant        = %6.2f + i %6.2f  magnitude =%6.2f\n',real(c),imag(c),abs(c))
            
            % tracking # of success
            if Supp_err < 0.2/N
                nSuccAltFW=nSuccAltFW+1;
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Method 4: PANM in MMV
            display('======================PANM in MMV========================')
            GainInv=1./Gain;
            try
                
                cvx_begin sdp quiet
                variable HH(N,N) complex hermitian
                variable u(N) complex
                TT=toeplitz(u);
                TT=0.5*(TT+TT'); % taking care of non-symmetric warning
                minimize real(trace(TT))
                subject to
                TT - Y_noisy_cov_rec.*HH >= 0;
                HH >= 0;
                HH(2,1) == 10;
                cvx_end
                [uu, hh]=eigs(HH,2);
                hh=diag(hh);
                
                GainInv_rec=uu(:,1);
                g_rec=1./GainInv_rec;
                Gain_rec=g_rec/norm(g_rec)*AlphaL2;
                % MUSIC and error
                ImagingFun=MUSIC(TT,sparsity,ell);
                Local Maximum extraction in MUSIC
                [Supp_rec1 , supprecind]=LM(ell,ImagingFun,RL,sparsity);                  % recovered support via our algorithm
                Supp_rec1      = rootmusic(TT,sparsity,'corr')/(2*pi);            % recovered support via rootmusic
                Compute error and adjust trivial ambiguities
                [Supp_rec_c2,Beta_rec_c12,Gain_rec_c,Supp_err1,Alpha_reerr,Beta_reerr,Gain_reerr,c0,c1,c2,c] = ComputeReconstructionError(Supp,Supp_rec1,Gain,Gain_rec,ell);
                fprintf('Support recovery error = %6.4f RL\n',Supp_err1)
                fprintf('Amplitude error mean =%6.4f 0.9 0.5 0.1 quantitle =  %6.4f %6.4f %6.4f \n',mean(Alpha_reerr),quantile(Alpha_reerr,0.9),median(Alpha_reerr),quantile(Alpha_reerr,0.1) );
                fprintf('Phase     error mean =%6.4f 0.9 0.5 0.1 quantitle = %6.4f %6.4f %6.4f \n',mean(Beta_reerr),quantile(Beta_reerr,0.9),median(Beta_reerr),quantile(Beta_reerr,0.1) );
                fprintf('Cali      error mean =%6.4f 0.9 0.5 0.1 quantitle = %6.4f %6.4f %6.4f \n',mean(Gain_reerr),quantile(Gain_reerr,0.9),median(Gain_reerr),quantile(Gain_reerr,0.1) );
                fprintf('cali constant        = %6.2f + i %6.2f  magnitude =%6.2f\n',real(c),imag(c),abs(c))
                
                % Imaging function
                figure; plot(mod(Grid+c2/(2*pi),1),ImagingFun,'b*-'); hold on;
                stem(Supp,zeros(sparsity,1),'ro','filled');
                % stem(Supp_rec_rootmusic,zeros(sparsity,1),'g+','filled');
                title('Imaging function in trace minimization');
                xlabel(['support error = ' num2str(Supp_err) ' RL']);
                
                % tracking # of success
                if Supp_err1 < 0.2/N
                    nSuccPANM=nSuccPANM+1;
                end
            catch
                Error occures do nothing!
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Method 5: Reweighted PANM in MMV
            display('======================RPANM in MMV========================')
            GainInv=1./Gain;
            try
                
                uOld=zeros(N,1);
                epsilon=10^-6;
                lambda=10;
                for itrRPANM=1:10
                    W=inv(toeplitz(uOld) + epsilon*eye(N));
                    cvx_begin sdp quiet
                    variable HH2(N,N) complex hermitian
                    variable u2(N) complex
                    TT2=toeplitz(u2);
                    TT2=0.5*(TT2+TT2'); % taking care of non-symmetric warning
                    minimize real(trace(W*TT2))
                    subject to
                    TT2 - Y_noisy_cov_rec.*HH2 >= 0;
                    HH2 >= 0;
                    HH2(2,1) == 10;
                    cvx_end
                    [uu, hh]=eigs(HH2,2);
                    hh=diag(hh);
                    
                    GainInv_rec=uu(:,1);
                    g_rec=1./GainInv_rec;
                    Gain_rec2=g_rec/norm(g_rec)*AlphaL2;
                    ImagingFun=MUSIC(TT2,sparsity,ell);
                    Local Maximum extraction in MUSIC
                    [Supp_rec2 , supprecind]=LM(ell,ImagingFun,RL,sparsity);                  % recovered support via our algorithm
                    Supp_rec2      = rootmusic(TT2,sparsity,'corr')/(2*pi);            % recovered support via rootmusic
                    Compute error and adjust trivial ambiguities
                    [Supp_rec_c2,Beta_rec_c12,Gain_rec_c,Supp_err2,Alpha_reerr,Beta_reerr,Gain_reerr,c0,c1,c2,c] = ComputeReconstructionError(Supp,Supp_rec2,Gain,Gain_rec2,ell);
                    if norm(HH2(2,1) - 10) > 1
                        TT2=TTold;
                        HH2=HHold;
                        break;
                    end
                    if norm(uOld-u2)/norm(u2)< 10^-5
                        break;
                    end
                    HHold=HH2;
                    TTold=TT2;
                    uOld=u2;
                end
                
                [uu, hh]=eigs(HH2,2);
                hh=diag(hh);
                GainInv_rec=uu(:,1);
                g_rec=1./GainInv_rec;
                Gain_rec=g_rec/norm(g_rec)*AlphaL2;
                
                
                % MUSIC and error
                ImagingFun=MUSIC(TT2,sparsity,ell);
                Local Maximum extraction in MUSIC
                [Supp_rec2 , supprecind]=LM(ell,ImagingFun,RL,sparsity);                  % recovered support via our algorithm
                Supp_rec2      = rootmusic(TT2,sparsity,'corr')/(2*pi);            % recovered support via rootmusic
                Compute error and adjust trivial ambiguities
                [Supp_rec_c2,Beta_rec_c12,Gain_rec_c,Supp_err2,Alpha_reerr,Beta_reerr,Gain_reerr,c0,c1,c2,c] = ComputeReconstructionError(Supp,Supp_rec2,Gain,Gain_rec,ell);
                fprintf('Support recovery error = %6.4f RL\n',Supp_err2)
                fprintf('Amplitude error mean =%6.4f 0.9 0.5 0.1 quantitle =  %6.4f %6.4f %6.4f \n',mean(Alpha_reerr),quantile(Alpha_reerr,0.9),median(Alpha_reerr),quantile(Alpha_reerr,0.1) );
                fprintf('Phase     error mean =%6.4f 0.9 0.5 0.1 quantitle = %6.4f %6.4f %6.4f \n',mean(Beta_reerr),quantile(Beta_reerr,0.9),median(Beta_reerr),quantile(Beta_reerr,0.1) );
                fprintf('Cali      error mean =%6.4f 0.9 0.5 0.1 quantitle = %6.4f %6.4f %6.4f \n',mean(Gain_reerr),quantile(Gain_reerr,0.9),median(Gain_reerr),quantile(Gain_reerr,0.1) );
                fprintf('cali constant        = %6.2f + i %6.2f  magnitude =%6.2f\n',real(c),imag(c),abs(c))
                
                % tracking # of success
                if Supp_err2 < 0.2/N
                    nSuccRPANM=nSuccRPANM+1;
                end
            catch
                Do nothing when error occured!
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Method 6: Non-convex algo.
            display('======================Non-convex algo.========================')
            
            [GainInv,TT]=Manopt_CalPara(Y_noisy_cov_rec);
            
            g_rec=1./GainInv_rec;
            Gain_rec=g_rec/norm(g_rec)*AlphaL2;
            
            %% MUSIC and error
            ImagingFun=MUSIC(TT,sparsity,ell);
            % Local Maximum extraction in MUSIC
            [Supp_rec , supprecind]=LM(ell,ImagingFun,RL,sparsity);                  % recovered support via our algorithm
            % Supp_rec      = rootmusic(TT,sparsity,'corr')/(2*pi);            % recovered support via rootmusic
            % Compute error and adjust trivial ambiguities
            [Supp_rec_c2,Beta_rec_c12,Gain_rec_c,Supp_err,Alpha_reerr,Beta_reerr,Gain_reerr,c0,c1,c2,c] = ComputeReconstructionError(Supp,Supp_rec,Gain,Gain_rec,ell);
            fprintf('Support recovery error = %6.4f RL\n',Supp_err)
            fprintf('Amplitude error mean =%6.4f 0.9 0.5 0.1 quantitle =  %6.4f %6.4f %6.4f \n',mean(Alpha_reerr),quantile(Alpha_reerr,0.9),median(Alpha_reerr),quantile(Alpha_reerr,0.1) );
            fprintf('Phase     error mean =%6.4f 0.9 0.5 0.1 quantitle = %6.4f %6.4f %6.4f \n',mean(Beta_reerr),quantile(Beta_reerr,0.9),median(Beta_reerr),quantile(Beta_reerr,0.1) );
            fprintf('Cali      error mean =%6.4f 0.9 0.5 0.1 quantitle = %6.4f %6.4f %6.4f \n',mean(Gain_reerr),quantile(Gain_reerr,0.9),median(Gain_reerr),quantile(Gain_reerr,0.1) );
            fprintf('cali constant        = %6.2f + i %6.2f  magnitude =%6.2f\n',real(c),imag(c),abs(c))

            % tracking # of success
            if Supp_err < 0.2/N
                nSuccNON=nSuccNON+1;
            end

        end
        resBuf = [resBuf;N,sparsity,RL,nlevel,L,nSuccAlge,nSuccOptWF,nSuccAltFW,nSuccPANM,nSuccRPANM,nSuccNON];
        % save data to xls file
        countXLS = countXLS+1;
        xlRange = sprintf('A%d',countXLS+1);
        xlsInsert = resBuf;
        xlswrite(filename, xlsInsert, 'Sheet1', xlRange);
    end
end




