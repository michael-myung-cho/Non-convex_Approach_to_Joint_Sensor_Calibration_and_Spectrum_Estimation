function [g_rec,Supp_rec,Outs] = FriedlanderWeiss(Y_noisy_cov_rec,sparsity,ell,Opts)
% function to recover the calibration parameters via the alternating method
% by Friedlander and Weiss
%
% Reference: B. Friedlander and Weiss, Eigenstructure methods for direction
% finding with sensor gain and phase uncertainties
%
% Input:      - Y_noisy_cov_rec:   covariance matrix of y
%             - sparsity:          sparsity
%             - Opts:              MaxIter:   maximal number of iterations
%                                  OptsFun:   objective value at true
%                                  epsilon:   stop if J_{k+1}-J(k)<=
%                                  epsilon
%                                  g0:        initial guess of g
%                                  Jtrue      J at the true solution
% Output:     - g_rec:             recovered calibration parameters
%             - ImagingFun:        recovered imaging function
%

N                 = size(Y_noisy_cov_rec,1);
RL                = 1/N;

[U , D]           = eig(Y_noisy_cov_rec);
D                 = diag(D);
[~ , MinInd]      = min(D);
if MinInd == 1
    NoiseSpace  = U(:,1:N-sparsity);
elseif MinInd == N
    NoiseSpace  = U(:,sparsity+1:end);
else
    NoiseSpace  = U(:,1:N-sparsity);  
end

Grid              = 0 : ell : 1;
GridLen           = length(Grid);

gk                = Opts.g0;                   % starting point of gain parameters
Outs.Jk           = zeros(Opts.MaxIter,1);
Jk                = inf;
w                 = zeros(N,1);
w(1)              = 1;

for iter = 1 : Opts.MaxIter
    % Step 1: MUSIC for spectral estimation based on gain gk
%     RFun              = zeros(GridLen,1);
%     ImagingFun        = zeros(GridLen,1);
%     for l = 1 : GridLen
%         phi               = gk.*exp(2*pi*1i*(0:N-1)'*Grid(l));
%         RFun(l)           = (norm(NoiseSpace'*phi,2))/norm(phi,2);
%         ImagingFun(l)     = 1/RFun(l);
%     end
    %[Supp_k , ~] = LM(ell,ImagingFun,RL,sparsity);
    %% rootmusic
    F_rec        = spdiags(1./gk,0,N,N)*Y_noisy_cov_rec*spdiags(1./conj(gk),0,N,N);
    f_rec        = zeros(N,1);
    for k = 1 : N
        f_rec(k) = mean(diag(F_rec,-(k-1)));
    end    
    Supp_k      = rootmusic(f_rec,sparsity)/(2*pi);            % MCho: recovered support via rootmusic
    
%     % plot imaging function
%     figure; %stem(Supp,Samp,'ro'); hold on;
%     plot(Grid,ImagingFun,'b*-')
%     title(['Imaging function at Iteration '  num2str(iter)])
    % Step 2: recover calibration parameters
    Qk  = zeros(N,N);
    for l =1:sparsity
        phi   = exp(2*pi*1i*(0:N-1)'*Supp_k(l));
        Qk    = Qk + diag(phi')*(NoiseSpace*NoiseSpace')*diag(phi);
    end
    Jkpre             = Jk;
%     %% Method 1 to find calibration parameters: the eigenvector of Qk associated with the minimum eigenvalue
%     [U , D]           = eig(Qk);
%     D                 = diag(D);
%     [Jk , MinIndex]   = min(D);
%     gk                = U(:,MinIndex)*GainL2;
%     Jk                = Jk^2*GainL2^2;
    %% Method 2 to find calibration parameters with the assumption that g_0 = 1
    Qkinvw            = Qk\w;
    gk                = Qkinvw/Qkinvw(1);
    Jk                = gk'*(Qk*gk);
    Outs.Jk(iter)     = Jk;
    if abs(Jk-Jkpre)< Opts.epsilon
        break;
    end
%    fprintf('iter = %6.0f  Jk = %e Opt = %e  \n',iter,Jk,Opts.Jtrue);
end

fprintf('iter = %6.0f  Jrec = %e Jtrue = %e  \n',iter,Outs.Jk(iter),Opts.Jtrue);
Outs.Iter = iter;


% plot imaging function
% figure; %stem(Supp,Samp,'ro'); hold on;
% plot(Grid,ImagingFun,'b*-')
% title('Imaging function by Friedlander and Weiss')

g_rec     = gk;
Supp_rec  = Supp_k;