function [Gain_rec , ImagingFun, f_rec] = AlgebraicMethod(Y_noisy_cov_rec,sparsity,ell,AlphaL2)
% Algebraic method for spectral estimation and sensor calibration 
% Input:   - Y_noisy_cov_rec: covariance matrix of y
%          - sparsity:        number of frequencies / objects to recover
% Output:  - Gain_rec:        recovered calibration parameters
%          - Supp_rec:        recovered support of the frequencies

N  = size(Y_noisy_cov_rec,1);
RL = 1/N;

%% Step 1: Estimate calibration amplitudes
Y_noisy_cov_diag  = diag(Y_noisy_cov_rec);
Alpha_rec         = sqrt(Y_noisy_cov_diag); 
Alpha_rec         = Alpha_rec/norm(Alpha_rec)*AlphaL2;

%% Step 2: Estimate the calibration phases
RHS          = diag(Y_noisy_cov_rec,-1);
RHS          = angle(RHS(2:end)./RHS(1:end-1));
RHS          =[RHS ; 0 ; 0];
Psi1         = spdiags([ones(N-2,1) -2*ones(N-2,1) ones(N-2,1)], 0:2, N-2, N);
Psi          = [Psi1 ; sparse(2,N)];
Psi(N-1,N-1) = 1;
Psi(N,N)     = 1;
Beta_rec     = mod(Psi\RHS,2*pi);
Gain_rec     = Alpha_rec.*exp(1i*Beta_rec);

LSres        = norm(Psi*Beta_rec-RHS)/norm(RHS);
fprintf('LS relative res = %6.4f\n',LSres)

%% Step 3: Apply MUSIC for spectral estimation
F_rec        = spdiags(1./Gain_rec,0,N,N)*Y_noisy_cov_rec*spdiags(1./conj(Gain_rec),0,N,N);
f_rec        = zeros(N,1);
for k = 1 : N
    f_rec(k) = mean(diag(F_rec,-(k-1)));
end
[U , D]      = eig(F_rec);
D            = diag(D);
if abs(D(1)) >= abs(D(end))
    NoiseSpace   = U(:,sparsity+1:end);
else
    NoiseSpace  = U(:,1:N-sparsity);
end
Grid         = 0 : ell : 1;
GridLen      = length(Grid);
RFun         = zeros(GridLen,1);
ImagingFun   = zeros(GridLen,1);
for k = 1 : GridLen
    phi                  = exp(2*pi*1i*(0:N-1)'*Grid(k));
    RFun(k)              = (norm(NoiseSpace'*phi,2))/norm(phi,2);
    ImagingFun(k)        = 1/RFun(k);
end

