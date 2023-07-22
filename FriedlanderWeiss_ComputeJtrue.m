function Jtrue = FriedlanderWeiss_ComputeJtrue(Y_noisy_cov_rec,Gain,Supp)

N           = size(Y_noisy_cov_rec,1);
sparsity    = length(Supp);
GainL2      = norm(Gain);


% noise space
[U , D]      = eig(Y_noisy_cov_rec);
D            = diag(D);
if abs(D(1)) >= abs(D(end))
    NoiseSpace    = U(:,sparsity+1:end);
else
    NoiseSpace    = U(:,1:N-sparsity);
end

% % Compute Jtrue as Eq. (8) in the paper by Friedlander and Weiss
% Jtrue = 0;
% for l =1:sparsity
%     phi      = Gain.*exp(2*pi*1i*(0:N-1)'*Supp(l));
%     Jtrue    = Jtrue + norm(NoiseSpace'*phi)^2;
% end
% fprintf('Jtrue =%e through Eq. (8)\n',Jtrue)
% 

Qtrue        = zeros(N,N);
for l =1:sparsity
    phi      = exp(2*pi*1i*(0:N-1)'*Supp(l));
    Qtrue    = Qtrue + diag(phi')*(NoiseSpace*NoiseSpace')*diag(phi);
end
% [U , D]           = eig(Qtrue);
% D                 = diag(D);
% [Jtrue , MinInd]  = min(D);
% Jtrue             = Jtrue^2*GainL2^2;
% fprintf('Jtrue =%e of mimimum eigenvalue of Q \n',Jtrue)


Jtrue             = Gain'*Qtrue*Gain;
fprintf('Jtrue =%e throught g^H*Q*g \n',Jtrue)

