function ImagingFun = MUSIC(f_recOrRx,sparsity,ell)
% function to apply MUSIC given the Toeplitz sequence f_rec or the
% covariance matrix Rx
% Input:        - f_rec:     recovered f
%               - Rx:        covariance matrix
%               - sparsity:  number of frequencies
%               - ell:       gives a grid of [0,1) with spacing ell
%               - Supp:      exact support
%               - Samp:      true amplitudes
% Output:       - Supp_rec   recovered support of frequencies


if min(size(f_recOrRx)) == 1   % the first input is the Toeplitz sequqnce
    f_rec = f_recOrRx;
    N     = length(f_rec);
    Rx    = diag(f_rec(1)*ones(N,1));
    for k = 2:N
        Rx    = Rx + diag(f_rec(k)*ones(N-k+1,1),-(k-1)) + diag(conj(f_rec(k))*ones(N-k+1,1),k-1);
    end
else                           % the first input in the square convariance matrix
    Rx    = f_recOrRx;
    N     = size(Rx,1);
end


[U , D]           = eig(Rx);
D                 = diag(D);
if abs(D(1)) >= abs(D(end))
    NoiseSpace    = U(:,sparsity+1:end);
else
    NoiseSpace    = U(:,1:N-sparsity);
end

% [U , D]           = eigs(Rx,N);
% D                 = diag(D);
% NoiseSpace        = U(:,sparsity+1:end);

    
Grid              = 0 : ell : 1;
GridLen           = length(Grid);
RFun              = zeros(GridLen,1);
ImagingFun        = zeros(GridLen,1);
for k = 1 : GridLen
    phi                  = exp(2*pi*1i*(0:N-1)'*Grid(k));
    RFun(k)              = (norm(NoiseSpace'*phi,2))/norm(phi,2);
    ImagingFun(k)        = 1/RFun(k);
end

