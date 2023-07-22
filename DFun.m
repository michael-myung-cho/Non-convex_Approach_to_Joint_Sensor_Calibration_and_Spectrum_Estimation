function D = DFun(GF,Ry)
% function to compute the gradient of the objective function ||GFG'-Ry||_F^2 where 
%                     G = diag(g) with g = gr + 1i * gi
%                     F = diag(f,1:-1:-N)+diag(conj(f(2:end)),2:N)

N     = size(GF,1);
gr    = GF(:,1);    gi = GF(:,2);     fr = GF(:,3);     fi = GF(:,4);

%% Form matrix G and F
g     = gr+1i*gi;   gbar = conj(g);
f     = fr + 1i * fi;
F     = diag(f(1)*ones(N,1));
for k = 2:N
    F = F + diag(f(k)*ones(N-k+1,1),-(k-1)) + diag(conj(f(k))*ones(N-k+1,1),k-1);
end

%% Compute the gradient
DGr   = zeros(N,1);   DGi   = zeros(N,1);  DFr   = zeros(N,1);  DFi   = zeros(N,1);
for m = 1 : N
    % diff with respect to gr and gi
    midseq = (F(m,1:m-1)).'.*gbar(1:m-1);
    DGr(m) = 2*real(midseq'*(g(m)*midseq-(Ry(m,1:m-1)).'));
    DGi(m) = 2*real((1i*midseq)'*(g(m)*midseq-(Ry(m,1:m-1)).'));
    %
    DGr(m) = DGr(m) + 2*real(conj(f(1)*2*gr(m))*(g(m)*gbar(m)*f(1)-Ry(m,m)));
    DGi(m) = DGi(m) + 2*real(conj(f(1)*2*gi(m))*(g(m)*gbar(m)*f(1)-Ry(m,m)));
    %
    midseq = g(m+1:end).*F(m+1:end,m);
    DGr(m) = DGr(m) + 2*real(midseq'*(midseq*gbar(m)-Ry(m+1:end,m)));
    DGi(m) = DGi(m) + 2*real((-1i*midseq)'*(midseq*gbar(m)-Ry(m+1:end,m)));
    % diff with respect to f
    midseq = g(m:end).*gbar(1:N-m+1);
    DFr(m) = 2*real(midseq'*(f(m)*midseq-diag(Ry,-m+1)));
    DFi(m) = 2*real((1i*midseq)'*(f(m)*midseq-diag(Ry,-m+1)));
end
    
D = [DGr DGi DFr DFi];
