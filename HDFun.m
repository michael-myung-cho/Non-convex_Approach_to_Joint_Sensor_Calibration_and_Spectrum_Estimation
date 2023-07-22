function H = HDFun(GF,Ry)
% function to compute the Hessian of the objective function at GF through
% Wirtinger derivative
% ||GfG'-Ry||_F^2 where G=diag(g) with g=gr+i*gi
%                       G = diag(g) with g = gr + 1i * gi
%                       F = diag(f,1:-1:-N)+diag(conj(f(2:end)),2:N)  
 
N   =size(GF,1);                                  % number of sensors 
 
%% Form matrix G and F                   
g   = GF(:,1);      gbar = conj(g);
f   = GF(:,2);      fbar = conj(f);

%% Compute the Hessian Matrix 
H=zeros(4*N,4*N); 


%% construct the upper diagonal entries 

for n=1:N
    H(2*n-1,2*n)=gbar(n)^2*f(1)^2;                                         % partial g_n g_n
    for m=n+1:N
        H(2*n-1,2*m-1)=fbar(m-n+1)*(g(m)*gbar(n)*f(m-n+1)-Ry(m,n));        % partial g_n gbar_m
        H(2*n-1,2*m)=gbar(m)*abs(f(m-n+1))^2*gbar(n);                      % partial g_n g_m
        %% Sui
        H(2*n,2*m-1)=H(2*n-1,2*m)';   
        H(2*n,2*m)=H(2*n-1,2*m-1)';
    end
end

for n=1:N
    for m=1:N
        if 1<=m && m<=n
            H(2*n-1,2*N+2*m-1)=H(2*n-1,2*N+2*m-1)+abs(g(n-m+1))^2*gbar(n)*f(m);                     % partial g_n barf_m
            H(2*n-1,2*N+2*m)=H(2*n-1,2*N+2*m)+gbar(n-m+1)*(g(n-m+1)*gbar(n)*fbar(m)-Ry(n,n-m+1)');  % partial g_n f_m
        end
        if 2<=m && m<=N-n+1
            H(2*n-1,2*N+2*m-1)=H(2*n-1,2*N+2*m-1)+gbar(n+m-1)*(g(n+m-1)*gbar(n)*f(m)-Ry(n+m-1,n));  % partial g_n barf_m
            H(2*n-1,2*N+2*m)=H(2*n-1,2*N+2*m)+abs(gbar(n+m-1))^2*fbar(m)*gbar(n);                   % partial g_n f_m
        end
    end
end




for n=1:N
    for m=1:N
        %% Sui
        H(2*n, 2*N+2*m-1) = H(2*n-1,2*N+2*m)';                               % partial barg_n barf_m
        H(2*n,2*N+2*m)    = H(2*n-1,2*N+2*m-1)';                                % partial barg_n f_m
    end
end


H=H+H';


%% construct the diagonal entries        
for n=1:N
%     H(2*n-1,2*n-1)  = ones(n,1)'*abs(g(1:n).*f(n:-1:1)).^2+ ones(N-n,1)'*abs(g(n+1:N).*f(2:N-n+1)).^2+f(1)*(abs(g(n))^2*fbar(1)-Ry(n,n));
%     H(2*n,2*n)      = ones(n,1)'*abs(g(1:n).*f(n:-1:1)).^2+ ones(N-n,1)'*abs(g(n+1:N).*f(2:N-n+1)).^2+f(1)*(abs(g(n))^2*fbar(1)-Ry(n,n));
    %% Sui
    H(2*n-1,2*n-1)  = sum(abs(g(1:n).*f(n:-1:1)).^2) + sum(abs(g(n+1:N).*f(2:N-n+1)).^2) + f(1)*(abs(g(n))^2*fbar(1)-Ry(n,n));
    %% Wenjing
    H(2*n-1,2*n-1)  = sum(abs(g(1:n).*f(n:-1:1)).^2) + sum(abs(g(n+1:N).*f(2:N-n+1)).^2) + f(1)*(abs(g(n))^2*fbar(1)-Ry(n,n)) + f(1).^2*g(n)*g(n)';
    H(2*n,2*n)      = H(2*n-1,2*n-1);
end

% H(2*N-1,2*N-1)=ones(N,1)'*abs(g(1:N).*f(N:-1:1)).^2;
% H(2*N,2*N)=ones(N,1)'*abs(g(1:N).*f(N:-1:1)).^2;


for n=1:N
%     H(2*N+2*n-1,2*N+2*n-1)  = ones(N-n+1,1)'*abs(g(1:N-n+1).*g(n:N)).^2;
%     H(2*N+2*n,2*N+2*n)      = ones(N-n+1,1)'*abs(g(1:N-n+1).*g(n:N)).^2;
    %% same
    H(2*N+2*n-1,2*N+2*n-1)  = sum(abs(g(1:N-n+1).*g(n:N)).^2);
    H(2*N+2*n,2*N+2*n)      = H(2*N+2*n-1,2*N+2*n-1);
end


end