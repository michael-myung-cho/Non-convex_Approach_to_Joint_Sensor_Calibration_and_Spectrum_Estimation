function [ B ] = toep( A )
    [n,m]=size(A);
    ind=1;
    for ii=0:-1:-n+1
        vA(ind)=sum(diag(A,ii))/size(diag(A,ii),1);
        ind=ind+1;
    end
    for ii=0:n-1
        vA(ind)=sum(diag(A,ii))/size(diag(A,ii),1);
        ind=ind+1;
    end
    B=toeplitz(vA(1:n),vA(n+1:end));
end

