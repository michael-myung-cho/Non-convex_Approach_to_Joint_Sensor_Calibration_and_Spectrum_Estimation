function [ gain,Tu ] = AltRecCalPara( Y_cov,r, oh )
%
% n: signal dimension
% r: rank
% Y_cov: covariance matrix
%
% -----------------------
oH=oh*oh';

%% initialization
n=size(Y_cov,1);
Tu=zeros(n,n);
rH=ones(n,n);
itr=1;
err(itr)=norm(Y_cov.*rH-Tu,'fro');
rz=zeros(n,1);
%% iterations
epsilon=10^-3;
stepsize=0.9;
stepsize2=0;
while ((err(itr)>epsilon) & (itr<n*100)) 
    %% display error
    fprintf('Recovery error at iteration %d is: %d \n',itr,err(itr));
    itr=itr+1; 
    %% Projection onto Toeplitz
    rHp=rH;
    Tup=Tu;
    Tu=stepsize.*Tup+(1-stepsize).*(Y_cov.*rH);
    rzp=rz;
    for ii=1:n
        rz(ii,1)=mean(diag(Tu,ii-1));
    end
    Tu=toeplitz(rz);
    %% normalize Toeplitz
%     wei=(Y_cov(2,1).*oH(2,1))/(Tu(2,1)+10^-6);
%     Tu=Tu./wei;

    %% rank of Toeplitz
    try 
        [U,S]=eig(Tu); 
    catch
        return;
    end
    [S,SInd]=sort(diag(S),'descend');    
%    Tu=U*diag(S.*(S>0))*U';
    Tu=U(:,SInd(1:r))*diag(S(1:r).*(S(1:r)>=0))*U(:,SInd(1:r))';
       
    %% recovery H
    rHp=rH;
    rH=stepsize2.*rHp+(1-stepsize2).*(Tu./Y_cov);

    
    %% rank 1 projection
    try
        [U,S]=eig(rH);
    catch
        return;
    end
    [S,SInd]=sort(diag(S),'descend');   
    rh=U(:,SInd(1,1))*sqrt(S(1,1));
    rH=rh*rh';
    rH(2,1)=oH(2,1);
    rH(1,2)=oH(1,2);
    
    err(itr) = norm(Tu-Tup,'fro')/norm(Tu,'fro');

end
gain=1./rh;
save('temp.mat')
end

