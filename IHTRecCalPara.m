function [rH,err]=IHTRecCalPara(Y_cov,r,oh,crange)
%% Iterative Hard-Thresholding Algorithm for Hankel Matrix Completion (HMC)

%ox: original (ground truth) signal x
%K: measurement index set
%N: signal dimension
%epsilon: tolerance
%r: rank
%-------------------------
oH=oh*oh';

%% initialization
n=size(Y_cov,1);
rX=zeros(n,n);
itr=1;
oX=Y_cov;
err(itr)=norm(rX-oX,'fro')/norm(oX,'fro');

%% iterations
Gra=zeros(n,n);
stepsize=0.1;
while ((err(itr)>epsilon) & (itr<n*100)) 
    itr=itr+1;
    %% display error
    %err(itr)=norm(RX-OX,'fro')/norm(OX,'fro');
    %fprintf('Recovery error at iteration %d is: %d (%f)\n',itr,Err(itr), Alp);
    
    %% gradient descent
    rXp=rX;
    Gra=(rX-oX);
    rX=rX-stepsize*Gra;

    %% hard thresholding
    try 
        [U,S]=eig(rX); 
    catch
        return;
    end
    [S,SInd]=sort(diag(S),'descend');    
    rX=U(:,SInd(1:r))*diag(S(1:r).*(S(1:r)>=0))*U(:,SInd(1:r))';
    
    %% devided into two factors
    for ii=1:n
        rz(ii,1)=mean(diag(rX,ii-1));
    end
    Tu=toeplitz(rz);
    rG=rX./Tu;
    %% rank 1 projection
    rGp=rG;
    try
        [U,S]=eig(rG);
    catch
        return;
    end
    [S,SInd]=sort(diag(S),'descend');   
    rg=U(:,SInd(1,1))*sqrt(S(1,1));
    rG=rg*rg';
    %%
    err(itr) = norm(rGp - rG,'fro')/ norm(rG,'fro');
end
rH=1./rG;
end