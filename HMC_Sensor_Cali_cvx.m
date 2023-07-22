function [rh,err]=HMC_Sensor_Cali_cvx(oy,oh,N)

Nc=(N+1)/2;
ry=oy;
%% Solve EMaC
cvx_solver sdpt3
cvx_begin sdp quiet
    variable Y(Nc,Nc) hermitian
    variable Z(Nc,Nc) hermitian
    variable rh(N,1) complex 
    rx=diag(rh)*ry;
    Q1=[Y, hankel(rx(1:Nc),rx(Nc:end)); hankel(rx(1:Nc),rx(Nc:end))', Z];
    minimize 0.5*trace(Y)+0.5*trace(Z)
    subject to
        Q1 >= 0;
        rh(1:Nc) == oh(1:Nc);
cvx_end
% estimate the frequencies using esprit
err = norm(rh-oh)/norm(oh)

save('hmc.mat');
end

