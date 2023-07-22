function [Supp_rec_c2,Beta_rec_c12,Gain_rec_c,Supp_err,Alpha_reerr,Beta_reerr,Gain_reerr,c0,c1,c2,c] = ComputeReconstructionError4RL(Supp,Supp_rec,Gain,Gain_rec,ell,RL)
% function to adjust trivial ambiguities and compute reconstruction error
% Input:     - Supp:         exact support
%            - Supp_rec:     recovered support
%            - Gain:         exact calibration parameters
%            - Gain_rec:     recovered calibration parameters
%            - ell:          gives a grid on [0,1) withs spacing ell
% Output:    - Supp_rec_c2:  recovered support shifted
%            - Beta_rec_c12: recovered beta with c1,c2 adjusted
%            - Supp_err:     support error mod a translation
%            - Alpha_reerr:  relative error of Alpha
%            - Beta_reerr:   relative error of Beta


Alpha        = abs(Gain);
Beta         = angle(Gain);
Alpha_rec    = abs(Gain_rec);
Beta_rec     = angle(Gain_rec);
Grid         = 0 : ell : 1;
GridLen      = length(Grid);
N            = length(Gain);
%RL           = 1/N;

%% Find the Hausdorff distance between Supp and Supp_rec % Changed for Resolution limit
Supp_rec=sort((Supp_rec < 0)+Supp_rec);

shift1=Supp_rec(1)-Supp(1);
dist2=Supp_rec(2)+shift1-Supp(2);
shift2=Supp_rec(2)-Supp(2);
dist1=Supp_rec(1)+shift2-Supp(1);

Supp_err=min(dist1,dist2)

% DistShiftGrid      = zeros(GridLen,1);
% for k = 1:GridLen
%     DistShiftGrid(k) = B_dist(Supp,mod(Supp_rec+Grid(k),1)); %/RL; % (MCho) why dividing B_Dist with RL? for unit
% end
% [Supp_err , Supp_shift] = min(DistShiftGrid);
% Supp_shift              = Grid(Supp_shift);
% Supp_rec_c2             = (mod(Supp_rec+Supp_shift,1))';

%% Find the calibration error
c2                      = Supp_shift*2*pi;
Beta_rec_c2             = wrapTo2Pi(Beta_rec-(0:N-1)'*c2);
Gain_rec_c2             = Alpha_rec.*exp(1i*Beta_rec_c2);

%% Amplitude scaling and phase shift  
c0                      = 1; % since we assume ||\alpha|| is given
c1                      = median(wrapToPi(Beta-Beta_rec_c2));
Beta_rec_c12            = wrapToPi(Beta_rec_c2+c1);
Beta_reerr              = abs(wrapToPi(Beta-Beta_rec_c12))/pi;
Alpha_reerr             = abs(Alpha_rec-Alpha)./Alpha;

%% Calibration scaling and phase shift: min ||c*Gain_rec_c2-Gain||^2 
cr                      = real(Gain'*Gain_rec_c2)/(Gain_rec_c2'*Gain_rec_c2);
ci                      = -imag(Gain'*Gain_rec_c2)/(Gain_rec_c2'*Gain_rec_c2);
c                       = cr + 1i * ci;
Gain_rec_c              = c*Gain_rec_c2;
Gain_reerr              = abs(Gain-Gain_rec_c)./abs(Gain);



