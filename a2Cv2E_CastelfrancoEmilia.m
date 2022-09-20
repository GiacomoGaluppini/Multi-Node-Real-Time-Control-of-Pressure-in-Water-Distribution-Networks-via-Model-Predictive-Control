function [E] = a2Cv2E_CastelfrancoEmilia(a,Cvmax,D)
%[E] = a2Cv2E_CastelfrancoEmilia(a,Cvmax,D)
% Compute valve curve (valve closure to local head loss coefficient) for
% Castelfranco Emilia WDN
% Inputs:
% a: valve closure
% Cvmax: valve maximum Cv
% D: valve diameter
% Outputs:
% E: local head loss coefficient
%%
xm=(1-a);
Cv=zeros(size(xm));
for i=1:length(xm)
    if xm(i)<=0.2
        Cv(i)=0.32.*xm(i).*Cvmax;
    else
        Cv(i)=(-3.25.*xm(i).^3+6.*xm(i).^2-2.*xm(i)+0.25).*Cvmax;
    end
end
E=2*9.81.*(3.14.*D.^2./4).^2./Cv.^2;
end

