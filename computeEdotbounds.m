function [Edotmaxplus,Edotmaxminus] = computeEdotbounds(a,adotabsmax,Cvmax,diam,Ts)
% [Edotmaxplus,Edotmaxminus] = computeEdotbounds(a,adotabsmax,Cvmax,diam,Ts)
% Convert bounds on valve operations rate to local head loss rate
% Inputs:
% a: current valve closure
% adotabsmax: maximum valve operation rate
% Cvmax: valve maximum Cv
% diam: valve diameter
% Ts: control sampling time
% Outputs:
% Edotmaxplus: maximum local head loss rate
% Edotmaxminus: minimum local head loss rate
%%
aplus=a+adotabsmax*Ts;
aminus=a-adotabsmax*Ts;
Edotmaxplus=(a2Cv2E(aplus,Cvmax,diam)-a2Cv2E(a,Cvmax,diam))/(Ts);
Edotmaxminus=(a2Cv2E(aminus,Cvmax,diam)-a2Cv2E(a,Cvmax,diam))/(Ts);
end

