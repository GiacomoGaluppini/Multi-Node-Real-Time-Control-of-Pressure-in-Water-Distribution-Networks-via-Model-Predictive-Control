function [X,Y,P,convMesure] = kalmd_upd(sysd,weights,X,U,Y,Ktype)
% [X,Y,P,convMesure] = kalmd_upd(sysd,weights,X,U,Y,Ktype)
% Get state estimation and update Kalman Filter/Predictor
% Inputs:
% sysd: discrete time state stace system object
% weights: structure (with fileds Q,R,P) with Kalman weights and estimate covariance matrix
% X: current state estimate
% U: input measurements
% Y: output measurements
% Ktype: string, 'filter' or 'predictor', defining the type of Kalman estimator
% Outputs:
% X: upated state estimate
% Y: upated output estimate
% P: updated estimate covariance matrix
% convMesure: measure of convergence for P (inf norm)
%%
Pold=weights.P;
switch lower(Ktype)
    case 'filter'
     
        weights.P=sysd.A*weights.P*sysd.A'+weights.Q;
        L=weights.P*sysd.C'*inv(sysd.C*weights.P*sysd.C'+weights.R);
        weights.P=weights.P-L*sysd.C*weights.P;
        
        X=sysd.A*X+sysd.B*U+L*(Y-sysd.C*(sysd.A*X+sysd.B*U));
        Y=sysd.C*X;
        P=weights.P;
        
    case 'predictor'
        
        weights.P=sysd.A*weights.P*sysd.A'+weights.Q+...
            -sysd.A*weights.P*sysd.C'*inv(sysd.C*weights.P*sysd.C'+weights.R)*sysd.C*weights.P*sysd.A';
        L=sysd.A*weights.P*sysd.C'*inv(sysd.C*weights.P*sysd.C'+weights.R);
        
        
        X=sysd.A*X+sysd.B*U+L*(Y-sysd.C*X);
        Y=sysd.C*X;
        P=weights.P;
        
end

convMesure=norm(weights.P-Pold,Inf);

end

