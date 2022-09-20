clear all
close all
clc

s=tf('s');
addpath(genpath('C:\Program Files\Mosek\9.2\toolbox\R2015a'))
%% Working Point and Bounds

%Indexes of controlled nodes
WP.criticalNodes=[1 24];

%Upper & lower bound for valve closure
WP.Bounds.amax=0.9;
WP.Bounds.amin=0;

%Upper & lower bound for pressure
WP.Bounds.hmin=[20 20]';
WP.Bounds.hmax=[40 40]';

%Upper & lower bound for valve operations rate
WP.Bounds.adotabsmax=1/300;

%Working point specification
WP.a=0.61;
WP.h=[25.3 25.7]';
WP.Q=0.06;

%Valve specifications for valve closure to local head loss conversion
ValveInfo.Cvmax=0.093;
ValveInfo.diam=0.25;

%Valve closure to local head loss conversion for WP
WP.Ebar=a2Cv2E_CastelfrancoEmilia(WP.a,ValveInfo.Cvmax,ValveInfo.diam);

%Valve closure bounds to local head loss bounds conversion for WP
WP.Bounds.Emax=a2Cv2E_CastelfrancoEmilia(WP.Bounds.amax,ValveInfo.Cvmax,ValveInfo.diam);
WP.Bounds.Emin=a2Cv2E_CastelfrancoEmilia(WP.Bounds.amin,ValveInfo.Cvmax,ValveInfo.diam);

%WP for local model (variation signals)
WP.Bounds.dEmin=WP.Bounds.Emin-WP.Ebar;
WP.Bounds.dEmax=WP.Bounds.Emax-WP.Ebar;

WP.Bounds.dhmin=WP.Bounds.hmin-WP.h;
WP.Bounds.dhmax=WP.Bounds.hmax-WP.h;

%Pressure setpoint bounds (variation signals)
WP.Bounds.minSetpoint=0.8*WP.Bounds.dhmin;
WP.Bounds.maxSetpoint=WP.Bounds.dhmax;
%% System Definition

%Control sampling time
MPC.Ts=3;

%Load a model (parametric model object from Matlab Identification Toolbox)
load('HOTFmodel_WDN_SIMO_2.mat')

%Discretization
Gz=c2d(tf(Gs),MPC.Ts);

%Realizaion
S_mpc=ss(Gz);

%% Order reduction (if required)

figure(1)
hsvd(S_mpc)

%Selected order for reduced order model
orders=[11];
newsys=balred(S_mpc,orders);

figure(2)
clf
impulse(S_mpc,newsys)
figure(3)
clf
bode(S_mpc,newsys)

keyboard

S_mpc=newsys;

%% System Definition (continued)

%Absorb delays
S_mpc=absorbDelay(S_mpc);

WDN=S_mpc;

nx=size(WDN.A,1);
ny=size(WDN.C,1);
nu=size(WDN.B,2);

%% System Definition (extension with integral action)

%add integral of control action for MPC
[As,Bs,Cs,Ds]=ssdata(S_mpc);

Abig=[As Bs; zeros(nu,nx) eye(nu)];
Bbig=[zeros(nx,1);MPC.Ts*ones(nu,1)];
Cbig=[Cs Ds];
Dbig=zeros(ny,nu);

MPC.sys=ss(Abig,Bbig,Cbig,Dbig,MPC.Ts);

%verify controllability
MR=ctrb(Abig,Bbig);
MPC.isReach=rank(MR)==min(size(MR));

nx=size(MPC.sys.A,1);
ny=size(MPC.sys.C,1);
nu=size(MPC.sys.B,2);
%% Filters for Q_v and Alpha (Tf=filter time constant)

Tfq=180;
FiltroQ=c2d(1/(1+s*Tfq),MPC.Ts);
FiltroQ=ss(FiltroQ);

Tfa=180;
FiltroA=c2d(1/(1+s*Tfa),MPC.Ts);
FiltroA=ss(FiltroA);

%% Kalman 

Ak=blkdiag(WDN.A,eye(ny));
Bk=[WDN.B; zeros(ny,1)];
Ck=[WDN.C eye(ny)];

Kalman.sys=ss(Ak,Bk,Ck,[],MPC.Ts);

Kalman.Q=blkdiag(eye(nx-nu),0.1*eye(ny));
Kalman.R=eye(ny);
Kalman.P = dare(Ak',Ck',Kalman.Q,Kalman.R,[],[]);


%Verify Convergence of P
convMesure=Inf;
while convMesure>1e-5
    [~,~,Kalman.P,convMesure] = kalmd_upd(Kalman.sys,Kalman,0,0,0,'predictor');
    convMesure
end

%Verify Observability
MO=obsv(Ak,Ck);
Kalman.isObsv=rank(MO)==min(size(MO));

%% Steady State Auxiliary Target Calculator

Qss=eye(ny);
Rss=1e-3*eye(nu);
H=0.5*ones(nu,ny);

Bt=sdpvar(nx-1,1,'full');
xt=sdpvar(nx-1,1,'full');
ut=sdpvar(nu,1,'full');
dt=sdpvar(ny,1,'full');

%Soft constraints slack variable
sCtol=sdpvar(1,1,'full');

constr=[ xt==WDN.A*xt+Bt*ut];
constr = [constr, WP.Bounds.dEmin <= ut<= WP.Bounds.dEmax];

constr = [constr,0<=sCtol];
constr = [constr,WP.Bounds.minSetpoint-sCtol<=WDN.C*xt+dt];
constr = [constr,WDN.C*xt+dt<=WP.Bounds.maxSetpoint+sCtol];

cst=(WDN.C*xt+dt)'*Qss*(WDN.C*xt+dt)+ut'*Rss*ut;
cst=cst+(H*(WDN.C*xt+dt))'*1e5*(H*(WDN.C*xt+dt));%For Offset Free Control
cst=cst+1e10*sCtol;

ops = sdpsettings('solver','mosek','verbose',0);

%Save optimiser object
SScalculator.Trimmer = optimizer(constr, cst,ops,{Bt,dt},{xt,ut});
SScalculator.sys=WDN;

%% Model Predictive Controller

MPC.Qy=blkdiag(1,1);

%Scheduling law for R
centro=0.3;
larghezza=0.035;
Rlow=10;
Rhigh=100;
MPC.R=@(x) scheduleR(centro,larghezza,Rlow,Rhigh,x);

MPC.N=60;
%Not used
% MPC.alpha=10;

U = sdpvar(repmat(nu,1,MPC.N),repmat(1,1,MPC.N));
X = sdpvar(repmat(nx,1,MPC.N+1),repmat(1,1,MPC.N+1));

%Optimisation parameters (steady state values)
xss=sdpvar(nx,1,'full');
yss=sdpvar(ny,1,'full');

%Optimisation parameters (Dynamic matrix A, Terminal weight Py, Input weight R)
A=sdpvar(size(MPC.sys.A,1),size(MPC.sys.A,2),'full');
Py=sdpvar(nx,nx,'symmetric');
R=sdpvar(nu,nu,'symmetric');

%Optimisation parameters (disturance value - constant along horizon)
d=sdpvar(ny,1);

%Optimisation parameters (bounds on local head loss rate)
Edotmaxminus=sdpvar(1,1);
Edotmaxplus=sdpvar(1,1);

%Soft constraints slack variable
softCtol=sdpvar(1,1);

objective = softCtol*1e5;
constraints = [0<=softCtol];
for k = 1:MPC.N
    objective = objective + (MPC.sys.C*X{k}+d-yss)'*MPC.Qy*(MPC.sys.C*X{k}+d-yss);
    objective = objective + U{k}'*R*U{k};
 constraints = [constraints, X{k+1} == A*X{k} + MPC.sys.B*U{k}];
 constraints = [constraints, WP.Bounds.dEmin <= X{k}(end-nu+1:end)<= WP.Bounds.dEmax];
 constraints = [constraints, Edotmaxminus<=U{k}<=Edotmaxplus];
 constraints = [constraints, WP.Bounds.dhmin-softCtol<=(MPC.sys.C*X{k}+d)];
 constraints = [constraints, (MPC.sys.C*X{k}+d)<=WP.Bounds.dhmax+softCtol];
 
end
constraints = [constraints, WP.Bounds.dhmin-softCtol<=(MPC.sys.C*X{MPC.N+1}+d)];
constraints = [constraints, (MPC.sys.C*X{MPC.N+1}+d)<=WP.Bounds.dhmax+softCtol];
constraints = [constraints,WP.Bounds.dEmin <= X{MPC.N+1}(end-nu+1:end)<= WP.Bounds.dEmax];

objective = objective + (X{MPC.N+1}-xss)'*Py*(X{MPC.N+1}-xss);

%Save optimiser object
parameters_in = {X{1},A,R,Py,d,Edotmaxminus,Edotmaxplus,xss,yss};
ops = sdpsettings('solver','mosek');
MPC.Controller = optimizer(constraints, objective,ops,parameters_in,[U{:}]);

%% Save objects

keyboard

save('ControlloreMPC-OffsetFree','WP','MPC','Kalman','SScalculator','FiltroQ','FiltroA','ValveInfo')
save('WDNlinearsimulator-OffsetFree','WDN')

