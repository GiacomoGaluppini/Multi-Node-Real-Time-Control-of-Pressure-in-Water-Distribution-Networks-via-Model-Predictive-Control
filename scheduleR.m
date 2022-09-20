function [R] = scheduleR(c,l,Rlow,Rhigh,a)
% [R] = scheduleR(c,l,Rlow,Rhigh,a)
% Scheduling law for R in MPC 
% R = (Rhigh-Rlow)*(-1./(1 + exp(-(a-c)/l))+1)+Rlow;)
% Inputs:
% c: center of scheduling function
% l: width of scheduling function
% Rlow: minimum R value
% Rhigh: maximum R value
% a: current valve closure value
% Outputs:
% R: R value for MPC
%%
h=Rhigh-Rlow;
R = h*(-1./(1 + exp(-(a-c)/l))+1)+Rlow;
end

