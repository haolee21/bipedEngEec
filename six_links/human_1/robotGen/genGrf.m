%% this script is for generating grf
% the equation we are using is Nonlinear spring and Damper expression, 
% from "Stabilization of a Dynamic Walking Gait Simulation (2005)"


%%
function [Fn,Fs]=genGrf(y_pos,y_vel,th)
% all y_pos,y_vel,th are symbolic expressions
k = 2e6;
e =2.2;
cmax = 1500;
dmax = 1e-3;
u_s = 0.8;
u_d = 0.2;




end