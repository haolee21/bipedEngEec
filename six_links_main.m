%% This script will try to calculate the optimized trajectory of given model and boundary conditions

addpath 'obj/'
% we will use Hermite-Simpson Collocation method to deal with the dynamic
% equation constraints

%% Define param
totT = 1;
dT = 0.02; % set sample frequency 50

time = 0:dT:totT;
gaitLen = 1;

qmax = 90/180*pi;

%% init
q = 0.05*qmax*sin(0.4*randn(5,1)+time); %just want to use some random but linear function for init
dq = qmax*sin(0.4*randn(5,1)+time);
u = 4*sin(0.4*randn(5,1)+time);
x = [q;dq;u];



%% Define object function

options = optimoptions('fmincon','SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true);


obj(x) = objFun(x,5);
