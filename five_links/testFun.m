%% Testing if the functions are functionable

%% generate some dummy q,dq,u

addpath 'obj/'
addpath 'robotGen'
addpath 'grad'
% we will use Hermite-Simpson Collocation method to deal with the dynamic
% equation constraints

% Define param
totT = 1;
dT = 0.02; % set sample frequency 50

time = 0:dT:totT;
gaitLen = 1;

qmax = 90/180*pi;

% init
q = 0.05*qmax*sin(0.4*randn(5,1)+time); %just want to use some random but linear function for init
dq = qmax*sin(0.4*randn(5,1)+time);
u = 4*sin(0.4*randn(5,1)+time);
x = [q;dq;u];


%% test betaFun
% betaFun(x(1,:));


% idxList = 1:1:size(x,1);
% fun = @(i)betaFun(x,i);
% out = arrayfun(fun,idxList,'UniformOutput',false);


%% test f_x

% [out_t,grad_t,grad_t_k] = f_x([x(2,:);x(1,:)],0.05);

%% test end_y functions

ypos = end_y_pos(x(:,2).');
ygrad = end_y_poss_grad(x(:,2).');