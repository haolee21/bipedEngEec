%% Calculate the optimized trajectories for 5 link biped model with GRF
clear;
addpath dyn/
addpath grad/
addpath robotGen/
addpath obj/
addpath gaitCon/
%% simulation parameter

param.numJ=5;
param.toe_th =- 1e2;

param.gaitT = 1;
param.sampT = 0.002;

param.gaitLen = 0.8;

time = 0:param.sampT:param.gaitT;

% set torque/angular velocity constraints
max_tau = 30;
max_vel = 30/180*pi;
%% initialize joint pos and torque
qmax = 170/180/pi;
% q = qmax*sin((2*time/param.gaitT+randn(param.jointNum,1))*pi);
% dq = qmax*sin((2*time/param.gaitT+randn(param.jointNum,1))*pi)*10;
q = [pi/2*ones(1,length(time));
     -pi/2*ones(1,length(time));
     zeros(1,length(time));
     -pi*ones(1,length(time));
     pi/2*ones(1,length(time))];
dq = zeros(param.numJ,length(time));


u = zeros(param.numJ,length(q));

ext_tau = zeros(size(time,2),param.numJ);

x0 = [q;dq;u];
prob.x0 = x0;


% dyndq = @(x,tauext)dyn_dq(x,tauext,jointNum,toe_th);
% rowfun(dyndq,table(x,ext_tau))

%% define constraints
%  dynamic+gait length + vertical pos
prob.nonlcon = @(x)five_link_nonlcon(x,param);

% linear constraints
prob.beq = zeros(param.numJ,1);
Aeq = zeros(param.numJ,size(x0,1)*size(x0,2));
Aeq(1:param.numJ,1:param.numJ)=[1,1,1,1,1;
                                0,0,0,0,1;
                                0,0,0,1,0;
                                0,0,1,0,0;
                                0,1,0,0,0;];
Aeq(1:param.numJ,end-param.numJ*3+1:end-param.numJ*2) = [-1,0,0,0,0;
                                                          0,1,0,0,0;
                                                          0,0,1,0,0;
                                                          0,0,0,1,0;
                                                          0,0,0,0,1;];
prob.Aeq = Aeq;

% upper limit and lower limit for each joints
prob.ub = [pi*ones(1,size(x0,2));
           zeros(1,size(x0,2));
           pi/2*ones(1,size(x0,2));
           -pi/2*ones(1,size(x0,2));
           pi*ones(1,size(x0,2));
           max_vel*ones(param.numJ,size(x0,2));
           max_tau*ones(param.numJ,size(x0,2))];
prob.lb = [zeros(1,size(x0,2));
           -pi*ones(1,size(x0,2));
           -pi/2*ones(1,size(x0,2));
           -3/2*pi*ones(1,size(x0,2));
           zeros(1,size(x0,2));
           -max_vel*ones(param.numJ,size(x0,2));
           -max_tau*ones(param.numJ,size(x0,2))];
           



%% define object function
prob.objective = @(x)objFun(x,param);




%% solve

options = optimoptions('fmincon','MaxIter',1e3,'MaxFunEvals',1e3,...
    'Display','iter','GradObj','on','TolCon',1e-3,'SpecifyConstraintGradient',true,...
    'SpecifyObjectiveGradient',true,'MaxIterations',100,'StepTolerance',1e-15);
prob.options = options;
prob.solver = 'fmincon';

sol = fmincon(prob);