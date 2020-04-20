%% Calculate the optimized trajectories for 5 link biped model with GRF
clear;
addpath dyn/
addpath grad/
addpath robotGen/
addpath obj/
addpath gaitCon/
addpath plotRobot/
%% simulation parameter

param.numJ=5;
param.toe_th = 1e-2;

param.gaitT = 0.8;
param.sampT = 0.005;

param.gaitLen = 1.8;

time = 0:param.sampT:param.gaitT;

% set torque/angular velocity constraints
max_tau = 5;
max_vel = 120/180*pi/param.sampT;
%% initialize joint pos and torque
qmax = 170/180/pi;
% q = qmax*sin((2*time/param.gaitT+randn(param.jointNum,1))*pi);
% dq = qmax*sin((2*time/param.gaitT+randn(param.jointNum,1))*pi)*10;

qStart=[70/180*pi,-10/180*pi,45/180*pi,-pi,20/180*pi];
qEnd = [pi+qStart(1)+qStart(2)+qStart(3)+qStart(4)+qStart(5),...
        -qStart(5),...
        qStart(4)+pi,...
        -pi-qStart(3),...
        -qStart(2)];

q = [linspace(qStart(1),qEnd(1),length(time));
     linspace(qStart(2),qEnd(2),length(time));
     linspace(qStart(3),qEnd(3),length(time));
     linspace(qStart(4),qEnd(4),length(time));
     linspace(qStart(5),qEnd(5),length(time))];
dq = [0,(q(1,2:end)-q(1,1:end-1))/param.sampT;
      0,(q(2,2:end)-q(2,1:end-1))/param.sampT;
      0,(q(3,2:end)-q(3,1:end-1))/param.sampT;
      0,(q(4,2:end)-q(4,1:end-1))/param.sampT;
      0,(q(5,2:end)-q(5,1:end-1))/param.sampT];


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
                                                          0,0,-1,0,0;
                                                          0,0,0,1,0;
                                                          0,0,0,0,1;];
prob.Aeq = Aeq;
prob.beq = [-pi;0;-pi;-pi;0];
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

options = optimoptions('fmincon','MaxIter',200,...
    'Display','iter','GradObj','on','TolCon',1e-3,'SpecifyConstraintGradient',true,...
    'SpecifyObjectiveGradient',true,'StepTolerance',1e-15,'UseParallel',true);
prob.options = options;
prob.solver = 'fmincon';

sol = fmincon(prob);