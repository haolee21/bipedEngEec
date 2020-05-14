%% Calculate the optimized trajectories for 5 link biped model with GRF
clear;
addpath dyn/
addpath robotGen/grad/
addpath robotGen/
addpath robotGen/posCons/
addpath robotGen/dyn/
addpath robotGen/obj/
addpath obj/
addpath gaitCon/
addpath plotRobot/
addpath robotGen/grf/
%% simulation parameter

param.numJ=6;
param.toe_th =-8e-2;
param.head_h = 1.4 ; %the head should be at least 1.6m
param.fri_coeff=50;
param.gaitT = 0.3;
param.sampT = 0.005;
param.init_y = -7e-2; %initial feet height
param.heel_h = 7e-2;
param.gaitLen = 1.7;
param.hipLen=0.3;
param.gndclear = 0;
param.jointW = [0.1,1,1,1,1,0.1];
time = 0:param.sampT:param.gaitT;

% set torque/angular velocity constraints
param.max_tau = 10;
param.max_vel = 360/180*pi;

param.max_front_ank_tau = param.max_tau*0.01;
%% initialize joint pos and torque
qmax = 170/180/pi;
% q = qmax*sin((2*time/param.gaitT+randn(param.jointNum,1))*pi);
% dq = qmax*sin((2*time/param.gaitT+randn(param.jointNum,1))*pi)*10;
q1 = 75;
q2 = -5;
q3 =45;
q4 = -200;
q5 =2;
q6 = -180-q1-q2-q3-q4-q5;
qStart=[q1/180*pi,q2/180*pi,q3/180*pi,q4/180*pi,q5/180*pi,q6/180*pi];
qEnd = [pi+qStart(1)+qStart(2)+qStart(3)+qStart(4)+qStart(5),...
        -qStart(5),...
        -qStart(4)-pi,...
        -pi-qStart(3),...
        -qStart(2),...
        -qStart(1)];

q = [linspace(qStart(1),qEnd(1),length(time));
     linspace(qStart(2),qEnd(2),length(time));
     linspace(qStart(3),qEnd(3),length(time));
     linspace(qStart(4),qEnd(4),length(time));
     linspace(qStart(5),qEnd(5),length(time));
     linspace(qStart(6),qEnd(6),length(time))];
dq = [0,(q(1,2:end)-q(1,1:end-1))/param.sampT;
      0,(q(2,2:end)-q(2,1:end-1))/param.sampT;
      0,(q(3,2:end)-q(3,1:end-1))/param.sampT;
      0,(q(4,2:end)-q(4,1:end-1))/param.sampT;
      0,(q(5,2:end)-q(5,1:end-1))/param.sampT;
      0,(q(6,2:end)-q(6,1:end-1))/param.sampT];


u = zeros(param.numJ,length(q));

ext_tau = zeros(size(time,2),param.numJ);

x0 = [q;dq;u];



% dyndq = @(x,tauext)dyn_dq(x,tauext,jointNum,toe_th);
% rowfun(dyndq,table(x,ext_tau))

%% define constraints
%  dynamic+gait length + vertical pos
prob.nonlcon = @(x)five_link_nonlcon(x,param);

% linear constraints

Aeq = zeros(param.numJ,size(x0,1)*size(x0,2)); %start-end joint condition 
Aeq(1:param.numJ,1:param.numJ)=[1,1,1,1,1,0;   %start frame 
                                0,0,0,0,1,0;
                                0,0,0,1,0,0;
                                0,0,1,0,0,0;
                                0,1,0,0,0,0;
                                1,0,0,0,0,0;];
Aeq(1:param.numJ,end-param.numJ*3+1:end-param.numJ*2) = [-1,0,0,0,0,0; % endframe
                                                          0,1,0,0,0,0;
                                                          0,0,1,0,0,0;
                                                          0,0,0,1,0,0;
                                                          0,0,0,0,1,0;
                                                          0,0,0,0,0,1];
prob.Aeq = Aeq;
prob.beq = [-pi;0;-pi;-pi;0;0];

% back never bend backward
% -1*q3 <-90 deg, 
Asamp = zeros(1,3*param.numJ); % create A for single frame
Asamp(1:3) = [-1,-1,-1];
Acell = repmat({Asamp},1,param.gaitT/param.sampT+1);
prob.Aineq = blkdiag(Acell{:});
prob.bineq = -pi/2*ones(param.gaitT/param.sampT+1,1); 

% upper limit and lower limit for each joints
prob.ub = [180/180*pi*ones(1,size(x0,2));
           zeros(1,size(x0,2))/180*pi;
           90/180*pi*ones(1,size(x0,2));
           -90/180*pi*ones(1,size(x0,2));
           180/180*pi*ones(1,size(x0,2));
           -10/180*pi*ones(1,size(x0,2));
           param.max_vel*ones(param.numJ,size(x0,2));
           param.max_front_ank_tau*ones(1,size(x0,2));
           param.max_tau*ones(param.numJ-1,size(x0,2))];
prob.lb = [ones(1,size(x0,2))/180*pi;
           -180/180*pi*ones(1,size(x0,2));
           -90/180*pi*ones(1,size(x0,2));
           -270/180*pi*ones(1,size(x0,2));
           zeros(1,size(x0,2))/180*pi;
           -145/180*pi*ones(1,size(x0,2));
           -param.max_vel*ones(param.numJ-1,size(x0,2));
           -param.max_front_ank_tau*ones(1,size(x0,2));
           -param.max_tau*ones(param.numJ,size(x0,2))];
           



%% define object function
prob.objective = @(x)objFun(x,param);




%% solve
iterTime = 2000;
options = optimoptions('fmincon','MaxIter',iterTime,...
    'Display','iter','GradObj','on','TolCon',1e-3,'SpecifyConstraintGradient',true,...
    'SpecifyObjectiveGradient',true,'StepTolerance',1e-10,'UseParallel',true);
prob.options = options;
prob.solver = 'fmincon';

exitflag=0;
x = x0;
prob.x0 = x0;
while(true)
     
    [x,fval,exitflag,output] = fmincon(prob);
    if(exitflag==0)
        prob.x0 = x;
       
    else
        break;
    end
    
    
end

    
[t1,~]=clock;
fileName = [num2str(t1(2),'%02d'),num2str(t1(3),'%02d'),num2str(t1(4),'%02d'),num2str(t1(5),'%02d')];
result.x = x;
result.fval=fval;
result.exitflag = exitflag;
result.output = output;
result.param = param;
result.x0=prob.x0;
result.set_iterTime = iterTime;

save(fileName,'result');
