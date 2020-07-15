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
addpath robotGen/knee_spring
%% simulation parameter
model = load('robotGen/model').model;
param.numJ=6;
param.toe_th =-model.l_heel+0.01;
param.head_h = 1.1 ; %the head should be at least 1.6m
param.fri_coeff=1000;
param.gaitT = 0.4;
param.sampT = 0.005;
param.init_y = -model.l_heel+0.005; %initial feet height
param.heel_h = model.l_heel; %this is fix in the model parameter
param.foot_l = model.l_foot;
param.hipLen=2*model.l_foot;
param.gndclear = -5e-5;
param.jointW = [1,1,1,1,1,1];
param.knee_stiff=500;
param.ank_stiff=70;
time = 0:param.sampT:param.gaitT;


% set torque/angular velocity constraints
param.max_tau = 200;
param.max_vel = 480/180*pi;



param.mass =model.totM;
% param.max_hip_tau =1*param.mass;
% param.min_hip_tau = 0.8*param.mass;
% param.max_kne_tau = 1.5*param.mass;
% param.min_kne_tau =0.5*param.mass;
% param.max_ank_tau =2*param.mass;
% param.min_ank_tau=0.01*param.mass;

param.max_hip_tau =param.mass;
param.min_hip_tau = param.mass;
param.max_kne_tau = param.mass;
param.min_kne_tau =param.mass;
param.max_ank_tau =param.mass;
param.min_ank_tau= param.mass;
%% initialize joint pos and torque
qmax = 170/180/pi;
% q = qmax*sin((2*time/param.gaitT+randn(param.jointNum,1))*pi);
% dq = qmax*sin((2*time/param.gaitT+randn(param.jointNum,1))*pi)*10;
q1 = 70;
q2 = -5;
q3 =45;
q4 = -180;
q5 =2;


q6 = -180-q1-q2-q3-q4-q5;%+atan2d(param.heel_h,0.26);%0.26 is feet length

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

% the A matrix is define in the following way:
%     [x(0),x(1),x(2).......x(end)], one condition, one row

numCond = param.numJ+param.numJ/2; %start-end pos conditions, velocity conditions

%start-end joint condition 
Aeq = zeros(numCond,size(x0,1)*size(x0,2)); 
Aeq(1:param.numJ,1:param.numJ)=[1,1,1,1,1,0;   %start frame 
                                0,0,0,0,1,0;
                                0,0,0,1,0,0;
                                0,0,1,0,0,0;
                                0,1,0,0,0,0;
                                1,0,0,0,0,0];
Aeq(1:param.numJ,end-param.numJ*3+1:end-param.numJ*2) = [-1,0,0,0,0,0; % endframe
                                                          0,1,0,0,0,0;
                                                          0,0,1,0,0,0;
                                                          0,0,0,1,0,0;
                                                          0,0,0,0,1,0;
                                                          0,0,0,0,0,1];
                                                      
Aeq(param.numJ+1:param.numJ+3,param.numJ+1:param.numJ*2) = [1,0,0,0,0,-1;
                                                            0,1,0,0,1,0;
                                                            0,0,1,1,0,0];
                                                      
                                                      
                                                      
prob.Aeq = Aeq;
prob.beq = [-pi;0;-pi;-pi;0;0;0;0;0];

% back never bend backward
% -1*(q1+q2+q3) <-88 deg, 
%    -q3 < -1 deg
Asamp = zeros(1,3*param.numJ); % create A for single frame
Asamp(:,1:3) = [-1,-1,-1];
                % 0, 0, 1];    

Acell = repmat({Asamp},1,param.gaitT/param.sampT+1);
prob.Aineq = blkdiag(Acell{:});
Bsamp = [-90/180*pi];
         % -1/180/pi];
prob.bineq = repmat(Bsamp,param.gaitT/param.sampT+1,1); 

% upper limit and lower limit for each joints
prob.ub = [179/180*pi*ones(1,size(x0,2));
           -ones(1,size(x0,2))/180*pi;
           89/180*pi*ones(1,size(x0,2));
           -91/180*pi*ones(1,size(x0,2));
           179/180*pi*ones(1,size(x0,2));
           -10/180*pi*ones(1,size(x0,2));
           param.max_vel*ones(param.numJ,size(x0,2));
           param.max_ank_tau*ones(1,size(x0,2));
           param.max_kne_tau*ones(1,size(x0,2));
           param.max_hip_tau*ones(1,size(x0,2));
           param.min_hip_tau*ones(1,size(x0,2));
           param.min_kne_tau*ones(1,size(x0,2));
           param.max_ank_tau*ones(1,size(x0,2))];
prob.lb = [ones(1,size(x0,2))/180*pi;
           -179/180*pi*ones(1,size(x0,2));
           -89/180*pi*ones(1,size(x0,2));
           -259/180*pi*ones(1,size(x0,2));
           -ones(1,size(x0,2))/180*pi;
           -125/180*pi*ones(1,size(x0,2));
           -param.max_vel*ones(param.numJ,size(x0,2));
           -param.min_ank_tau*ones(1,size(x0,2));
           -param.min_kne_tau*ones(1,size(x0,2));
           -param.min_hip_tau*ones(1,size(x0,2));
           -param.max_hip_tau*ones(1,size(x0,2));
           -param.max_kne_tau*ones(1,size(x0,2));
           -param.min_ank_tau*ones(1,size(x0,2))];


%% define object function
prob.objective = @(x)objFun(x,param);




%% solve
iterTime = 2000;
options = optimoptions('fmincon','MaxIter',iterTime,'MaxFunctionEvaluations',iterTime*5,...
    'Display','iter','GradObj','on','TolCon',1e-8,'SpecifyConstraintGradient',true,...
    'SpecifyObjectiveGradient',true,'StepTolerance',1e-10,'UseParallel',true);
prob.options = options;
prob.solver = 'fmincon';

% exitflag=0;
% x = x0;
prob.x0 = x0;
% while(true)
%      
%     [x,fval,exitflag,output] = fmincon(prob);
%     if(exitflag==0)
%         prob.x0 = x;
%         
%        
%     else
%         break;
%     end
%     
%     
% end
[x,fval,exitflag,output] = fmincon(prob);
    
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
msgbox('optimization done');