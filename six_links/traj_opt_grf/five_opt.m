%% Calculate the optimized trajectories for 5 link biped model with GRF
clear;
modelName='human_3';
warning on verbose
%add share functions
addpath dyn/
addpath obj/
addpath gaitCon/
addpath plotRobot/
dbstop if error
addpath hessian/

addpath (['../',modelName,'/robotGen/'])
addpath (['../',modelName,'/robotGen/grad/'])
addpath (['../',modelName,'/robotGen/posCons/'])
addpath (['../',modelName,'/robotGen/dyn/'])
addpath (['../',modelName,'/robotGen/obj/'])
addpath (['../',modelName,'/robotGen/grf/'])
addpath (['../',modelName,'/robotGen/knee_spring/'])

%% simulation parameter
model = load(['../',modelName,'/robotGen/model']).model;
param.numJ=6;
param.toe_th =-model.l_heel+0.01;
param.head_h = 1.1 ; %the head should be at least 1.6m
param.fri_coeff=0;
param.gaitT = 0.5;
param.sampT = 0.01;
%param.init_y = -model.l_heel+0.01; %initial feet height
param.heel_h = model.l_heel; %this is fix in the model parameter
param.foot_l = model.l_foot;
param.dmax =1e-3;
param.cmax=1250;
param.k=2e6;
param.us=0.8;
param.ud=0.6;

param.joint_fri = 0.003;

param.hip_feet_ratio = 2;
param.hipLen=param.hip_feet_ratio*model.l_foot;
param.gaitLen = model.l_foot*3;


param.gndclear = -model.l_heel+0.02;
param.jointW = [1,1,1,1,1,1];
% param.knee_stiff=1906.2/2; %patellar tendon, since we have two in series at the knee, the spring constant is half
param.knee_stiff =76.325; % I use max moment (MVC/angle), since the stiffness of the paper is too high

param.ank_stiff=408.65;

time = 0:param.sampT:param.gaitT;


% set torque/angular velocity constraints

% 
% param.max_ank_vel = 2*360/180*pi;
% param.max_kne_vel = 2*360/180*pi;
% param.max_hip_vel = 2*270/180*pi;

param.max_ank_vel = Inf;
param.max_kne_vel = Inf;
param.max_hip_vel = Inf;


param.mass =model.totM;
% param.max_hip_tau =2*param.mass;
% param.min_hip_tau = 2*param.mass;
% param.max_kne_tau = 2*param.mass;
% param.min_kne_tau =1.5*param.mass;
% param.max_ank_tau =2*param.mass;
% param.min_ank_tau=0.1*param.mass;

param.max_Fy = Inf;
param.max_Fx = Inf;
param.min_Fx = Inf;
param.min_Fy = Inf;
% 
param.max_hip_tau =Inf;
param.min_hip_tau = Inf;
param.max_kne_tau = Inf;
param.min_kne_tau =Inf;
param.max_ank_tau =Inf;
param.min_ank_tau= Inf;

% param.max_hip_tau =param.mass;
% param.min_hip_tau = param.mass;
% param.max_kne_tau = param.mass*2;
% param.min_kne_tau =param.mass*2;
% param.max_ank_tau =param.mass*2;
% param.min_ank_tau= param.mass*2;

% param.max_hip_tau =Inf;
% param.min_hip_tau = Inf;
% param.max_kne_tau = Inf;
% param.min_kne_tau = Inf;
% param.max_ank_tau = Inf;
% param.min_ank_tau= Inf;
%% initialize joint pos and torque

% q = qmax*sin((2*time/param.gaitT+randn(param.jointNum,1))*pi);
% dq = qmax*sin((2*time/param.gaitT+randn(param.jointNum,1))*pi)*10;
% q1 = 70;
% q2 = -5;
% q3 =45;
% q4 = -180;
% q5 =10;

num_1 = floor(length(time)/4);
num_2 = floor((length(time)-num_1)/2);
num_3 = length(time)-num_1-num_2;

q1 =70;
q2 = -1;
q3 =90-q1-q2;
q4 = -80-q1;
q5 = 3;

% 
% q1 =90;
% q2 = 0;
% q3 =90-q1-q2;
% q4 = -180;
% q5 = 0;
% q6=-90;

q6 = -180-q1-q2-q3-q4-q5;%+atan2d(param.heel_h,0.26);%0.26 is feet length

qStart=[q1/180*pi,q2/180*pi,q3/180*pi,q4/180*pi,q5/180*pi,q6/180*pi];
% qEnd = [pi+qStart(1)+qStart(2)+qStart(3)+qStart(4)+qStart(5),...
%         -qStart(5),...
%         -qStart(4)-pi,...
%         -pi-qStart(3),...
%         -qStart(2),...
%         -qStart(1)];
q1_mid_1 = 110;
q2_mid_1 = -15;
q3_mid_1 = 90-q1_mid_1-q2_mid_1;
q4_mid_1 = -220;
q5_mid_1 = 90;
q6_mid_1 = -90;

qMid_1 = [q1_mid_1,q2_mid_1,q3_mid_1,q4_mid_1,q5_mid_1,q6_mid_1]*pi/180;

q1_mid_2 = 100;
q2_mid_2 = -5;
q3_mid_2 = 90-q1_mid_2-q2_mid_2;
q4_mid_2 = -210;
q5_mid_2 = 45;
q6_mid_2 = -90;
qMid_2 = [q1_mid_2,q2_mid_2,q3_mid_2,q4_mid_2,q5_mid_2,q6_mid_2]*pi/180;


q1_end = 180+q1+q2+q3+q4+q5;
q2_end = -q5;
q3_end = -180-q4;
q4_end = -q1_end-90;
q5_end = -q2;
q6_end = -180-q1_end-q2_end-q3_end-q4_end-q5_end;
qEnd = [q1_end,q2_end,q3_end,q4_end,q5_end,q6_end]*pi/180;

q = [linspace(qStart(1),qMid_1(1),num_1),linspace(qMid_1(1),qMid_2(1),num_2),linspace(qMid_2(1),qEnd(1),num_3);
     linspace(qStart(2),qMid_1(2),num_1),linspace(qMid_1(2),qMid_2(2),num_2),linspace(qMid_2(2),qEnd(2),num_3);
     linspace(qStart(3),qMid_1(3),num_1),linspace(qMid_1(3),qMid_2(3),num_2),linspace(qMid_2(3),qEnd(3),num_3);
     linspace(qStart(4),qMid_1(4),num_1),linspace(qMid_1(4),qMid_2(4),num_2),linspace(qMid_2(4),qEnd(4),num_3);
     linspace(qStart(5),qMid_1(5),num_1),linspace(qMid_1(5),qMid_2(5),num_2),linspace(qMid_2(5),qEnd(5),num_3);
     linspace(qStart(6),qMid_1(6),num_1),linspace(qMid_1(6),qMid_2(6),num_2),linspace(qMid_2(6),qEnd(6),num_3)];

% q = [linspace(91/180*pi,91/180*pi,length(time));
%      linspace(-1/180*pi,-1/180*pi,length(time));
%      linspace(0/180*pi,0/180*pi,length(time));
%      linspace(-181/180*pi,-181/180*pi,length(time));
%      linspace(1/180*pi,1/180*pi,length(time));
%      linspace(-90/180*pi,-90/180*pi,length(time))];

dq = [0,(q(1,2:end)-q(1,1:end-1))/param.sampT;
      0,(q(2,2:end)-q(2,1:end-1))/param.sampT;
      0,(q(3,2:end)-q(3,1:end-1))/param.sampT;
      0,(q(4,2:end)-q(4,1:end-1))/param.sampT;
      0,(q(5,2:end)-q(5,1:end-1))/param.sampT;
      0,(q(6,2:end)-q(6,1:end-1))/param.sampT];

% the initial guess of u are calculated by forward dynamics
u = zeros(size(q,1),size(q,2));
for i=1:size(u,2)
    cur_q = q(:,i);
    cur_dq =dq(:,i);
    % although not exactly, but for here I assume ddq is zero (dq is almost
    % steady)
    V = six_V(cur_q(2),cur_q(3),cur_q(4),cur_q(5),cur_q(6),cur_dq(1),cur_dq(2),cur_dq(3),cur_dq(4),cur_dq(5),cur_dq(6));
    G = six_G(cur_q(1),cur_q(2),cur_q(3),cur_q(4),cur_q(5),cur_q(6));
    u(:,i) = V.'+G.';
    
end



Fext_toe = [zeros(1,length(q));ones(1,length(q))];
Fext_heel = [zeros(1,length(q));ones(1,length(q))];

slack_var = zeros(2,length(q));
x0 = [q;dq;u;Fext_toe;Fext_heel;slack_var];

% tori= linspace(0,1,size(x0,2));
% x0_temp = load('x0_val5').x;
% t_temp = linspace(0,1,size(x0_temp,2));
% x0=interp1(t_temp,x0_temp.',tori).';
% dyndq = @(x,tauext)dyn_dq(x,tauext,jointNum,toe_th);
% rowfun(dyndq,table(x,ext_tau))

%% define constraints
%  dynamic+gait length + vertical pos
prob.nonlcon = @(x)five_link_nonlcon(x,param);

% linear constraints

% the A matrix is define in the following way:
%     [x(0),x(1),x(2).......x(end)], one condition, one row

numCond = 19; %start-end pos conditions, velocity conditions
numS = param.numJ*3+4+2;
%start-end joint condition 
%position
Aeq = zeros(numCond,size(x0,1)*size(x0,2)); 
Aeq(1:param.numJ+1,1:param.numJ*3)=[1,1,1,1,1,0, 0,0,0,0,0,0,0,0,0,0,0,0;   %start frame 
                                    0,0,0,0,1,0, 0,0,0,0,0,0,0,0,0,0,0,0;
                                    0,0,0,1,0,0, 0,0,0,0,0,0,0,0,0,0,0,0;
                                    0,0,1,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0;
                                    0,1,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0;
                                    1,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0;
                                    1,1,1,1,1,1, 0,0,0,0,0,0,0,0,0,0,0,0];
% endframe
Aeq(1:param.numJ,end-numS+1:end-numS+param.numJ*3) =   [-1,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0; 
                                                         0,1,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0;
                                                         0,0,1,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0;
                                                         0,0,0,1,0,0, 0,0,0,0,0,0,0,0,0,0,0,0;
                                                         0,0,0,0,1,0, 0,0,0,0,0,0,0,0,0,0,0,0;
                                                         0,0,0,0,0,1, 0,0,0,0,0,0,0,0,0,0,0,0];
%velocity    
%start frame
Aeq(param.numJ+2:param.numJ*2+1,1:param.numJ*3) = [0,0,0,0,0,0,  1,0,0,0,0,0  ,0,0,0,0,0,0;
                                                 0,0,0,0,0,0,  0,1,0,0,0,0  ,0,0,0,0,0,0;
                                                 0,0,0,0,0,0,  0,0,1,0,0,0  ,0,0,0,0,0,0;
                                                 0,0,0,0,0,0,  0,0,0,1,0,0  ,0,0,0,0,0,0;
                                                 0,0,0,0,0,0,  0,0,0,0,1,0  ,0,0,0,0,0,0;
                                                 0,0,0,0,0,0,  0,0,0,0,0,1  ,0,0,0,0,0,0];
%end frame
Aeq(param.numJ+2:param.numJ*2+1,end-numS+1:end-numS+param.numJ*3)=[0,0,0,0,0,0,  0,0,0,0,0,1  ,0,0,0,0,0,0;
                                                       0,0,0,0,0,0,  0,0,0,0,1,0  ,0,0,0,0,0,0;
                                                       0,0,0,0,0,0,  0,0,0,1,0,0  ,0,0,0,0,0,0;
                                                       0,0,0,0,0,0,  0,0,1,0,0,0  ,0,0,0,0,0,0;
                                                       0,0,0,0,0,0,  0,1,0,0,0,0  ,0,0,0,0,0,0;
                                                       0,0,0,0,0,0,  1,0,0,0,0,0  ,0,0,0,0,0,0];
%toeque
Aeq(param.numJ*2+2:param.numJ*3+1,1:param.numJ*3) = [0,0,0,0,0,0,  0,0,0,0,0,0  ,1,0,0,0,0,0;
                                                   0,0,0,0,0,0,  0,0,0,0,0,0  ,0,1,0,0,0,0;
                                                   0,0,0,0,0,0,  0,0,0,0,0,0  ,0,0,1,0,0,0;
                                                   0,0,0,0,0,0,  0,0,0,0,0,0  ,0,0,0,1,0,0;
                                                   0,0,0,0,0,0,  0,0,0,0,0,0  ,0,0,0,0,1,0;
                                                   0,0,0,0,0,0,  0,0,0,0,0,0  ,0,0,0,0,0,1];

Aeq(param.numJ*2+2:param.numJ*3+1,end-numS+1:end-numS+param.numJ*3)=[0,0,0,0,0,0,  0,0,0,0,0,0  ,0,0,0,0,0,1;
                                                         0,0,0,0,0,0,  0,0,0,0,0,0  ,0,0,0,0,1,0;
                                                         0,0,0,0,0,0,  0,0,0,0,0,0  ,0,0,0,1,0,0;
                                                         0,0,0,0,0,0,  0,0,0,0,0,0  ,0,0,1,0,0,0;
                                                         0,0,0,0,0,0,  0,0,0,0,0,0  ,0,1,0,0,0,0;
                                                         0,0,0,0,0,0,  0,0,0,0,0,0  ,1,0,0,0,0,0];     
                                                                         
                                                                         
% initial feet must be flat

                                     
                                                      
                                                      
prob.Aeq = Aeq;
prob.beq = [-pi;0;-pi;-pi;0;0;-pi;0;0;0;0;0;0;0;0;0;0;0;0];

% back never bend backward
% -1*(q1+q2+q3) <-88 deg, 
%    -q3 < -1 deg
% abs(dq1+dq2+dq3) <10 deg/sec
% Fy_toe>0, Fy_heel>0
% s_toe >0, s_heel>0

Asamp = zeros(8,numS); % create A for single frame
Asamp(1:2,1:3) = [-1,-1,-1;
                   1, 1, 1];
Asamp(3:4,1+param.numJ:3+param.numJ)=[-1,-1,-1;
                                       1, 1, 1];
Asamp(5:6,param.numJ*3+1:param.numJ*3+4) = [0,-1,0,0;
                                            0,0,0,-1];
Asamp(7:8,param.numJ*3+5:param.numJ*3+6)=[-1,0;
                                          0,-1];

Acell = repmat({Asamp},1,floor(param.gaitT/param.sampT+1));
prob.Aineq = blkdiag(Acell{:});
Bsamp = [-90/180*pi;
          110/180*pi;
          10/180*pi/param.sampT;
          10/180*pi/param.sampT;
          0;
          0;
          0;
          0];
prob.bineq = repmat(Bsamp,floor(param.gaitT/param.sampT+1),1); 

% upper limit and lower limit for each joints
prob.ub = [179/180*pi*ones(1,size(x0,2));
           0*ones(1,size(x0,2))/180*pi;
           75/180*pi*ones(1,size(x0,2));
           -100/180*pi*ones(1,size(x0,2));
           179/180*pi*ones(1,size(x0,2));
           -45/180*pi*ones(1,size(x0,2));
           param.max_ank_vel*ones(1,size(x0,2));
           param.max_kne_vel*ones(1,size(x0,2));
           param.max_hip_vel*ones(1,size(x0,2));
           param.max_hip_vel*ones(1,size(x0,2));
           param.max_kne_vel*ones(1,size(x0,2));
           param.max_ank_vel*ones(1,size(x0,2));
           param.min_ank_tau*ones(1,size(x0,2));
           param.max_kne_tau*ones(1,size(x0,2));
           param.max_hip_tau*ones(1,size(x0,2));
           param.min_hip_tau*ones(1,size(x0,2));
           param.min_kne_tau*ones(1,size(x0,2));
           param.max_ank_tau*ones(1,size(x0,2));
           param.max_Fx*ones(1,size(x0,2));
           param.max_Fy*ones(1,size(x0,2));
           param.max_Fx*ones(1,size(x0,2));
           param.max_Fy*ones(1,size(x0,2));
           ones(2,size(x0,2))];
prob.lb = [ones(1,size(x0,2))/180*pi;
           -179/180*pi*ones(1,size(x0,2));
           -75/180*pi*ones(1,size(x0,2));
           -260/180*pi*ones(1,size(x0,2));
           -0*ones(1,size(x0,2))/180*pi;
           -135/180*pi*ones(1,size(x0,2));
           -param.max_ank_vel*ones(1,size(x0,2));
           -param.max_kne_vel*ones(1,size(x0,2));
           -param.max_hip_vel*ones(1,size(x0,2));
           -param.max_hip_vel*ones(1,size(x0,2));
           -param.max_kne_vel*ones(1,size(x0,2));
           -param.max_ank_vel*ones(1,size(x0,2));
           -param.max_ank_tau*ones(1,size(x0,2));
           -param.min_kne_tau*ones(1,size(x0,2));
           -param.min_hip_tau*ones(1,size(x0,2));
           -param.max_hip_tau*ones(1,size(x0,2));
           -param.max_kne_tau*ones(1,size(x0,2));
           -param.min_ank_tau*ones(1,size(x0,2));
           -param.max_Fx*ones(1,size(x0,2));
           -0*ones(1,size(x0,2)); %I just add some values to make it easier to calculate, after all optimal solution should gave us zero
           -param.max_Fx*ones(1,size(x0,2));
           -0*ones(1,size(x0,2));
            -0.0001*ones(2,size(x0,2))];


%% define object function
prob.objective = @(x)objFun(x,param);




%% solve

prob.x0 = x0;

iterTime =10000;
% options = optimoptions('fmincon','MaxIter',iterTime,'MaxFunctionEvaluations',iterTime*5,...
%     'Display','iter','GradObj','on','TolCon',1e-8,'SpecifyConstraintGradient',true,...
%     'SpecifyObjectiveGradient',true,'StepTolerance',1e-10,'UseParallel',true,'SubproblemAlgorithm' ,'factorization')%,'HessianFcn',@(x,lambda)hessianfcn(x,lambda,param));%,'ScaleProblem',true);%'HessianFcn',@hessianfcn


options = optimoptions('fmincon','Algorithm','interior-point','MaxIter',iterTime,'MaxFunctionEvaluations',iterTime*5,...
    'Display','iter','GradObj','on','TolCon',1e-8,'SpecifyConstraintGradient',true,...
    'SpecifyObjectiveGradient',true,'StepTolerance',1e-15,'UseParallel',true,'DiffMinChange',0,'ScaleProblem',true);%,'HessianFcn',@(x,lambda)hessianfcn(x,lambda,param));%,'ScaleProblem',true);%'HessianFcn',@hessianfcn


prob.options = options;
prob.solver = 'fmincon';

% exitflag=0;
% x = x0;

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

% ms = MultiStart('FunctionTolerance',2e-4,'UseParallel',true);
% gs = GlobalSearch(ms);
% [x,fval,exitflag,output,solutions] =run(gs,prob);


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

save(['../',modelName,'/',fileName],'result');
disp(['file name: ',modelName,'-',fileName]);
msgbox(['optimization done',num2str(exitflag)]);