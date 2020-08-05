%% optimal gait generating with contact invarient method
% states to solve: q, c (contact indicator)
clear;
clc;
modelName='human_7';
addpath dyn/
addpath obj/
addpath gaitCon/
addpath plotRobot/
addpath (['../',modelName,'/robotGen/'])
addpath (['../',modelName,'/robotGen/grad/'])
addpath (['../',modelName,'/robotGen/posCons/'])
addpath (['../',modelName,'/robotGen/dyn/'])
addpath (['../',modelName,'/robotGen/obj/'])
addpath (['../',modelName,'/robotGen/grf/'])
addpath (['../',modelName,'/robotGen/knee_spring/'])

%% set up the parameters
model = load(['../',modelName,'/robotGen/model']).model;
param.toe_th =-model.l_heel+0.01;
param.numJ=6;
param.gaitT = 0.5;
param.sampT = 0.01;
param.dmax =1e-3;
param.cmax=100;
param.k=2e6;
param.miu =0.8;
param.floor_h = -model.l_heel;
param.jointW = [1,1,1,1,1,1]*10;
dataLen = param.gaitT/param.sampT+1;
param.dataLen=dataLen;
k1 = 1;
k0 = 10;
param.c_matFun = @(c)eye(2)*(k0/(k1+c^2));

param.M1 = zeros(1,2*dataLen);
param.M1(1,1) = 1;
param.M1(1,end-1) =-1;
param.M2 = zeros(1,2*dataLen);
param.M2(1,2)=1;
param.M2(1,end)=-1;

fri_mat = [1, -param.miu;
    -1, -param.miu];
fri_mat = repmat({fri_mat},1,dataLen);
fri_mat = blkdiag(fri_mat{:});
param.fri_mat = fri_mat;





param.knee_stiff =76.325; % I use max moment (MVC/angle), since the stiffness of the paper is too high
param.ank_stiff=408.65;

param.joint_fri = 0.003;
param.gndclear = -model.l_heel+0.02;

param.head_h = 1.1 ; %the head should be at least 1.6m
param.hip_feet_ratio = 2;
param.hipLen=param.hip_feet_ratio*model.l_foot;
param.gaitLen = model.l_foot*3;

param.dynCoeff=0.001;
param.gndCoeff = 0;
param.max_ank_vel = Inf;
param.max_kne_vel = Inf;
param.max_hip_vel = Inf;
%% Generate nominal trajectory

q1 =60;
q2 = -10;
q3 =90-q1-q2;
q4 = -180+q3;
q5 = 1;
q6 = -180-q1-q2-q3-q4-q5;
qStart=[q1/180*pi,q2/180*pi,q3/180*pi,q4/180*pi,q5/180*pi,q6/180*pi];
q1_end = -q6;
q2_end = -q5;
q3_end = -180-q4;
q4_end = -180-q3;
q5_end = -q2;
q6_end = -q1;
qEnd = [q1_end,q2_end,q3_end,q4_end,q5_end,q6_end]*pi/180;


q = [linspace(qStart(1),qEnd(1),dataLen);
    linspace(qStart(2),qEnd(2),dataLen);
    linspace(qStart(3),qEnd(3),dataLen);
    linspace(qStart(4),qEnd(4),dataLen);
    linspace(qStart(5),qEnd(5),dataLen);
    linspace(qStart(6),qEnd(6),dataLen)];
dq = gradient(q)/param.sampT;
ddq = gradient(dq)/param.sampT;
% calculate the c value (is there contact)
c_toe = 0*ones(1,dataLen);
c_heel = 0*ones(1,dataLen);
J.toe = cell(1,dataLen);
J.heel = cell(1,dataLen);
tau= zeros(param.numJ,dataLen);
for i =1:dataLen
    cur_q = q(:,i);
    cur_dq = dq(:,i);
    cur_ddq = ddq(:,i);
    % check toe pos
    toe_h = end_y_pos(cur_q.');
    heel_h = heel_y_pos(cur_q.');
    J_toe = six_J(cur_q(1),cur_q(2),cur_q(3),cur_q(4),cur_q(5),cur_q(6));
    J_heel = six_J2(cur_q(1),cur_q(2),cur_q(3),cur_q(4),cur_q(5),cur_q(6));
    
    J.toe{1,i} = J_toe(1:2,:);
    J.heel{1,i} = J_heel(1:2,:);
    
    tau(:,i) = u_no_ext([cur_q;cur_dq],ddq(:,i),param);
    
    toe_vel = J_toe*dq(:,i);
    heel_vel = J_heel*dq(:,i);
    if(toe_h<param.toe_th)
        cur_c=param.k*(param.toe_th-toe_h)^2-param.cmax*(param.toe_th-toe_h)/param.dmax*toe_vel(2);
        if(cur_c>0)
            c_toe(1,i) = cur_c;
        end
    end
    if(heel_h<param.toe_th)
        cur_c = param.k*(param.toe_th-heel_h)^2-param.cmax*(param.toe_th-heel_h)/param.dmax*heel_vel(2);
        if(cur_c>0)
            c_heel(1,i)=cur_c;
        end
    end
end


Aeq = zeros(param.numJ+1,(param.numJ+2)*dataLen);

Aeq(1:param.numJ,1:param.numJ)=[0,0,0,0,0,1;   %start frame
    0,0,0,0,1,0;
    0,0,0,1,0,0;
    0,0,1,0,0,0;
    0,1,0,0,0,0;
    1,0,0,0,0,0];

% endframe
Aeq(1:param.numJ,end-(param.numJ+2)+1:end-2) =[1,0,0,0,0,0;
    0,1,0,0,0,0;
    0,0,1,0,0,0;
    0,0,0,1,0,0;
    0,0,0,0,1,0;
    0,0,0,0,0,1];
% initial feet is flat
Aeq(param.numJ+1,1:param.numJ)=[1,1,1,1,1,1];


prob.Aeq = Aeq;
prob.beq = [ 0;0;-pi;-pi;0;0;-pi];

%torso angle limit
Asamp = zeros(2,param.numJ+2); % create A for single frame
Asamp(1:2,1:3) = [-1,-1,-1;
    1, 1, 1];

Acell = repmat({Asamp},1,floor(dataLen));
prob.Aineq = blkdiag(Acell{:});
Bsamp = [-90/180*pi;
    110/180*pi];
prob.bineq = repmat(Bsamp,floor(param.gaitT/param.sampT+1),1);

iterTime =200;

options = optimoptions('fmincon','Algorithm','interior-point','MaxIter',floor(iterTime),'MaxFunctionEvaluations',floor(iterTime)*1000,...
    'Display','iter','TolCon',1e-8,...
    'StepTolerance',1e-15,'UseParallel',true,'DiffMinChange',0,'ScaleProblem',true,'ObjectiveLimit',1e-4);
prob.solver = 'fmincon';
prob.options = options;

prob.ub = [179/180*pi*ones(1,dataLen);
    0*ones(1,dataLen)/180*pi;
    75/180*pi*ones(1,dataLen);
    -100/180*pi*ones(1,dataLen);
    179/180*pi*ones(1,dataLen);
    -45/180*pi*ones(1,dataLen);
    Inf*ones(1,dataLen);
    Inf*ones(1,dataLen)];

prob.lb = [ones(1,dataLen)/180*pi;
    -179/180*pi*ones(1,dataLen);
    -75/180*pi*ones(1,dataLen);
    -260/180*pi*ones(1,dataLen);
    -0*ones(1,dataLen)/180*pi;
    -135/180*pi*ones(1,dataLen);
    -0.0000000001*ones(1,dataLen);
    -0.0000000001*ones(1,dataLen)];




prob.objective = @(x)Contact_inv_loss(x,param);
prob.x0 = [q;c_toe;c_heel];


[x,fval,exitflag,output] = fmincon(prob);