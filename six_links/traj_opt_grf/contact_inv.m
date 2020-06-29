%% Optimal gait trajectory generating with contact invarient method
clear;
modelName='human_3';
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
param.k1 = 1e-3;
param.k0 = 1e-2;
param.jointW = [1,1,1,1,1,1];

param.knee_stiff =76.325; % I use max moment (MVC/angle), since the stiffness of the paper is too high

param.ank_stiff=408.65;
param.joint_fri = 0.003;

param.head_h = 1.1 ; %the head should be at least 1.6m
param.hip_feet_ratio = 2;
param.hipLen=param.hip_feet_ratio*model.l_foot;
param.gaitLen = model.l_foot*3;

param.max_ank_vel = Inf;
param.max_kne_vel = Inf;
param.max_hip_vel = Inf;
%% Generate nominal trajectory

q1 =70;
q2 = -1;
q3 =90-q1-q2;
q4 = -80-q1;
q5 = 3;
q6 = -180-q1-q2-q3-q4-q5;
qStart=[q1/180*pi,q2/180*pi,q3/180*pi,q4/180*pi,q5/180*pi,q6/180*pi];
q1_end = 180+q1+q2+q3+q4+q5;
q2_end = -q5;
q3_end = -180-q4;
q4_end = -q1_end-90;
q5_end = -q2;
q6_end = -180-q1_end-q2_end-q3_end-q4_end-q5_end;
qEnd = [q1_end,q2_end,q3_end,q4_end,q5_end,q6_end]*pi/180;

dataLen = param.gaitT/param.sampT;
q = [linspace(qStart(1),qEnd(1),dataLen+1);
    linspace(qStart(2),qEnd(2),dataLen+1);
    linspace(qStart(3),qEnd(3),dataLen+1);
    linspace(qStart(4),qEnd(4),dataLen+1);
    linspace(qStart(5),qEnd(5),dataLen+1);
    linspace(qStart(6),qEnd(6),dataLen+1)];

dq =[(q(1,2)-q(1,1))/param.sampT*ones(1,dataLen+1);
    (q(2,2)-q(2,1))/param.sampT*ones(1,dataLen+1);
    (q(3,2)-q(3,1))/param.sampT*ones(1,dataLen+1);
    (q(4,2)-q(4,1))/param.sampT*ones(1,dataLen+1);
    (q(5,2)-q(5,1))/param.sampT*ones(1,dataLen+1);
    (q(6,2)-q(6,1))/param.sampT*ones(1,dataLen+1)];
ddq = zeros(6,dataLen+1);

u = zeros(param.numJ,dataLen);
f = zeros(4,dataLen);
tau = zeros(param.numJ,dataLen);

% we build friction matrix outside since it will never change
% constraints:
% fri_mat*f<=0
fri_mat = [1, -param.miu;
    -1, -param.miu];
fri_mat = repmat({fri_mat},1,dataLen);
fri_mat = blkdiag(fri_mat{:});

% function we need for generating W matrix from c values
c_matFun = @(c)eye(2)*c;


% creating fmincon solver parameters
prob.ub = [179/180*pi*ones(1,dataLen+1);
    0*ones(1,dataLen+1)/180*pi;
    75/180*pi*ones(1,dataLen+1);
    -100/180*pi*ones(1,dataLen+1);
    179/180*pi*ones(1,dataLen+1);
    -45/180*pi*ones(1,dataLen+1);
    param.max_ank_vel*ones(1,dataLen+1);
    param.max_kne_vel*ones(1,dataLen+1);
    param.max_hip_vel*ones(1,dataLen+1);
    param.max_hip_vel*ones(1,dataLen+1);
    param.max_kne_vel*ones(1,dataLen+1);
    param.max_ank_vel*ones(1,dataLen+1)];

prob.lb = [ones(1,dataLen+1)/180*pi;
    -179/180*pi*ones(1,dataLen+1);
    -75/180*pi*ones(1,dataLen+1);
    -260/180*pi*ones(1,dataLen+1);
    -0*ones(1,dataLen+1)/180*pi;
    -135/180*pi*ones(1,dataLen+1);
    -param.max_ank_vel*ones(1,dataLen+1);
    -param.max_kne_vel*ones(1,dataLen+1);
    -param.max_hip_vel*ones(1,dataLen+1);
    -param.max_hip_vel*ones(1,dataLen+1);
    -param.max_kne_vel*ones(1,dataLen+1);
    -param.max_ank_vel*ones(1,dataLen+1)];

% conditions for start/end angles, velocity
Aeq = zeros(param.numJ*2,param.numJ*2*(dataLen+1));

Aeq(1:param.numJ,1:param.numJ*2)=[1,1,1,1,1,0, 0,0,0,0,0,0;   %start frame
                                    0,0,0,0,1,0, 0,0,0,0,0,0;
                                    0,0,0,1,0,0, 0,0,0,0,0,0;
                                    0,0,1,0,0,0, 0,0,0,0,0,0;
                                    0,1,0,0,0,0, 0,0,0,0,0,0;
                                    1,0,0,0,0,0, 0,0,0,0,0,0];
                                    
% endframe
Aeq(1:param.numJ,end-param.numJ*2+1:end) =  [-1,0,0,0,0,0, 0,0,0,0,0,0;
                                              0,1,0,0,0,0, 0,0,0,0,0,0;
                                              0,0,1,0,0,0, 0,0,0,0,0,0;
                                              0,0,0,1,0,0, 0,0,0,0,0,0;
                                              0,0,0,0,1,0, 0,0,0,0,0,0;
                                              0,0,0,0,0,1, 0,0,0,0,0,0];
%velocity
%start frame
Aeq(param.numJ+1:param.numJ*2,1:param.numJ*2) = [0,0,0,0,0,0,  1,0,0,0,0,0;
                                                 0,0,0,0,0,0,  0,1,0,0,0,0;
                                                 0,0,0,0,0,0,  0,0,1,0,0,0;
                                                 0,0,0,0,0,0,  0,0,0,1,0,0;
                                                 0,0,0,0,0,0,  0,0,0,0,1,0;
                                                 0,0,0,0,0,0,  0,0,0,0,0,1];
%end frame
Aeq(param.numJ+1:param.numJ*2,end-param.numJ*2+1:end)=[0,0,0,0,0,0,  0,0,0,0,0,1;
                                                       0,0,0,0,0,0,  0,0,0,0,1,0;
                                                       0,0,0,0,0,0,  0,0,0,1,0,0;
                                                       0,0,0,0,0,0,  0,0,1,0,0,0;
                                                       0,0,0,0,0,0,  0,1,0,0,0,0;
                                                       0,0,0,0,0,0,  1,0,0,0,0,0];
% initial feet is flat
Aeq(param.numJ*2+1,1:param.numJ*2)=[1,1,1,1,1,1, 0,0,0,0,0,0];
                                                   
                                                   
prob.Aeq = Aeq;
prob.beq = [-pi;0;-pi;-pi;0;0;0;0;0;0;0;0;-pi];                                                
   
%torso angle limit
Asamp = zeros(4,param.numJ*2); % create A for single frame
Asamp(1:2,1:3) = [-1,-1,-1;
                   1, 1, 1];
% torso speed limit
Asamp(3:4,1+param.numJ:3+param.numJ)=[-1,-1,-1;
                                       1, 1, 1];

Acell = repmat({Asamp},1,floor(dataLen+1));
prob.Aineq = blkdiag(Acell{:});
Bsamp = [-90/180*pi;
    110/180*pi;
    10/180*pi/param.sampT;
    10/180*pi/param.sampT];
prob.bineq = repmat(Bsamp,floor(param.gaitT/param.sampT+1),1);

iterTime=15000;
options = optimoptions('fmincon','Algorithm','interior-point','MaxIter',iterTime,'MaxFunctionEvaluations',iterTime*5,...
    'Display','iter','GradObj','on','TolCon',1e-3,'SpecifyConstraintGradient',true,...
    'SpecifyObjectiveGradient',true,'StepTolerance',1e-15,'UseParallel',true,'DiffMinChange',0,'ScaleProblem',true,'Display','final');
prob.options = options;
prob.solver = 'fmincon';
% optimization started
optCont = true;
tot_count =1;
while(optCont)
    % calculate the c value (is there contact)
    c.c_toe = param.k0/param.k1*ones(1,dataLen);
    c.c_heel = param.k0/param.k1*ones(1,dataLen);
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
        
        G = six_G(cur_q(1),cur_q(2),cur_q(3),cur_q(4),cur_q(5),cur_q(6));
        V = six_V(cur_q(2),cur_q(3),cur_q(4),cur_q(5),cur_q(6),cur_dq(1),cur_dq(2),cur_dq(3),cur_dq(4),cur_dq(5),cur_dq(6));
        M = six_M(cur_q(2),cur_q(3),cur_q(4),cur_q(5),cur_q(6));
        tau(:,i) = M*ddq(:,i)+V.'+G.';
        
        toe_vel = J_toe*dq(:,i);
        heel_vel = J_heel*dq(:,i);
        if(toe_h<param.toe_th)
            cur_c=param.k*(param.toe_th-toe_h)^2-param.cmax*(param.toe_th-toe_h)/param.dmax*toe_vel(2);
            if(cur_c>0)
                c.c_toe(1,i) = param.k0/(cur_c^2+param.k1);
            end
            
        end
        if(heel_h<param.toe_th)
            cur_c = param.k*(param.toe_th-heel_h)^2-param.cmax*(param.toe_th-heel_h)/param.dmax*heel_vel(2);
            if(cur_c>0)
                c.c_heel(1,i)=param.k0/(cur_c^2+param.k1);
                
            end
        end
        
    end
    %% solve u,f based on q,dq,ddq
    
    tau = reshape(tau,[size(tau,1)*size(tau,2),1]);
    JT_toe_cell = cellfun(@transpose,J.toe,'UniformOutput',false);
    JT_heel_cell = cellfun(@transpose,J.heel,'UniformOutput',false);
    JT_toe_big = blkdiag(JT_toe_cell{:});
    JT_heel_big = blkdiag(JT_heel_cell{:});
    R = diag(repmat(param.jointW,1,dataLen));
    
    W_toe_cell = arrayfun(c_matFun,c.c_toe,'UniformOutput',false);
    W_toe = blkdiag(W_toe_cell{:});
    W_heel_cell = arrayfun(c_matFun,c.c_heel,'UniformOutput',false);
    W_heel = blkdiag(W_heel_cell{:});
    
    cvx_begin quiet
        variables f_toe(2*dataLen,1) f_heel(2*dataLen,1) u(param.numJ*dataLen,1)
        minimize( norm(JT_toe_big*f_toe+JT_heel_big*f_heel+u-tau)+f_toe.'*W_toe*f_toe+f_heel.'*W_heel*f_heel+u.'*R*u)
        subject to
            fri_mat*f_toe<=0;
            fri_mat*f_heel<=0;
    cvx_end
    % reshape to fit the format later
    f_toe = reshape(f_toe,[2,dataLen]);
    f_heel= reshape(f_heel,[2,dataLen]);
    u = reshape(u,[param.numJ,dataLen]);
    % duplicate the first column to the last with changing direction
    
    f_toe = [f_toe,f_toe(:,1)];
    f_heel = [f_heel,f_heel(:,1)];
    u = [u, [-u(4:6,1);-u(1:3,1)]];
    
    %% with the previously got u, f_toe, f_heel, now we solve for the trajectory
    
    % since this problem will be non-convex and nonlinear, we use fmincon
    % to solve this problem
    prob.objective = @(x)dynObj(x,param,u,f_toe,f_heel);
    prob.x0 = [q;dq];
   
    
    
    [x,fval,exitflag,output] = fmincon(prob);
    
    q=x(1:param.numJ,:);
    dq = x(1:param.numJ,:);
    ddq = dq(:,2:end)-dq(:,1:end-1)/param.sampT;
    
    
    tot_count = tot_count+1;
    if(tot_count==50)
        break;
    end
    optCont = true;
end



























