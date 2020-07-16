%% Optimal gait trajectory generating with contact invarient method
clear all;
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
param.init_y = -model.l_heel;
param.numJ=6;
param.gaitT = 0.5;
param.sampT = 0.01;
param.dmax =1e-3;
param.cmax=1500;
param.k=2e6;
param.miu =0.8;
param.k1 = 0.001;
param.k0 = 0.01;
param.jointW = [100,100,100,100,100,100];

param.knee_stiff =76.325; % I use max moment (MVC/angle), since the stiffness of the paper is too high
param.ank_stiff=408.65;
param.joint_fri = 0.003;


% 
% param.knee_stiff =0; % I use max moment (MVC/angle), since the stiffness of the paper is too high
% param.ank_stiff=0;
% param.joint_fri = 0;

param.head_h = 1.1 ; %the head should be at least 1.6m
param.hip_feet_ratio = 2;
param.hipLen=param.hip_feet_ratio*model.l_foot;
param.gaitLen = model.l_foot*3;


% obj coeff
param.dynCoeff=1;
param.gndCoeff = 0;
param.initPos_w = 1;
param.init_vel_w = 1;
param.end_pos_w = 1;


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

dataLen = param.gaitT/param.sampT+1;
param.dataLen=dataLen;
q = [linspace(qStart(1),qEnd(1),dataLen);
    linspace(qStart(2),qEnd(2),dataLen);
    linspace(qStart(3),qEnd(3),dataLen);
    linspace(qStart(4),qEnd(4),dataLen);
    linspace(qStart(5),qEnd(5),dataLen);
    linspace(qStart(6),qEnd(6),dataLen)];


prob.x0=q;
% function we need for generating W matrix from c values
c_matFun = @(c)eye(2)*c;
param.c_matFun = c_matFun;
param.smooth_vel =0.01;
% creating fmincon solver parameters
prob.ub = [135/180*pi*ones(1,dataLen);
    0*ones(1,dataLen)/180*pi;
    75/180*pi*ones(1,dataLen);
    -100/180*pi*ones(1,dataLen);
    179/180*pi*ones(1,dataLen);
    -45/180*pi*ones(1,dataLen)];

prob.lb = [45*ones(1,dataLen)/180*pi;
    -179/180*pi*ones(1,dataLen);
    -75/180*pi*ones(1,dataLen);
    -260/180*pi*ones(1,dataLen);
    -0*ones(1,dataLen)/180*pi;
    -135/180*pi*ones(1,dataLen)];

% conditions for start/end angles, velocity, and initial feet is flat
Aeq = zeros(7,param.numJ*(dataLen));

Aeq(1:param.numJ+1,1:param.numJ)=[0,0,0,0,0,1;   %start frame
                                0,0,0,0,1,0;
                                0,0,0,1,0,0;
                                0,0,1,0,0,0;
                                0,1,0,0,0,0;
                                1,0,0,0,0,0;
                                1,1,1,1,1,1];
                                    
% endframe
Aeq(1:param.numJ,end-param.numJ+1:end) =  [1,0,0,0,0,0;
                                           0,1,0,0,0,0;
                                           0,0,1,0,0,0;
                                           0,0,0,1,0,0;
                                           0,0,0,0,1,0;
                                           0,0,0,0,0,1];

                                                   

param.Aeq = Aeq;
param.beq = [0;0;-pi;-pi;0;0;-pi];
   
%torso angle limit
Asamp = zeros(2,param.numJ); % create A for single frame
Asamp(1:2,1:3) = [-1,-1,-1;
                   1, 1, 1];

Acell = repmat({Asamp},1,floor(dataLen));
prob.Aineq = blkdiag(Acell{:});
Bsamp = [-90/180*pi;
        110/180*pi];
prob.bineq = repmat(Bsamp,dataLen,1);

iterTime=1000;
stepTol = 1e-10;
options = optimoptions('fmincon','Algorithm','interior-point','MaxIter',iterTime,'MaxFunctionEvaluations',iterTime*1000,...
    'Display','iter','TolCon',1e-3,'GradObj','on','GradConstr','on',...
    'StepTolerance',stepTol,'UseParallel',true,'DiffMinChange',0,'ScaleProblem',true,'ObjectiveLimit',1e-4);
prob.options = options;
prob.solver = 'fmincon';


% optimization started
optCont = true;
tot_count =1;
t_ext = 0:param.sampT:param.gaitT;
% t_over_samp = 0:param.sampT/10:param.gaitT;
% over_samp_len = size(t_over_samp,2);
% cvx_solver sedumi
% we build friction matrix outside since it will never change
% constraints:
% fri_mat*f<=0
fri_mat = [1, -param.miu;
    -1, -param.miu];
fri_mat = repmat({fri_mat},1,dataLen);
fri_mat = blkdiag(fri_mat{:});

% while the above constraint is applied on each time step, when calculating
% optimization, we need to work on (f1+2*f2+f3/4, thus, we need a matrix to
% scale f_toe, f_heel
scal_mat_cell = 0.25*[1,0,2,0,1,0;
                      0,1,0,2,0,1];

scal_mat = zeros(dataLen-2,dataLen);
for i=1:dataLen-2
    scal_mat((i-1)*2+1:(i-1)*2+2,(i-1)*2+1:(i-1)*2+6)=scal_mat_cell;
end

scal_mat_u_cell = 0.25*[1,0,0,0,0,0,  2,0,0,0,0,0, 1,0,0,0,0,0;
                        0,1,0,0,0,0,  0,2,0,0,0,0, 0,1,0,0,0,0;
                        0,0,1,0,0,0,  0,0,2,0,0,0, 0,0,1,0,0,0;
                        0,0,0,1,0,0,  0,0,0,2,0,0, 0,0,0,1,0,0;
                        0,0,0,0,1,0,  0,0,0,0,2,0, 0,0,0,0,1,0;
                        0,0,0,0,0,1,  0,0,0,0,0,2, 0,0,0,0,0,1];
scal_mat_u = zeros((dataLen-2)*6,dataLen*6);
for i=1:dataLen-2
    scal_mat_u((i-1)*6+1:i*6,(i-1)*6+1:(i-1)*6+18)=scal_mat_u_cell;
end

tau= zeros(param.numJ,dataLen-2);
p=param;
while(optCont)
    % calculate the c value (is there contact)
     
    % in order to increase the accuracy, we over sampled the original
    % states to 10 times more. Since the problem we are solving here is
    % convex, adding more variables will not increase the computation too
    % much. 
    c.c_toe = zeros(1,dataLen-2);
    c.c_heel = zeros(1,dataLen-2);
    J.toe = cell(1,dataLen-2);
    J.heel = cell(1,dataLen-2);
    
    for i =1:dataLen-2
        q1 =prob.x0(1:p.numJ,i);
        q2 =prob.x0(1:p.numJ,i+1);
        q3 =prob.x0(1:p.numJ,i+2);
        
        
        
        % check toe pos
        
        
   

        
        [tau(:,i),J.toe{1,i},J.heel{1,i},c.c_toe(1,i),c.c_heel(1,i)] = u_no_ext(q1,q2,q3,param);
        
        
        
        
    end
    
    
    %% solve u,f based on q,dq,ddq
    
    tau = reshape(tau,[size(tau,1)*size(tau,2),1]);
    JT_toe_cell = cellfun(@transpose,J.toe,'UniformOutput',false);
    JT_heel_cell = cellfun(@transpose,J.heel,'UniformOutput',false);
    JT_toe_big = blkdiag(JT_toe_cell{:})*scal_mat;
    JT_heel_big = blkdiag(JT_heel_cell{:})*scal_mat;
    R = diag(repmat(param.jointW,1,dataLen-2));
    R = scal_mat_u.'*R*scal_mat_u;
    
    W_toe_cell = arrayfun(c_matFun,c.c_toe,'UniformOutput',false);
    W_toe = blkdiag(W_toe_cell{:});
    W_toe = scal_mat.'*W_toe*scal_mat;
    
    W_heel_cell = arrayfun(c_matFun,c.c_heel,'UniformOutput',false);
    W_heel = blkdiag(W_heel_cell{:});
    W_heel = scal_mat.'*W_heel*scal_mat;
    
    %the front ankle cannot create positive torque
%     u_con = repmat({[1,0,0,0,0,0]},[1,over_samp_len]);
%     u_con = blkdiag(u_con{:});
    

%     u_limit = 0;
%     try_time = 0;
%     while(~pass_err_check)
%         cvx_precision high
        cvx_begin
        variables f_toe(2*dataLen,1) f_heel(2*dataLen,1) u(param.numJ*dataLen,1)
        minimize( norm(JT_toe_big*f_toe+JT_heel_big*f_heel+scal_mat_u*u-tau,2)+u.'*R*u+f_toe.'*W_toe*f_toe+f_heel.'*W_heel*f_heel)    %
        subject to
        fri_mat*f_toe<=0;
        fri_mat*f_heel<=0;
%         u_con*u<=u_limit;
        
        cvx_end
        
%         if norm(JT_toe_big*f_toe+JT_heel_big*f_heel+u-tau)<700
%             pass_err_check = true;
%         else
%             u_limit = u_limit+7.5;
%             try_time=try_time+1;
%             if try_time==20
%                 error('exceed max cvx trial');
%                 
%             end
%         end
%     end
    

    
    
    
    % reshape to fit the format later
    f_toe = reshape(f_toe,[2,dataLen]);
    f_heel= reshape(f_heel,[2,dataLen]);
    u = reshape(u,[param.numJ,dataLen]);
    
    
    
    
    tau = reshape(tau,[param.numJ,dataLen-2]);

    
    
    
    %% with the previously got u, f_toe, f_heel, now we solve for the trajectory
    
    % since this problem will be non-convex and nonlinear, we use fmincon
    % to solve this problem
%     prob.nonlcon = @(x)gaitConst_inv(x,param);
    prob.objective = @(x)dynObj(x,param,u,f_toe,f_heel);
%     prob.objective = @(x)dynObj(x,param);
    
   
    
    
    [x,fval,exitflag,output] = fmincon(prob);
    
   


 
    

    res.x = x;
    res.u =u;
    res.f_toe = f_toe;
    res.f_heel = f_heel;
    save(['res_loop_',num2str(tot_count)],'res');
    
    tot_count = tot_count+1;
    if(tot_count==100)
        break;
    else
        prob.x0=x;
    end
    optCont = true;
%     param.dynCoeff = param.dynCoeff*1.03;
    iterTime = 1.003*iterTime;
    stepTol = 0.97*stepTol;
    options = optimoptions('fmincon','Algorithm','interior-point','MaxIter',floor(iterTime),'MaxFunctionEvaluations',floor(iterTime)*1000,...
     'Display','iter','TolCon',stepTol,'GradObj','on',...
    'StepTolerance',1e-5,'UseParallel',true,'DiffMinChange',0,'ScaleProblem',true,'ObjectiveLimit',1e-4);
    prob.options = options;
    disp(['current loop num: ',num2str(tot_count)]);
end



























