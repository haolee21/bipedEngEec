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
param.jointW = [10,1,1,1,1,1];

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
param.dynCoeff=100;
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
% q_ext = [[-q(6,end-1);-q(5,end-1);-pi-q(4,end-1);-pi-q(3,end-1);-q(2,end-1);-q(1,end-1)],q,[-q(6,2);-q(5,2);-pi-q(4,2);-pi-q(3,2);-q(2,2);-q(1,2)]];
% t_ext2 = -param.sampT:param.sampT:dataLen*param.sampT;
% t_samp = t_ext2+0.5*param.sampT;
% q_samp = interp1(t_ext2,q_ext.',t_samp(1:end-1),'spline').';
% dq =(q_samp(:,2:end)-q_samp(:,1:end-1))/param.sampT;
% dq_ext = [[-dq(6,end-1);-dq(5,end-1);-dq(4,end-1);-dq(3,end-1);-dq(2,end-1);-dq(1,end-1)],dq,[-dq(6,2);-dq(5,2);-dq(4,2);-dq(3,2);-dq(2,2);-dq(1,2)]];
% dq_samp = interp1(t_ext2,dq_ext.',t_samp(1:end-1),'spline').';
% ddq = (dq_samp(:,2:end)-dq_samp(:,1:end-1))/param.sampT;

dq = gradient(q)/param.sampT;
ddq = gradient(dq)/param.sampT;
f = zeros(4,dataLen);



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
    -45/180*pi*ones(1,dataLen);
    param.max_ank_vel*ones(1,dataLen);
    param.max_kne_vel*ones(1,dataLen);
    param.max_hip_vel*ones(1,dataLen);
    param.max_hip_vel*ones(1,dataLen);
    param.max_kne_vel*ones(1,dataLen);
    param.max_ank_vel*ones(1,dataLen)];

prob.lb = [45*ones(1,dataLen)/180*pi;
    -179/180*pi*ones(1,dataLen);
    -75/180*pi*ones(1,dataLen);
    -260/180*pi*ones(1,dataLen);
    -0*ones(1,dataLen)/180*pi;
    -135/180*pi*ones(1,dataLen);
    -param.max_ank_vel*ones(1,dataLen);
    -param.max_kne_vel*ones(1,dataLen);
    -param.max_hip_vel*ones(1,dataLen);
    -param.max_hip_vel*ones(1,dataLen);
    -param.max_kne_vel*ones(1,dataLen);
    -param.max_ank_vel*ones(1,dataLen)];

% conditions for start/end angles, velocity
Aeq = zeros(param.numJ*2,param.numJ*2*(dataLen));

Aeq(1:param.numJ,1:param.numJ*2)=[0,0,0,0,0,1, 0,0,0,0,0,0;   %start frame
                                  0,0,0,0,1,0, 0,0,0,0,0,0;
                                  0,0,0,1,0,0, 0,0,0,0,0,0;
                                  0,0,1,0,0,0, 0,0,0,0,0,0;
                                  0,1,0,0,0,0, 0,0,0,0,0,0;
                                  1,0,0,0,0,0, 0,0,0,0,0,0];
                                    
% endframe
Aeq(1:param.numJ,end-param.numJ*2+1:end) =  [ 1,0,0,0,0,0, 0,0,0,0,0,0;
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
                                                   
                                                   
% prob.Aeq = Aeq;
% prob.beq = [ 0;0;-pi;-pi;0;0;
%              0;0;0;0;0;0;-pi];    

param.Aeq = Aeq;
param.beq =  [ 0;0;-pi;-pi;0;0;
             0;0;0;0;0;0;-pi];    

   
%torso angle limit
Asamp = zeros(4,param.numJ*2); % create A for single frame
Asamp(1:2,1:3) = [-1,-1,-1;
                   1, 1, 1];
% torso speed limit
Asamp(3:4,1+param.numJ:3+param.numJ)=[-1,-1,-1;
                                       1, 1, 1];

Acell = repmat({Asamp},1,floor(dataLen));
prob.Aineq = blkdiag(Acell{:});
Bsamp = [-90/180*pi;
    110/180*pi;
    10/180*pi/param.sampT;
    10/180*pi/param.sampT];
prob.bineq = repmat(Bsamp,floor(param.gaitT/param.sampT+1),1);

iterTime=5;
stepTol = 1e-5;
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
cvx_solver sedumi
% we build friction matrix outside since it will never change
% constraints:
% fri_mat*f<=0
fri_mat = [1, -param.miu;
    -1, -param.miu];
fri_mat = repmat({fri_mat},1,dataLen);
fri_mat = blkdiag(fri_mat{:});
param.fri_mat = fri_mat;


tau= zeros(param.numJ,dataLen);

d1 = designfilt('lowpassiir','FilterOrder',12, ...
    'HalfPowerFrequency',10*param.sampT/10,'DesignMethod','butter');

while(optCont)
    % calculate the c value (is there contact)
     
    % in order to increase the accuracy, we over sampled the original
    % states to 10 times more. Since the problem we are solving here is
    % convex, adding more variables will not increase the computation too
    % much. 
    c.c_toe = param.k0/param.k1*ones(1,dataLen);
    c.c_heel = param.k0/param.k1*ones(1,dataLen);
    J.toe = cell(1,dataLen);
    J.heel = cell(1,dataLen);
    
%     q = interp1(t_ext,q.',t_over_samp,'spline').';
    
    % since there might be some jittering in the original trajectory, we
    % need to filter them
%     q = filtfilt(d1,q.').';
%     % however, filtered q may surpasses the upper/lower bound, we can set
%     % them to be max/min here
%     q(1,q(1,:)>prob.ub(1,1)) = prob.ub(1,1)-0.0017;
%     q(2,q(2,:)>prob.ub(2,1)) = prob.ub(2,1)-0.0017;
%     q(3,q(3,:)>prob.ub(3,1)) = prob.ub(3,1)-0.0017;
%     q(4,q(4,:)>prob.ub(4,1)) = prob.ub(4,1)-0.0017;
%     q(5,q(5,:)>prob.ub(5,1)) = prob.ub(5,1)-0.0017;
%     q(6,q(6,:)>prob.ub(6,1)) = prob.ub(6,1)-0.0017;
%     
%     q(1,q(1,:)<prob.lb(1,1)) = prob.lb(1,1)+0.0017;
%     q(2,q(2,:)<prob.lb(2,1)) = prob.lb(2,1)+0.0017;
%     q(3,q(3,:)<prob.lb(3,1)) = prob.lb(3,1)+0.0017;
%     q(4,q(4,:)<prob.lb(4,1)) = prob.lb(4,1)+0.0017;
%     q(5,q(5,:)<prob.lb(5,1)) = prob.lb(5,1)+0.0017;
%     q(6,q(6,:)<prob.lb(6,1)) = prob.lb(6,1)+0.0017;
%     
%     q = interp1(t_over_samp,q.',t_ext,'spline').';
%     q = interp1(t_ext,q.',t_over_samp,'spline').';
    
    
    
    
    dq = gradient(q)/param.sampT;
    ddq = gradient(dq)/param.sampT;
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
                c.c_toe(1,i) = param.k0/(cur_c^2*0.1+param.k1);
            end
            
        end
        if(heel_h<param.toe_th)
            cur_c = param.k*(param.toe_th-heel_h)^2-param.cmax*(param.toe_th-heel_h)/param.dmax*heel_vel(2);
            if(cur_c>0)
                c.c_heel(1,i)=param.k0/(cur_c^2*0.1+param.k1);
                
            end
        end
        
    end
    
    % however, the above tau(:,i) is often off, for here I will use
    % fmincon to apporach it
    
%     prob2.objective=@(u)inverseDyn(u,q,dq,param);
%     prob2.x0 = tau;
%     options2 = optimoptions('fmincon','Algorithm','interior-point','MaxIter',1000,'MaxFunctionEvaluations',1e6,...
%     'Display','iter','TolCon',1e-3,'GradObj','on',...
%     'StepTolerance',1e-8,'UseParallel',true,'DiffMinChange',0,'ScaleProblem',true,'ObjectiveLimit',1e-4);
%     prob2.options = options2;
%     prob2.solver = 'fmincon';
%     prob2.ub = Inf*ones(6,dataLen);
%     prob2.lb = -Inf*ones(6,dataLen);
%     [tau2,~,~,~] = fmincon(prob2);
    
    
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
    
    %the front ankle cannot create positive torque
%     u_con = repmat({[1,0,0,0,0,0]},[1,over_samp_len]);
%     u_con = blkdiag(u_con{:});
    

%     u_limit = 0;
%     try_time = 0;
%     while(~pass_err_check)
        cvx_precision high
        cvx_begin
        variables f_toe(2*dataLen,1) f_heel(2*dataLen,1) u(param.numJ*dataLen,1)
        minimize( norm(JT_toe_big*f_toe+JT_heel_big*f_heel+u-tau,2)+u.'*R*u+f_toe.'*W_toe*f_toe+f_heel.'*W_heel*f_heel)    %
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
    

    %down-sample the u, f_toe, f_heel
    
    
    
%     
% 
%     
%     %% simulate the trajectory with forward dynamic, with q(:,1), dq(:,1) as initial conditions
%     ics = [q(:,1);dq(:,1)];
%     opts = odeset('RelTol',1e-3,'AbsTol',1e-4);
%     [t_sim,x_sim] = ode45(@(t,x)robot_dyn(t,x,t_over_samp,u.',f_toe.',f_heel.',param),[0,param.gaitT],ics,opts);
%         
%     x_show = interp1(t_sim,x_sim,t_ext);
%     test = [x_show.';u;f_toe;f_heel];
    tau = reshape(tau,[param.numJ,dataLen]);
%     tau = interp1(t_over_samp,tau.',t_ext,'spline').';
%     u = interp1(t_over_samp,u.',t_ext,'spline').';
%     f_toe = interp1(t_over_samp,f_toe.',t_ext,'spline').';
%     f_heel = interp1(t_over_samp,f_heel.',t_ext,'spline').';
%     q = interp1(t_sim,x_sim(:,1:param.numJ),t_ext,'spline').';
%     dq = interp1(t_sim,x_sim(:,param.numJ+1:end),t_ext,'spline').';
%     q = interp1(t_over_samp,q.',t_ext).';
%     dq = interp1(t_over_samp,dq.',t_ext).';
    %% with the previously got u, f_toe, f_heel, now we solve for the trajectory
    
    % since this problem will be non-convex and nonlinear, we use fmincon
    % to solve this problem
%     prob.nonlcon = @(x)gaitConst_inv(x,param);
    prob.objective = @(x)dynObj(x,param,u,f_toe,f_heel);
%     prob.objective = @(x)dynObj(x,param);
    prob.x0 = [q;dq];
   
    
    
    [x,fval,exitflag,output] = fmincon(prob);
    
    q=x(1:param.numJ,:);
    dq = (gradient(q)/param.sampT);

 
    

    res.x = x;
    res.u =u;
    res.f_toe = f_toe;
    res.f_heel = f_heel;
    save(['res_loop_',num2str(tot_count)],'res');
    
    tot_count = tot_count+1;
    if(tot_count==10)
        break;
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



























