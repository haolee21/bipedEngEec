%% gradient checker

%baseline x

clear;
clc;
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
addpath robotGen/knee_spring/
%% simulation parameter

param.numJ=6;
numJ = param.numJ;
param.numJ=6;
param.toe_th = 5e-4;
param.head_h = 1.1 ; %the head should be at least 1.6m

param.gaitT = 0.5;
param.sampT = 0.001;
param.init_y = 1e-3; %initial feet height
param.gaitLen = 1.8;
param.hipLen=0.67;
param.jointW=[1,2,1,1,1,1];
param.gndclear = 5e-2;
param.fri_coeff=5;
param.floor_stiff=0.4;
param.ank_stiff=70;
time = 0:param.sampT:param.gaitT;
param.knee_stiff=500;
% set torque/angular velocity constraints
max_tau = 30;
max_vel = 30/180*pi;
%% initialize joint pos and torque
qmax = 170/180/pi;
% q = qmax*sin((2*time/param.gaitT+randn(param.jointNum,1))*pi);
% dq = qmax*sin((2*time/param.gaitT+randn(param.jointNum,1))*pi)*10;


% add some noise to the states
q = [pi/2*ones(1,length(time))+randn(1,length(time))*0.0001*pi/2;
     0/2*ones(1,length(time))+randn(1,length(time))*0.0001*pi/2;
     zeros(1,length(time))+randn(1,length(time))*0.0001*pi/2;
     -pi*ones(1,length(time))+randn(1,length(time))*0.0001*pi/2;
     pi/2*ones(1,length(time))+randn(1,length(time))*0.0001*pi/2;
     pi/2*ones(1,length(time))+randn(1,length(time))*0.0001*pi/2];
dq = zeros(param.numJ,length(time))+randn(param.numJ,length(time))*0.01*pi/2;


u = zeros(param.numJ,length(q))+randn(param.numJ,length(q))*0.01*pi/2;

ext_tau = zeros(size(time,2),param.numJ);


x1 = [q;dq;u];

x2 = [q+randn(size(q,1),size(q,2))*0.000001; dq+randn(size(q,1),size(q,2))*0.000001;u+randn(size(q,1),size(q,2))*0.000001];


dx = reshape(x2-x1,[size(x1,1)*size(x1,2),1]);


%% check gaitLenCon
[c1,grad1] = hipCon(x1,param);
[c2,grad2]=hipCon(x2,param);

err = c2-c1-0.5*(grad1+grad2).'*dx;

disp(['hip len gradient err:',string(norm(err)/norm(c2-c1))]);

clear c1 grad1 c2 grad2

%% check yposCon

[c1,grad1]=yposCon(x1,param);
[c2,grad2]=yposCon(x2,param);
err = c2-c1-0.5*(grad1+grad2).'*dx;


disp(['yposCon gradient err:',string(gather(norm(err)/norm(c2-c1)))]);
clear c1 grad1 c2 grad2

%% check betaFun
% there are multiple functions that return gradient in betaFun, 
% dV_dx, dG_dx,

% these functions takes row vector inputs

% check dV_dx
x1_row = x1(:,30).';
x2_row = x2(:,30).';
dx_row = x2_row-x1_row;

V_val1 = six_V(x1_row(2),x1_row(3),x1_row(4),x1_row(5),x1_row(6),x1_row(7),x1_row(8),x1_row(9),x1_row(10),x1_row(11),x1_row(12));
V_val2 = six_V(x2_row(2),x2_row(3),x2_row(4),x2_row(5),x2_row(6),x2_row(7),x2_row(8),x2_row(9),x2_row(10),x2_row(11),x2_row(12));
diff_V = V_val2-V_val1;
grad_V1 = dV_dx(x1_row(7),x1_row(8),x1_row(9),x1_row(10),x1_row(11),x1_row(12),x1_row(2),x1_row(3),x1_row(4),x1_row(5),x1_row(6));
grad_V2 = dV_dx(x2_row(7),x2_row(8),x2_row(9),x2_row(10),x2_row(11),x2_row(12),x2_row(2),x2_row(3),x2_row(4),x2_row(5),x2_row(6));
err = diff_V.' - 0.5*(grad_V1+grad_V2).'*dx_row.';
disp(['dV_dx gradient err:',string(gather(norm(err)/norm(diff_V)))]);

% check dG_dx
G_val1 = six_G(x1_row(1),x1_row(2),x1_row(3),x1_row(4),x1_row(5),x1_row(6));
G_val2 = six_G(x2_row(1),x2_row(2),x2_row(3),x2_row(4),x2_row(5),x1_row(6));
diff_G = G_val2-G_val1;
grad_G1 = dG_dx(x1_row(1),x1_row(2),x1_row(3),x1_row(4),x1_row(5),x1_row(6));
grad_G2 = dG_dx(x2_row(1),x2_row(2),x2_row(3),x2_row(4),x2_row(5),x1_row(6));
err = diff_G.' - 0.5*(grad_G1+grad_G2).'*dx_row.';
disp(['dG_dx gradient err:',string(gather(norm(err)/norm(diff_G)))]);





% M is a tensor, it is a bit hard to check, we check the overall betaFun
% directly

% [c1,grad1]=betaFun(x1_row);
% [c2,grad2]=betaFun(x2_row);
% 
% err = (c2-c1).'-0.5*(grad1+grad2).'*(x2(:,30)-x1(:,30));
% disp(['betaFun gradient err:',string(gather(norm(err)/norm(c2-c1)))]);
% clear c1 grad1 c2 grad2

%% check beta_out (function for ground touching)
% beta1 = beta_grf(x1_row(1),x1_row(2),x1_row(3),x1_row(4),x1_row(5),x1_row(6),param.toe_th);
% beta2 = beta_grf(x2_row(1),x2_row(2),x2_row(3),x2_row(4),x2_row(5),x2_row(6),param.toe_th);
% dBeta1_x1 = dbeta_dx1(x1_row(1),x1_row(2),x1_row(3),x1_row(4),x1_row(5),x1_row(6),param.toe_th);
% dBeta1_x2 = dbeta_dx2(x1_row(1),x1_row(2),x1_row(3),x1_row(4),x1_row(5),x1_row(6),param.toe_th);
% dBeta1_x3 = dbeta_dx3(x1_row(1),x1_row(2),x1_row(3),x1_row(4),x1_row(5),x1_row(6),param.toe_th);
% dBeta1_x4 = dbeta_dx4(x1_row(1),x1_row(2),x1_row(3),x1_row(4),x1_row(5),x1_row(6),param.toe_th);
% dBeta1_x5 = dbeta_dx5(x1_row(1),x1_row(2),x1_row(3),x1_row(4),x1_row(5),x1_row(6),param.toe_th);
% dBeta1_x6 = dbeta_dx6(x1_row(1),x1_row(2),x1_row(3),x1_row(4),x1_row(5),x1_row(6),param.toe_th);
% 
% dBeta2_x1 = dbeta_dx1(x2_row(1),x2_row(2),x2_row(3),x2_row(4),x2_row(5),x2_row(6),param.toe_th);
% dBeta2_x2 = dbeta_dx2(x2_row(1),x2_row(2),x2_row(3),x2_row(4),x2_row(5),x2_row(6),param.toe_th);
% dBeta2_x3 = dbeta_dx3(x2_row(1),x2_row(2),x2_row(3),x2_row(4),x2_row(5),x2_row(6),param.toe_th);
% dBeta2_x4 = dbeta_dx4(x2_row(1),x2_row(2),x2_row(3),x2_row(4),x2_row(5),x2_row(6),param.toe_th);
% dBeta2_x5 = dbeta_dx5(x2_row(1),x2_row(2),x2_row(3),x2_row(4),x2_row(5),x2_row(6),param.toe_th);
% dBeta2_x6 = dbeta_dx6(x2_row(1),x2_row(2),x2_row(3),x2_row(4),x2_row(5),x2_row(6),param.toe_th);
% 
% dBeta_dx =zeros(param.numJ,param.numJ);
% dx = x2_row(1:param.numJ)-x1_row(1:param.numJ);
% dBeta_dx=dBeta_dx+dx(1)*0.5*(dBeta1_x1+dBeta2_x1);
% dBeta_dx=dBeta_dx+dx(2)*0.5*(dBeta1_x2+dBeta2_x2);
% dBeta_dx=dBeta_dx+dx(3)*0.5*(dBeta1_x3+dBeta2_x3);
% dBeta_dx=dBeta_dx+dx(4)*0.5*(dBeta1_x4+dBeta2_x4);
% dBeta_dx=dBeta_dx+dx(5)*0.5*(dBeta1_x5+dBeta2_x5);
% dBeta_dx=dBeta_dx+dx(6)*0.5*(dBeta1_x6+dBeta2_x6);
% err=beta2-beta1-dBeta_dx;
% disp(['beta gradient err:',string(gather(norm(err)/norm(beta2-beta1)))]);
% 
% %% check beta*(u-G)
% u1 = x1_row(numJ*2+1:numJ*3);
% u2 = x2_row(numJ*2+1:numJ*3);
% G1 = six_G(x1_row(1),x1_row(2),x1_row(3),x1_row(4),x1_row(5),x1_row(6));
% G2 = six_G(x2_row(1),x2_row(2),x2_row(3),x2_row(4),x2_row(5),x2_row(6));
% out1 = beta1*(-G1).';
% out2 = beta2*(-G2).';
% dBeta1_dx =zeros(3*param.numJ,param.numJ);
% dBeta1_dx(1,:) = (-G1)*dbeta_dx1(x1_row(1),x1_row(2),x1_row(3),x1_row(4),x1_row(5),x1_row(6),param.toe_th).';
% dBeta1_dx(2,:) = (-G1)*dbeta_dx2(x1_row(1),x1_row(2),x1_row(3),x1_row(4),x1_row(5),x1_row(6),param.toe_th).';
% dBeta1_dx(3,:) = (-G1)*dbeta_dx3(x1_row(1),x1_row(2),x1_row(3),x1_row(4),x1_row(5),x1_row(6),param.toe_th).';
% dBeta1_dx(4,:) = (-G1)*dbeta_dx4(x1_row(1),x1_row(2),x1_row(3),x1_row(4),x1_row(5),x1_row(6),param.toe_th).';
% dBeta1_dx(5,:) = (-G1)*dbeta_dx5(x1_row(1),x1_row(2),x1_row(3),x1_row(4),x1_row(5),x1_row(6),param.toe_th).';
% dGdx1 = dG_dx(x1_row(1),x1_row(2),x1_row(3),x1_row(4),x1_row(5),x1_row(6));
% dTaudx = [zeros(numJ);zeros(numJ);eye(numJ)];
% grad1 = dBeta1_dx+dTaudx*beta1.'-dGdx1*beta1.';
% 
% dBeta2_dx =zeros(3*param.numJ,param.numJ);
% dBeta2_dx(1,:) = (-G2)*dbeta_dx1(x2_row(1),x2_row(2),x2_row(3),x2_row(4),x2_row(5),x2_row(6),param.toe_th).';
% dBeta2_dx(2,:) = (-G2)*dbeta_dx2(x2_row(1),x2_row(2),x2_row(3),x2_row(4),x2_row(5),x2_row(6),param.toe_th).';
% dBeta2_dx(3,:) = (-G2)*dbeta_dx3(x2_row(1),x2_row(2),x2_row(3),x2_row(4),x2_row(5),x2_row(6),param.toe_th).';
% dBeta2_dx(4,:) = (-G2)*dbeta_dx4(x2_row(1),x2_row(2),x2_row(3),x2_row(4),x2_row(5),x2_row(6),param.toe_th).';
% dBeta2_dx(5,:) = (-G2)*dbeta_dx5(x2_row(1),x2_row(2),x2_row(3),x2_row(4),x2_row(5),x2_row(6),param.toe_th).';
% dGdx2 = dG_dx(x2_row(1),x2_row(2),x2_row(3),x2_row(4),x2_row(5),x2_row(6));
% grad2 = dBeta1_dx-dGdx2*beta2.';
% 
% dx = x2_row-x1_row;
% err = out2-out1-0.5*(grad1+grad2).'*dx.';
% disp(['beta(u-G) gradient err:',string(gather(norm(err)/norm(out2-out1)))]);

%% check f_x

% first we treat f_x as one input vector function

[out_t1,grad_t1,grad_t_k1] = f_x(x1_row,param,850);
[out_t2,grad_t2,grad_t_k2] = f_x(x2_row,param,850);

err = out_t2-out_t1-dx_row*0.5*(grad_t1+grad_t2);
disp(['f_x (single var) gradient err:',string(gather(norm(err)/norm(out_t2-out_t1)))]);

% test f_x as two input vectos function
% x1_row_2 = x1(:,31).';
% x2_row_2 = x2(:,31).';
% dx_row_2 = x2_row_2-x1_row_2;
% [out_t1,grad_t1,grad_t_k1] = f_x([x1_row,x1_row_2],param);
% [out_t2,grad_t2,grad_t_k2] = f_x([x2_row,x2_row_2],param);
% 
% err = out_t2-out_t1-dx_row*0.5*(grad_t1+grad_t2);%-dx_row_2*0.5*(grad_t_k1+grad_t_k2);
% disp(['f_x (double var) gradient err:',string(gather(norm(err)/(norm(x2_row-x1_row)+norm(x2_row_2-x1_row_2))))]);

%% check dynCon
[c1,grad1] = dynConst(x1,param);
[c2,grad2]=dynConst(x2,param);
dx = reshape(x2-x1,[size(x1,1)*size(x1,2),1]);
err = c2-c1-0.5*(grad1+grad2).'*dx;


disp(['dyn gradient err:',string(gather(norm(err)/norm(c2-c1)))]);
% clear c1 grad1 c2 grad2

%% check nonlinear constraint
[c1,ceq1,gradc1,gradceq1] = five_link_nonlcon(x1,param);
[c2,ceq2,gradc2,gradceq2] = five_link_nonlcon(x2,param);

err_c = c2-c1-0.5*(gradc1+gradc2).'*dx;
disp(['c gradient err:',string(gather(norm(err_c)/norm(c2-c1)))]);

err_ceq = ceq2-ceq1 - 0.5*(gradceq1+gradceq2).'*dx;
disp(['ceq gradient err:',string(gather(norm(err_ceq)/norm(ceq2-ceq1)))]);


%% check object function 
[dObj1,dObjGrad1]=objFun(x1,param);
[dObj2,dObjGrad2]=objFun(x2,param);

err = dObj2-dObj1 -sum(0.5*(dObjGrad1+dObjGrad2).*(x2-x1),'all');
disp(['obj gradient err:',string(gather(err/(dObj2-dObj1)))]);