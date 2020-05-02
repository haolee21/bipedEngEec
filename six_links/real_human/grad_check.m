%% gradient checker

%baseline x

clear;
clc;
addpath dyn/
addpath robotGen/
addpath robotGen/grad
addpath obj/
addpath gaitCon/
%% simulation parameter

param.numJ=5;
param.toe_th =1e-2;

param.gaitT = 0.08;
param.sampT = 0.002;

param.gaitLen = 1;

time = 0:param.sampT:param.gaitT;

% set torque/angular velocity constraints
max_tau = 30;
max_vel = 30/180*pi;
%% initialize joint pos and torque
qmax = 170/180/pi;
% q = qmax*sin((2*time/param.gaitT+randn(param.jointNum,1))*pi);
% dq = qmax*sin((2*time/param.gaitT+randn(param.jointNum,1))*pi)*10;


% add some noise to the states
q = [pi/2*ones(1,length(time))+randn(1,length(time))*0.01*pi/2;
     -pi/2*ones(1,length(time))+randn(1,length(time))*0.01*pi/2;
     zeros(1,length(time))+randn(1,length(time))*0.01*pi/2;
     -pi*ones(1,length(time))+randn(1,length(time))*0.01*pi/2;
     pi/2*ones(1,length(time))+randn(1,length(time))*0.01*pi/2];
dq = zeros(param.numJ,length(time))+randn(param.numJ,length(time))*0.1*pi/2;


u = zeros(param.numJ,length(q))+randn(param.numJ,length(q))*0.1*pi/2;

ext_tau = zeros(size(time,2),param.numJ);


x1 = [q;dq;u];

x2 = [q+randn(size(q,1),size(q,2))*0.001; dq+randn(size(q,1),size(q,2))*0.01;u+randn(size(q,1),size(q,2))*0.01];


dx = reshape(x2-x1,[size(x1,1)*size(x1,2),1]);


%% check gaitLenCon
[c1,grad1] = gaitLenCon(x1,param);
[c2,grad2]=gaitLenCon(x2,param);

err = c2-c1-0.5*(grad1+grad2).'*dx;

disp(['gait len gradient err:',string(norm(err)/norm(c2-c1))]);

clear c1 grad1 c2 grad2

%% check yposCon

[c1,grad1]=yposCon(x1);
[c2,grad2]=yposCon(x2);
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

V_val1 = five_V(x1_row(2),x1_row(3),x1_row(4),x1_row(5),x1_row(6),x1_row(7),x1_row(8),x1_row(9),x1_row(10))*x1_row(6:10).';
V_val2 = five_V(x2_row(2),x2_row(3),x2_row(4),x2_row(5),x2_row(6),x2_row(7),x2_row(8),x2_row(9),x2_row(10))*x2_row(6:10).';
diff_V = V_val2-V_val1;
grad_V1 = dV_dx(x1_row(6),x1_row(7),x1_row(8),x1_row(9),x1_row(10),x1_row(2),x1_row(3),x1_row(4),x1_row(5));
grad_V2 = dV_dx(x2_row(6),x2_row(7),x2_row(8),x2_row(9),x2_row(10),x2_row(2),x2_row(3),x2_row(4),x2_row(5));
err = diff_V - 0.5*(grad_V1+grad_V2).'*dx_row.';
disp(['dV_dx gradient err:',string(gather(norm(err)/norm(diff_V)))]);

% check dG_dx
G_val1 = five_G(x1_row(1),x1_row(2),x1_row(3),x1_row(4),x1_row(5));
G_val2 = five_G(x2_row(1),x2_row(2),x2_row(3),x2_row(4),x2_row(5));
diff_G = G_val2-G_val1;
grad_G1 = dG_dx(x1_row(1),x1_row(2),x1_row(3),x1_row(4),x1_row(5));
grad_G2 = dG_dx(x2_row(1),x2_row(2),x2_row(3),x2_row(4),x2_row(5));
err = diff_G.' - 0.5*(grad_G1+grad_G2)*dx_row.';
disp(['dG_dx gradient err:',string(gather(norm(err)/norm(diff_G)))]);

% M is a tensor, it is a bit hard to check, we check the overall betaFun
% directly

% [c1,grad1]=betaFun(x1_row);
% [c2,grad2]=betaFun(x2_row);
% 
% err = (c2-c1).'-0.5*(grad1+grad2).'*(x2(:,30)-x1(:,30));
% disp(['betaFun gradient err:',string(gather(norm(err)/norm(c2-c1)))]);
% clear c1 grad1 c2 grad2



%% check f_x

% first we treat f_x as one input vector function

[out_t1,grad_t1,grad_t_k1] = f_x(x1_row,param.toe_th);
[out_t2,grad_t2,grad_t_k2] = f_x(x2_row,param.toe_th);

err = out_t2-out_t1-dx_row*0.5*(grad_t1+grad_t2);
disp(['f_x (single var) gradient err:',string(gather(norm(err)/norm(x2_row-x1_row)))]);

% test f_x as two input vectos function
x1_row_2 = x1(:,31).';
x2_row_2 = x2(:,31).';
dx_row_2 = x2_row_2-x1_row_2;
[out_t1,grad_t1,grad_t_k1] = f_x([x1_row,x1_row_2],param.toe_th);
[out_t2,grad_t2,grad_t_k2] = f_x([x2_row,x2_row_2],param.toe_th);

err = out_t2-out_t1-dx_row*0.5*(grad_t1+grad_t2)-dx_row_2*0.5*(grad_t_k1+grad_t_k2);
disp(['f_x (double var) gradient err:',string(gather(norm(err)/(norm(x2_row-x1_row)+norm(x2_row_2-x1_row_2))))]);

%% check dynCon
[c1,grad1] = dynConst(x1,param);
[c2,grad2]=dynConst(x2,param);

err = c2-c1-0.5*(grad1+grad2).'*dx;


disp(['dyn gradient err:',string(gather(norm(err)/norm(c2-c1)))]);
% clear c1 grad1 c2 grad2

