% Gradient check for contact invarient problems

clear;
clc;
addpath 'obj';
%load test data
u = load('u_contact_test.mat').u;
x1 = load('x_contact_test.mat').x;
f_toe = load('f_toe_test.mat').f_toe;
f_heel = load('f_heel_test.mat').f_heel;
param = load('param_contact_test.mat').param;
x2 = x1+randn(size(x1,1),size(x1,2))*0.001;

param.smooth_vel = 0.1;
param.initPos_w = 100;
param.init_vel_w = 1;
param.end_pos_w = 10;
%check f_dyn
err_rate =0;
for i=1:size(x1,2)
    x1_cur = x1(:,i);
    x2_cur = x2(:,i);
    
    [f1,df1x]=f_dyn(x1_cur,param,f_toe(:,i),f_heel(:,i),u(:,i));
    [f2,df2x]=f_dyn(x2_cur,param,f_toe(:,i),f_heel(:,i),u(:,i));
    dx = x2_cur-x1_cur;
    err = f2-f1-0.5*(df1x+df2x).'*dx;
    err_rate_temp = norm(err)/norm(f2-f1);
    if(err_rate_temp>err_rate)
        err_rate = err_rate_temp;
    end
    
end
disp([' max f_dyn error: ',num2str(err_rate)]);


%check dynObj
[out1,grad1]=dynObj(x1,param,u,f_toe,f_heel);
[out2,grad2]=dynObj(x2,param,u,f_toe,f_heel);
dx = x2-x1;
err = out2-out1-0.5*(grad1+grad2).'*reshape(dx,[size(dx,1)*size(dx,2),1]);

disp(['dynObj error',num2str(err/(out2-out1))]);