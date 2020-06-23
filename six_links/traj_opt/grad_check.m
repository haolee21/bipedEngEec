%% gradient checker

%baseline x

clear;
modelName='human_3';
warning on verbose
%add share functions
addpath dyn/
addpath obj/
addpath gaitCon/
addpath plotRobot/
dbstop if error


addpath (['../',modelName,'/robotGen/'])
addpath (['../',modelName,'/robotGen/grad/'])
addpath (['../',modelName,'/robotGen/posCons/'])
addpath (['../',modelName,'/robotGen/dyn/'])
addpath (['../',modelName,'/robotGen/obj/'])
addpath (['../',modelName,'/robotGen/grf/'])
addpath (['../',modelName,'/robotGen/knee_spring/'])

% load data, here I choose an infeasiable data set
data = load('grad_check_data').result;
model = load(['../',modelName,'/robotGen/model']).model;
%% simulation parameter
param = data.param;
x1 = data.x;

param.dmax =1e-3;
param.cmax=500;
param.k=0;
param.us=0.8;
param.ud=0.6;
param.joint_fri=0.003;
param.ank_stiff=408.65;
param.gaitLen = model.l_foot*3;
amp = max(abs(x1),[],2)*0.00001;%use 0.1% of the max value of each row as amp
x2 = x1+randn(size(x1,1),size(x1,2)).*amp;

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

% we check all x 
err_rate =0;
for i=1:size(x1,2)
    x1_row = x1(1:param.numJ*2,i).';
    x2_row = x2(1:param.numJ*2,i).';
    dx_row = x2_row-x1_row;
    
    V_val1 = six_V(x1_row(2),x1_row(3),x1_row(4),x1_row(5),x1_row(6),x1_row(7),x1_row(8),x1_row(9),x1_row(10),x1_row(11),x1_row(12));
    V_val2 = six_V(x2_row(2),x2_row(3),x2_row(4),x2_row(5),x2_row(6),x2_row(7),x2_row(8),x2_row(9),x2_row(10),x2_row(11),x2_row(12));
    diff_V = V_val2-V_val1;
    grad_V1 = dV_dx(x1_row(2),x1_row(3),x1_row(4),x1_row(5),x1_row(6),x1_row(7),x1_row(8),x1_row(9),x1_row(10),x1_row(11),x1_row(12));
    grad_V2 = dV_dx(x2_row(2),x2_row(3),x2_row(4),x2_row(5),x2_row(6),x2_row(7),x2_row(8),x2_row(9),x2_row(10),x2_row(11),x2_row(12));
    err = diff_V.' - 0.5*(grad_V1+grad_V2).'*dx_row.';
    err_rate_temp = norm(err)/norm(diff_V);
    if(err_rate_temp>err_rate)
        err_rate=err_rate_temp;
    end
    
end
disp(['dV_dx gradient err:',string(err_rate)]);

% check dG_dx
err_rate =0;
for i=1:size(x1,2)
    x1_row = x1(1:param.numJ,i).';
    x2_row = x2(1:param.numJ,i).';
    dx_row = x2_row-x1_row;
    
    G_val1 = six_G(x1_row(1),x1_row(2),x1_row(3),x1_row(4),x1_row(5),x1_row(6));
    G_val2 = six_G(x2_row(1),x2_row(2),x2_row(3),x2_row(4),x2_row(5),x1_row(6));
    diff_G = G_val2-G_val1;
    grad_G1 = dG_dx(x1_row(1),x1_row(2),x1_row(3),x1_row(4),x1_row(5),x1_row(6));
    grad_G2 = dG_dx(x2_row(1),x2_row(2),x2_row(3),x2_row(4),x2_row(5),x1_row(6));
    err = diff_G.' - 0.5*(grad_G1+grad_G2).'*dx_row.';
    
    err_rate_temp = norm(err)/norm(diff_G);
    if(err_rate_temp>err_rate)
        err_rate=err_rate_temp;
    end
end
disp(['dG_dx gradient err:',string(err_rate)]);
%% check Tau_Fext
err_rate =0;
err_sum =0;
count = 0;
for i=1:size(x1,2)
    
    x1_row = x1(:,i).';
    x2_row = x2(:,i).';
    dx_row = x2_row-x1_row;
    [tau_ext1,dTau_ext1,Fn1,Fs1,dFn1,dFs1] = toe_grf(x1_row,param);
    [tau_ext2,dTau_ext2,Fn2,Fs2,dFn2,dFs2] = toe_grf(x2_row,param);
    
    diff_Fn = Fn2-Fn1;
    err_Fn = diff_Fn - dx_row(1:12)*0.5*(dFn1+dFn2);
    if(norm(err_Fn)>1e-2)
        err_rate_Fn = norm(err_Fn)/norm(diff_Fn);
    end
    
    diff_Fs = Fs2-Fs1;
    err_Fs = diff_Fs - dx_row(1:12)*0.5*(dFs1+dFs2);
    if(norm(err_Fs)>0)
        err_rate_Fs = norm(err_Fs)/norm(diff_Fs);
    end
    
    
    
    diff = tau_ext2-tau_ext1;
    err = diff - dx_row*0.5*(dTau_ext1+dTau_ext2);
    err_rate_temp = norm(err)/norm(diff);
    if(~isnan(err_rate_temp))
        err_sum = err_sum+err_rate_temp;
        count = count+1;
    end
    if(err_rate_temp>err_rate)
        err_rate = err_rate_temp;
    end
    
end
disp(['dToe_ext gradient err:',string(err_rate)]);
disp(['dToe_ext avg gradient err:',string(err_sum/count)]);

err_rate =0;
err_sum =0;
count = 0;
for i=floor(data.param.gaitT/data.param.sampT/2)-1:size(x1,2)
    
    x1_row = x1(:,i).';
    x2_row = x2(:,i).';
    dx_row = x2_row-x1_row;
    [tau_ext1,dTau_ext1] = heel_grf(x1_row,param);
    [tau_ext2,dTau_ext2] = heel_grf(x2_row,param);
    
    diff = tau_ext2-tau_ext1;
    err = abs(diff - dx_row*0.5*(dTau_ext1+dTau_ext2));
    err_rate_temp = norm(err)/norm(diff);
    if(~isnan(err_rate_temp))
        err_sum = err_sum+err_rate_temp;
        count = count +1;
    end
    if(err_rate_temp>err_rate)
        err_rate = err_rate_temp;
    end
    
end
disp(['dHeel_ext gradient err:',string(err_rate)]);
disp(['avg dHeel_ext gradient err:', string(err_sum/count)]);
% 

%% check f_x

% first we treat f_x as one input vector function
% err_rate =0;
% for i=1:size(x1,2)
%     x1_row = x1(:,i).';
%     x2_row = x2(:,i).';
%     dx_row = x2_row-x1_row;
%     [out_t1,grad_t1,grad_t_k1] = f_x(x1_row,param,1);
%     [out_t2,grad_t2,grad_t_k2] = f_x(x2_row,param,1);
%     
%     err = out_t2-out_t1-dx_row*0.5*(grad_t1+grad_t2);
%     err_rate_temp = norm(err)/norm(out_t2-out_t1);
% %     if(err_rate_temp>err_rate)
% %         err_rate=err_rate_temp;
% %     end
%     err_rate = err_rate + err_rate_temp;
% end
% disp(['f_x (single var) avg gradient err:',string(err_rate/size(x1,2))]);


err_rate =0;
for i=1:size(x1,2)
    x1_row = x1(:,i);
    x2_row = x2(:,i);
    dx_row = x2_row-x1_row;
    [f1,df1s,df1u] = f_x2(x1_row,param);
    [f2,df2s,df2u] = f_x2(x2_row,param);
    
    err = f2-f1-0.5*([df1s;df1u]+[df2s;df2u]).'*dx_row;
    err_rate_temp = norm(err)/norm(f2-f1);
%     if(err_rate_temp>err_rate)
%         err_rate=err_rate_temp;
%     end
    err_rate = err_rate + err_rate_temp;
end
disp(['f_x2 (single var) avg gradient err:',string(err_rate/size(x1,2))]);




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
% [c1,grad1] = dynConst(x1,param);
% [c2,grad2]=dynConst(x2,param);
% dx = reshape(x2-x1,[size(x1,1)*size(x1,2),1]);
% err = c2-c1-0.5*(grad1+grad2).'*dx;
% 
% 
% disp(['dyn gradient err:',string(gather(norm(err)/norm(c2-c1)))]);

%check on u and fext
x2(1:param.numJ*2,:) = x1(1:param.numJ*2,:);

% check only position
% x2(param.numJ+1:end,:)=x1(param.numJ+1:end,:);

% check only velocity 
% x2(2*param.numJ+1:end,:)=x1(2*param.numJ+1:end,:);
% x2(1:param.numJ,:)=x1(1:param.numJ,:);

% %check on states
% x2(param.numJ*2+1:end,:) = x1(param.numJ*2+1:end,:);

[c1,grad1] = dynConst2(x1,param);
[c2,grad2]=dynConst2(x2,param);
dx = reshape(x2-x1,[size(x1,1)*size(x1,2),1]);
err = c2-c1-0.5*(grad1+grad2).'*dx;


disp(['dyn2 gradient err:',string(gather(norm(err)/norm(c2-c1)))]);
% clear c1 grad1 c2 grad2

%% check nonlinear constraint
[c1,ceq1,gradc1,gradceq1] = five_link_nonlcon(x1,param);
[c2,ceq2,gradc2,gradceq2] = five_link_nonlcon(x2,param);

% err_c = c2-c1-0.5*(gradc1+gradc2).'*dx;
% disp(['c gradient err:',string(gather(norm(err_c)/norm(c2-c1)))]);

err_ceq = ceq2-ceq1 - 0.5*(gradceq1+gradceq2).'*dx;
disp(['ceq gradient err:',string(gather(norm(err_ceq)/norm(ceq2-ceq1)))]);


%% check object function 
[dObj1,dObjGrad1]=objFun(x1,param);
[dObj2,dObjGrad2]=objFun(x2,param);

err = dObj2-dObj1 -sum(0.5*(dObjGrad1+dObjGrad2).*(x2-x1),'all');
disp(['obj gradient err:',string(gather(err/(dObj2-dObj1)))]);