%% This will calculate the hessian of the object function and nonlinear constraint
clear;
modelName='human_3';
warning on verbose
%add share functions
addpath dyn/
addpath obj/
addpath gaitCon/
addpath plotRobot/
dbstop if error
addpath hessian

addpath (['../',modelName,'/robotGen/'])
addpath (['../',modelName,'/robotGen/grad/'])
addpath (['../',modelName,'/robotGen/posCons/'])
addpath (['../',modelName,'/robotGen/dyn/'])
addpath (['../',modelName,'/robotGen/obj/'])
addpath (['../',modelName,'/robotGen/grf/'])
addpath (['../',modelName,'/robotGen/knee_spring/'])

syms numJ toe_th sampT real
syms dmax cmax k us ud ank_stiff knee_stiff joint_fri hipLen real

param.numJ=6;
param.toe_th =toe_th;

param.sampT = sampT;

param.dmax =dmax;
param.cmax=cmax;
param.k=k;
param.us=us;
param.ud=ud;
param.joint_fri = joint_fri;
param.knee_stiff = knee_stiff;
param.ank_stiff= ank_stiff;
param.hipLen=hipLen;
%% hipVelCon

syms q1 q2 q3 q4 q5 q6 qd1 qd2 qd3 qd4 qd5 qd6 u1 u2 u3 u4 u5 u6 fx_toe fy_toe fx_heel fy_heel real
xmat = [q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6,u1,u2,u3,u4,u5,u6];
xmat2 = [q1, q2, q3, q4, q5, q6, qd1, qd2, qd3, qd4, qd5, qd6, u1, u2, u3, u4, u5, u6, fx_toe, fy_toe, fx_heel, fy_heel];
[~,gradc1] = hipVelCon(xmat.',param);
[~,~,gradc_grf,gradceq_grf] = grf_cons(xmat2.',param);
task_i = 1;
hess_hipVelCon = sym(zeros(3*param.numJ+4,3*param.numJ+4));
hess_grf_ceq_toe = sym(zeros(3*param.numJ+4,3*param.numJ+4));
hess_grf_c_toe = sym(zeros(3*param.numJ+4,3*param.numJ+4));
hess_grf_ceq_heel = sym(zeros(3*param.numJ+4,3*param.numJ+4));
hess_grf_c_heel = sym(zeros(3*param.numJ+4,3*param.numJ+4));
for i1=1:length(gradc1)
    for i2=i1:length(gradc1)% we only solve the upper triangle part to save resources, will generate the other half when use it.                         
        hess_hipVelCon(i1,i2) = diff(gradc1(i1,1),xmat(i2));
    end
end

for i1=1:length(gradc_grf)
    for i2=i1:length(gradc_grf) %this also only solve the upper triangle part
%         hess_grf_c_toe(i1,i2)=diff(gradc_grf(i1,1),xmat2(i2));
%         hess_grf_ceq_toe(i1,i2)=diff(gradceq_grf(i1,1),xmat2(i2));
%         hess_grf_c_heel(i1,i2)=diff(gradc_grf(i1,2),xmat2(i2));
        task{1,task_i}=@()matlabFunction(diff(gradceq_grf(i1,2),xmat2(i2)),'file',['hessian/hess_grf_ceq_heel',num2str(i1),num2str(i2)],'vars',{xmat2,toe_th,dmax, cmax, k, us, ud}); task_i=task_i+1;
%         hess_grf_ceq_heel(i1,i2)=diff(gradceq_grf(i1,2),xmat2(i2));
    end
end

task{1,task_i}=@()matlabFunction(hess_hipVelCon,'file','hessian/h_hipVelCon','vars',{xmat2,toe_th,}); task_i=task_i+1;
task{1,task_i}=@()matlabFunction(hess_grf_c_toe,'file','hessian/hess_grf_c_toe','vars',{xmat2,toe_th,dmax, cmax, k, us, ud}); task_i=task_i+1;
task{1,task_i}=@()matlabFunction(hess_grf_ceq_toe,'file','hessian/hess_grf_ceq_toe','vars',{xmat2,toe_th,dmax, cmax, k, us, ud}); task_i=task_i+1;
task{1,task_i}=@()matlabFunction(hess_grf_c_heel,'file','hessian/hess_grf_c_heel','vars',{xmat2,toe_th,dmax, cmax, k, us, ud}); task_i=task_i+1;
task{1,task_i}=@()matlabFunction(hess_grf_ceq_heel,'file','hessian/hess_grf_ceq_heel','vars',{xmat2,toe_th,dmax, cmax, k, us, ud}); task_i=task_i+1;
%% dynConst
% since it involves M^-1 in the equations, the previous method is not
% applicable, we need to solve it with f_x2 in hessian_fcn

% for here we need the tools such as dV_dxx, dG_dxx, dM_dxx

dVdx = dV_dx(q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6);
dGdx = dG_dx(q1,q2,q3,q4,q5,q6);
xmat = [q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6];
dV_dxx = sym(zeros(2*param.numJ,2*param.numJ,param.numJ));
dG_dxx = sym(zeros(param.numJ,param.numJ,param.numJ));
dM_dxx = sym(zeros(param.numJ,param.numJ,param.numJ,param.numJ));
dJ_toe_dxx = sym(zeros(param.numJ,param.numJ,6,param.numJ));
dJ_heel_dxx = sym(zeros(param.numJ,param.numJ,6,param.numJ));
for i=1:param.numJ*2
    dV_dxx(i,:,:)=diff(dVdx,xmat(i));
    if(i<=param.numJ)
        dG_dxx(i,:,:)=diff(dGdx,xmat(i));
        dM_dxx(2,i,:,:)=diff(dM_dx2(q2,q3,q4,q5,q6),xmat(i));
        dM_dxx(3,i,:,:)=diff(dM_dx3(q2,q3,q4,q5,q6),xmat(i));
        dM_dxx(4,i,:,:)=diff(dM_dx4(q2,q3,q4,q5,q6),xmat(i));
        dM_dxx(5,i,:,:)=diff(dM_dx5(q2,q3,q4,q5,q6),xmat(i));
        dM_dxx(6,i,:,:)=diff(dM_dx6(q2,q3,q4,q5,q6),xmat(i));
        dJ_toe_dxx(1,i,:,:)=diff(dJ_q1([q1,q2,q3,q4,q5,q6]),xmat(i));
        dJ_toe_dxx(2,i,:,:)=diff(dJ_q2([q1,q2,q3,q4,q5,q6]),xmat(i));
        dJ_toe_dxx(3,i,:,:)=diff(dJ_q3([q1,q2,q3,q4,q5,q6]),xmat(i));
        dJ_toe_dxx(4,i,:,:)=diff(dJ_q4([q1,q2,q3,q4,q5,q6]),xmat(i));
        dJ_toe_dxx(5,i,:,:)=diff(dJ_q5([q1,q2,q3,q4,q5,q6]),xmat(i));
        dJ_toe_dxx(6,i,:,:)=diff(dJ_q6([q1,q2,q3,q4,q5,q6]),xmat(i));
        
        dJ_heel_dxx(1,i,:,:) = diff(dJ2_q1([q1,q2,q3,q4,q5,q6]),xmat(i));
        dJ_heel_dxx(2,i,:,:) = diff(dJ2_q2([q1,q2,q3,q4,q5,q6]),xmat(i));
        dJ_heel_dxx(3,i,:,:) = diff(dJ2_q3([q1,q2,q3,q4,q5,q6]),xmat(i));
        dJ_heel_dxx(4,i,:,:) = diff(dJ2_q4([q1,q2,q3,q4,q5,q6]),xmat(i));
        dJ_heel_dxx(5,i,:,:) = diff(dJ2_q5([q1,q2,q3,q4,q5,q6]),xmat(i));
        dJ_heel_dxx(6,i,:,:) = diff(dJ2_q6([q1,q2,q3,q4,q5,q6]),xmat(i));
    end
end


task{1,task_i} = @()matlabFunction(dV_dxx,'file','hessian/dV_dxx','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6]});task_i=task_i+1;
task{1,task_i} = @()matlabFunction(dG_dxx,'file','hessian/dG_dxx','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6]});task_i=task_i+1;
task{1,task_i} = @()matlabFunction(dM_dxx,'file','hessian/dM_dxx','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6]});task_i=task_i+1;
task{1,task_i} = @()matlabFunction(dJ_toe_dxx,'file','hessian/dJ_toe_dxx','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6]});task_i=task_i+1;
task{1,task_i} = @()matlabFunction(dJ_heel_dxx,'file','hessian/dJ_heel_dxx','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6]});task_i=task_i+1;  

syms q11 q21 q31 q41 q51 q61 qd11 qd21 qd31 qd41 qd51 qd61 u11 u21 u31 u41 u51 u61 fx_toe1 fy_toe1 fx_heel1 fy_heel1 real
syms q12 q22 q32 q42 q52 q62 qd12 qd22 qd32 qd42 qd52 qd62 u12 u22 u32 u42 u52 u62 fx_toe2 fy_toe2 fx_heel2 fy_heel2 real

% xmat = [q11, q21, q31, q41, q51, q61, qd11, qd21, qd31, qd41, qd51, qd61, u11, u21, u31, u41, u51, u61, fx_toe1, fy_toe1, fx_heel1, fy_heel1;
%         q12, q22, q32, q42, q52, q62, qd12, qd22, qd32, qd42, qd52, qd62, u12, u22, u32, u42, u52, u62, fx_toe2, fy_toe2, fx_heel2, fy_heel2].'; 
    
xmat2 = [q11, q21, q31, q41, q51, q61, qd11, qd21, qd31, qd41, qd51, qd61, u11, u21, u31, u41, u51, u61, fx_toe1, fy_toe1, fx_heel1, fy_heel1,... % this condition involves two time steps
        q12, q22, q32, q42, q52, q62, qd12, qd22, qd32, qd42, qd52, qd62, u12, u22, u32, u42, u52, u62, fx_toe2, fy_toe2, fx_heel2, fy_heel2];     
hess_hipLen = sym(zeros(3*param.numJ+4,3*param.numJ+4));
[~,gradceq_hipLen]=hipCon(xmat,param);

for i1=1:length(gradceq_hipLen)
    for i2=1:length(gradceq_hipLen) %this also only solve the upper triangle part
        hess_hipLen(i1,i2)=diff(gradceq_hipLen(i1,1),xmat2(i2));
    end
end
task{1,task_i} = @()matlabFunction(hess_hipLen,'file','hessian/h_hipLen','vars',{[q11,q21,q31,q41,q51,q61,qd11,qd21,qd31,qd41,qd51,qd61,u11,u21,u31,u41,u51,u61,fx_toe1,fy_toe1,fx_heel1,fy_heel1],...
                                                                                 [q12,q22,q32,q42,q52,q62,qd12,qd22,qd32,qd42,qd52,qd62,u12,u22,u32,u42,u52,u62,fx_toe2,fy_toe2,fx_heel2,fy_heel2],...
                                                                                  hipLen});task_i=task_i+1;
                                                                              
%% I believe it is too difficult to solve df_dxx, we should do something simpler
xmat = [q1, q2, q3, q4, q5, q6, qd1, qd2, qd3, qd4, qd5, qd6, u1, u2, u3, u4, u5, u6, fx_toe, fy_toe, fx_heel, fy_heel];
Mddq = getMddq(xmat.',param);
Mddq_dxx =sym(zeros(param.numJ*3+4,param.numJ*3+4,param.numJ));
for i1=1:param.numJ*3+4  
    for i2=1:param.numJ*3+4
        for i3=1:param.numJ
            dx = diff(Mddq(i3),xmat(i1));
            dxx = diff(dx,xmat(i2)); 
            Mddq_dxx(i1,i2,i3)=dxx;
        end
    end
end
task{1,task_i}=@()matlabFunction(Mddq_dxx,'file','hessian/Mddq_dxx','vars',{xmat,ank_stiff,knee_stiff}); task_i=task_i+1;



parfor i=1:length(task)
    task{1,i}();
end  


