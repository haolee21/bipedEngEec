%% This will calculate the hessian of the object function and nonlinear constraint

syms numJ toe_th head_h fri_coeff gaitT sampT init_y heel_h
syms foot_l dmax cmax k us ud

param.numJ=6;
param.toe_th =toe_th;
param.head_h = head_h ; %the head should be at least 1.6m
param.fri_coeff=fri_coeff;
param.gaitT = gaitT;
param.sampT = sampT;
param.init_y = init_y; %initial feet height
param.heel_h = heel_h; %this is fix in the model parameter
param.foot_l = foot_l;
param.dmax =dmax;
param.cmax=cmax;
param.k=k;
param.us=us;
param.ud=ud;

%% hipVelCon

syms q1 q2 q3 q4 q5 q6 qd1 qd2 qd3 qd4 qd5 qd6 u1 u2 u3 u4 u5 u6
xmat = [q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6,u1,u2,u3,u4,u5,u6];
cHipVelCon = hipVelCon(xmat.',param);
hessian_hipVelCon = sym(zeros(18,18));
for i1=1:length(xmat)
    for i2=1:length(xmat)
        hessian_hipVelCon(i1,i2) = diff(diff(cHipVelCon,xmat(i1)),xmat(i2));
        
        
    end
end
matlabFunction(hessian_hipVelCon,'file','h_hipVelCon','vars',{xmat});


%% dynConst

cDynConst = dynConst(xmat.',param);








