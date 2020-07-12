% check the gradient for discrete Lagrangian

clear;
clc;
addpath 'dyn'

x = load('x_contact_test.mat').x;
u = load('u_contact_test.mat').u;
p = load('param_contact_test.mat').param;
p.us = 0.8;
p.gndclear =-0.03;
f_toe = load('f_toe_test.mat').f_toe;
f_heel = load('f_heel_test.mat').f_heel;
x1 = [x(1:p.numJ,:);u;f_toe;f_heel;zeros(2,size(x,2))];
x2 = x1+rand(size(x1,1),size(x1,2))*0.0001;
% x2 = x1+1;
% x2(1:p.numJ,:)=x1(1:p.numJ,:);
% x2(p.numJ+1:end,:)=x1(p.numJ+1:end,:);
% x2(2*p.numJ+1:end,:) = x1(2*p.numJ+1:end,:);

[c1,gradc1] =dynConst_discrete(x1,p);
[c2,gradc2] = dynConst_discrete(x2,p);
dx = reshape(x2-x1,[size(x1,1)*size(x1,2),1]);
err = c2-c1-0.5*(gradc1+gradc2).'*dx;

disp(['dynConst err= ',num2str(norm(err)/norm(c2-c1))]); 


% check grf constraint
[c1,ceq1,gradc1,gradceq1] = grf_cons_d(x1,p);
[c2,ceq2,gradc2,gradceq2] = grf_cons_d(x2,p);

dx = reshape(x2-x1,[size(x1,1)*size(x1,2),1]);
err1 = c2-c1-0.5*(gradc1+gradc2).'*dx;
err2 = ceq2-ceq1-0.5*(gradceq1+gradceq2).'*dx;

disp(['grf c err=',num2str(norm(err1)/norm(c2-c1))]);
disp(['grf ceq err=',num2str(norm(err2)/norm(ceq2-ceq1))]);

% check start/end height constraint

[c1,gradc1]=initYPosCons(x1,p);
[c2,gradc2]=initYPosCons(x2,p);
dx = reshape(x2-x1,[size(x1,1)*size(x1,2),1]);
err = c2-c1-0.5*(gradc1+gradc2).'*dx;
disp(['initY ceq err=',num2str(norm(err)/norm(c2-c1))]);

% check yPosCon

[c1,gradc1]=yposCon(x1,p);
[c2,gradc2]=yposCon(x2,p);
dx = reshape(x2-x1,[size(x1,1)*size(x1,2),1]);
err = c2-c1-0.5*(gradc1+gradc2).'*dx;
disp(['yposCon err=',num2str(norm(err)/norm(c2-c1))]);

% check initYpos
[c1,gradc1]=initYPosCons(x1,p);
[c2,gradc2]=initYPosCons(x2,p);
dx = reshape(x2-x1,[size(x1,1)*size(x1,2),1]);
err = c2-c1-0.5*(gradc1+gradc2).'*dx;
disp(['initYpos err=',num2str(norm(err)/norm(c2-c1))]);

% check f_dyn
u =x1(p.numJ+1:2*p.numJ,:);
f_toe = x1(p.numJ*2+1:p.numJ*2+2,:);
f_heel= x1(p.numJ*2+3:p.numJ*2+4,:);



[f1,dfx1]=f_dyn(x1,p,f_toe,f_heel,u);
[f2,dfx2]=f_dyn(x2,p,f_toe,f_heel,u);
dx = reshape(x2(1:p.numJ,:)-x1(1:p.numJ,:),[p.numJ*size(x1,2),1]);
err = f2-f1-0.5*(dfx1+dfx2).'*dx;
disp(['f_dyn(dynLoss) err=',num2str(norm(err)/norm(f2-f1))]);


% test dynObj


[f1,dfx1]=dynObj(x1(1:p.numJ,:),p,u,f_toe,f_heel);
[f2,dfx2]=dynObj(x2(1:p.numJ,:),p,u,f_toe,f_heel);
dx = reshape(x2(1:p.numJ,:)-x1(1:p.numJ,:),[p.numJ*size(x1,2),1]);
err = f2-f1-0.5*(dfx1+dfx2).'*dx;
disp(['dynObj err=',num2str(norm(err)/norm(f2-f1))]);




