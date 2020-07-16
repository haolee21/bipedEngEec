%% analysis the data of this model
clear;
addpath robotGen/
close all;
%%

data1 = load('06171151.mat').result;
numJ = data1.param.numJ;
u1 = data1.x(numJ*2+1:numJ*3,:);
eng1 = sum(u1.^2,'all');


time = 0:data1.param.sampT:data1.param.gaitT;


qKne = [-1*data1.x(5,:),data1.x(2,2:end)];
qAnk = [-1*data1.x(6,:),data1.x(1,2:end)];
qHip = [data1.x(4,:)+pi,-data1.x(3,2:end)];

uKne = [-1*data1.x(5+2*numJ,:),data1.x(2+2*numJ,2:end)];
uAnk = [-1*data1.x(6+2*numJ,:),data1.x(1+2*numJ,2:end)];
uHip = [data1.x(4+2*numJ,:),-data1.x(3+2*numJ,2:end)];



gait_duty = linspace(0,100,size(data1.x,2)*2-1);

figure(1);
subplot(2,1,1);
plot(gait_duty,qKne*180/pi);
ylabel('Angle (deg)');
title('Knee');
subplot(2,1,2);
plot(gait_duty,uKne/75);
ylabel('Torque/kg (Nm/kg)');

figure(2);
subplot(2,1,1);
plot(gait_duty,qAnk*180/pi-105);
title('Ankle');
ylabel('Angle (deg)');
subplot(2,1,2);
plot(gait_duty,uAnk/75);
ylabel('Torque/kg (Nm/kg)');

figure(3);
subplot(2,1,1);

plot(gait_duty,qHip*180/pi);
title('Hip');
ylabel('Angle (deg)');
subplot(2,1,2);
plot(gait_duty,uHip/75);
ylabel('Torque/kg (Nm/kg)');

