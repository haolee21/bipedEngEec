%% analysis the data of this model
clear;
addpath robotGen/
close all;
%%

data1 = load('06090244.mat').result;
numJ = data1.param.numJ;
u1 = data1.x(numJ*2+1:numJ*3,:);
eng1 = sum(u1.^2,'all');


time = 0:data1.param.sampT:data1.param.gaitT;

% find heel strike point
strike_time =[0,0];
toe_strike = zeros(1,length(time));
heel_strike = zeros(1,length(time));
for i=1:floor(size(data1.x,2)/2)
    curX = data1.x(:,i);
    grf_toe =sigma_toe(curX.',data1.param.toe_th);
    if(grf_toe>0.001)
        toe_strike(i)=1;
    end
end
for i=floor(size(data1.x,2)/2):size(data1.x,2)
    curX = data1.x(:,i);
    grf_heel = sigma_heel(curX.',data1.param.toe_th);
    if(grf_heel>0.001)
        heel_strike(i)=1;
    end
end

heel_s_idx = find(heel_strike);
qKne = [-1*data1.x(5,heel_s_idx(1):end),data1.x(2,:),-1*data1.x(5,1:heel_s_idx(1)-1)];
qAnk = [-1*data1.x(6,heel_s_idx(1):end),data1.x(1,:),-1*data1.x(6,1:heel_s_idx(1)-1)];
qHip = [data1.x(4,heel_s_idx(1):end)+pi,-data1.x(3,:),data1.x(4,1:heel_s_idx(1)-1)+pi];

uKne = [-1*data1.x(5+2*numJ,heel_s_idx(1):end),data1.x(2+2*numJ,:),-1*data1.x(5+2*numJ,1:heel_s_idx(1)-1)];
uAnk = [-1*data1.x(6+2*numJ,heel_s_idx(1):end),data1.x(1+2*numJ,:),-1*data1.x(6+2*numJ,1:heel_s_idx(1)-1)];
uHip = [data1.x(4+2*numJ,heel_s_idx(1):end),-data1.x(3+2*numJ,:),data1.x(4+2*numJ,1:heel_s_idx(1)-1)];


dq = [-1*data1.x(4+numJ:6+numJ,heel_s_idx(end):end),data1.x(1+numJ:3+numJ,:),data1.x(4+numJ:6+numJ,1:heel_s_idx(end-1))];
gait_duty = linspace(0,100,size(data1.x,2)*2);

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

