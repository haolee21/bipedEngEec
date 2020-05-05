%% data analyzing 

data1 = load('05041150.mat').result;
data2 = load('05041200.mat').result;
data3 = load('05041233.mat').result;
data4 = load('05041325.mat').result;
% drawRobot_video(data1,'05041150');
% drawRobot_video(data2,'05041200');
% drawRobot_video(data2,'05041233');
% drawRobot_video(data4,'05041325');


numJ = data1.param.numJ;
u1 = data1.x(numJ*2+1:numJ*3,:);
eng1 = sum(u1.^2,'all');
time = 0:data1.param.sampT:data1.param.gaitT;



u2 = data2.x(numJ*2+1:numJ*3,:);
eng2 = sum(u2.^2,'all');

figure(2);

subplot(2,1,1);

plot(time,u1(numJ,:));
title('With Ankle Push');
ylim([-2,1]);
ylabel('Torque (N.m)');
subplot(2,1,2);
plot(time,u2(numJ,:));
title('Without Ankle Push');
ylim([-2,1]);
ylabel('Torque (N.m)');
xlabel('Time (sec)');



u3 = data3.x(numJ*2+1:numJ*3,:);
figure(3);

subplot(2,1,1);

plot(time,u3(numJ,:));
title('With Ankle Push');
ylim([-2,1]);
ylabel('Torque (N.m)');
subplot(2,1,2);
plot(time,u2(numJ,:));
title('Without Ankle Push');
ylim([-2,1]);
ylabel('Torque (N.m)');
xlabel('Time (sec)');