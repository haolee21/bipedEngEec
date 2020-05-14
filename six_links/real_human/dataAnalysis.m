%% data analyzing 
close all;
data1 = load('05122331.mat').result;

% data3 = load('05041233.mat').result;
% data4 = load('05041325.mat').result;
% drawRobot_video(data1,'05041150');
% drawRobot_video(data2,'05041200');
% drawRobot_video(data2,'05041233');
% drawRobot_video(data4,'05041325');


numJ = data1.param.numJ;
u1 = data1.x(numJ*2+1:numJ*3,:);
eng1 = sum(u1.^2,'all');
u1_no_ank = u1(1:end-1,:);
eng1_no_ank = sum(u1_no_ank.^2,'all');

time = 0:data1.param.sampT:data1.param.gaitT;

figure(2);
plot(time,u1);
legend('front ank','fronk knee','front hip','back hip','back knee','back ank');
ylim([-20,20]);


% 
% data2 = load('05111318.mat').result;
% numJ = data2.param.numJ;
% u2 = data2.x(numJ*2+1:numJ*3,:);
% eng2 = sum(u2.^2,'all');
% u2_no_ank = u2(1:end-1,:);
% eng2_no_ank = sum(u2_no_ank.^2,'all');
% time2 = 0:data2.param.sampT:data2.param.gaitT;
% figure(3);
% plot(time,u2);
% legend('front ank','fronk knee','front hip','back hip','back knee','back ank');
% ylim([-20,20]);



figure(4);
time1_all = 0:data1.param.sampT:data1.param.gaitT*2;
hold on;

plot(time,data1.x(2,:)/pi*180);
plot(time,-data1.x(5,:)/pi*180);
legend('p1','p2');
hold off;
figure(5);

comb = [-data1.x(5,:),data1.x(2,1:end-1)]*180/pi;

% comb = [comb(159:end),comb(1:158)];
time2 = 1:1:size(comb,2);
plot(time2,comb);

% 
% figure(2);
% 
% subplot(2,1,1);
% 
% plot(time,u1(numJ,:));
% title('With Ankle Push');
% ylim([-2,1]);
% ylabel('Torque (N.m)');
% subplot(2,1,2);
% plot(time,u2(numJ,:));
% title('Without Ankle Push');
% ylim([-2,1]);
% ylabel('Torque (N.m)');
% xlabel('Time (sec)');
% 
% 
% 
% u3 = data3.x(numJ*2+1:numJ*3,:);
% figure(3);
% 
% subplot(2,1,1);
% 
% plot(time,u3(numJ,:));
% title('With Ankle Push');
% ylim([-2,1]);
% ylabel('Torque (N.m)');
% subplot(2,1,2);
% plot(time,u2(numJ,:));
% title('Without Ankle Push');
% ylim([-2,1]);
% ylabel('Torque (N.m)');
% xlabel('Time (sec)');