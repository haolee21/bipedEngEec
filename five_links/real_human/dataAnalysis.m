addpath result/
%% data analyzing 

% data1 does not have grf
% data2 has grf

data1 = load('04271301.mat').result;
data2 = load('04271302.mat').result;
time1 = 0:data1.param.sampT:data1.param.gaitT;
time2 = 0:data2.param.sampT:data2.param.gaitT;

figure(1);
subplot(2,2,1);
plot(time1,data1.x(15,:),time2,data2.x(15,:));
legend('with ankle-push','without');
title('rear knee');
ylabel('Torque (N.m)');
xlabel('Time (sec)');
subplot(2,2,3);
plot(time1,data1.x(14,:),time2,data2.x(14,:));
legend('with ankle-push','without');
title('rear hip');
ylabel('Torque (N.m)');
xlabel('Time (sec)');
subplot(2,2,2);
plot(time1,data1.x(12,:),time2,data2.x(12,:));
legend('with ankle-push','without');
title('front knee');
ylabel('Torque (N.m)');
xlabel('Time (sec)');
subplot(2,2,4);
plot(time1,data1.x(13,:),time2,data2.x(13,:));
legend('with ankle-push','without');
title('front hip');
ylabel('Torque (N.m)');
xlabel('Time (sec)');


%% test did we always touch the floor
figure(2);
sig = zeros(1,size(data1.x,2));
for i=1:size(data1.x,2)
    sig(i) = sigma_out(data1.x(1,i),data1.x(2,i),data1.x(3,i),data1.x(4,i),data1.x(5,i),data1.param.toe_th);
    
end
plot(time1,sig);
title('Ground Touching');
xlabel('Time (sec)');
ylabel('Ground Contact (%)');

