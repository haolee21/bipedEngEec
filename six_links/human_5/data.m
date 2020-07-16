clear;
close all;
% data analysis for this model

% baseline
base = load('07121516').result.x;
base_v = gradient(base);
% data1: knee efficient (5 times more)
data1 = load('07121703').result.x;
data1_v = gradient(data1);

% data2: knee efficient (10 times more)
data2 = load('07122303').result.x;
data2_v = gradient(data2);

% data3: efficient knees, back ankle (10 times)
data3 = load('07131038').result.x;
data3_v = gradient(data3);

time = linspace(0,100,size(base,2)*2);
% knee joint
figure(1);

subplot(3,1,1);
hold on;
plot(time,[base(2,:),-base(5,:)]*180/pi);
plot(time,[data1(2,:),-data1(5,:)]*180/pi);
plot(time,[data2(2,:),-data2(5,:)]*180/pi);
% plot(time,[data3(2,:),-data3(5,:)]*180/pi);
legend('base line','knee eff 5','knee eff 10');
title('knee angle');
hold off;
subplot(3,1,2);
hold on;
plot(time,[base(8,:),-base(11,:)]);
plot(time,[data1(8,:),-data1(11,:)]);
plot(time,[data2(8,:),-data2(11,:)]);
% plot(time,[data3(8,:),-data3(11,:)]);
% legend('base line','knee eff 5','knee eff 10');
title('knee torque');
hold off;
subplot(3,1,3);
work_base = [base(8,:),-base(11,:)].*[base_v(2,:),-base(5,:)];
work_data1 = [data1(8,:),-data1(11,:)].*[data1_v(2,:),-data1_v(5,:)];
work_data2 = [data2(8,:),-data2(11,:)].*[data2_v(2,:),-data2_v(5,:)];
% work_data3 = [data3(8,:),-data3(11,:)].*[data3_v(2,:),-data3_v(5,:)];
title('Knee Work');
hold on;
plot(time,work_base);
plot(time,work_data1);
plot(time,work_data2);
% plot(time,work_data3);
% legend('base line','knee eff 5','knee eff 10');

hold off;

% ankle
figure(2);
subplot(3,1,1);
hold on;
plot(time,[base(1,:),-base(6,:)]*180/pi-90);
plot(time,[data1(1,:),-data1(6,:)]*180/pi-90);
plot(time,[data2(1,:),-data2(6,:)]*180/pi-90);
% plot(time,[data3(1,:),-data3(6,:)]*180/pi);
legend('base line','knee eff 5','knee eff 10');
title('Ankle angle');
hold off;
subplot(3,1,2);
hold on;
plot(time,[base(7,:),-base(12,:)]);
plot(time,[data1(7,:),-data1(12,:)]);
plot(time,[data2(7,:),-data2(12,:)]);
% plot(time,[data3(7,:),-data3(12,:)]);
% legend('base line','knee eff 5','knee eff 10');
title('Ankle torque');
hold off;
subplot(3,1,3);
work_base = [base(7,:),-base(12,:)].*[base_v(2,:),-base(5,:)];
work_data1 = [data1(7,:),-data1(12,:)].*[data1_v(1,:),-data1_v(6,:)];
work_data2 = [data2(7,:),-data2(12,:)].*[data2_v(1,:),-data2_v(6,:)];
% work_data3 = [data3(7,:),-data3(12,:)].*[data3_v(1,:),-data3_v(6,:)];
title('Ankle Work');
hold on;
plot(time,work_base);
plot(time,work_data1);
plot(time,work_data2);
% plot(time,work_data3);
% legend('base line','knee eff 5','knee eff 10');

% hip
figure(3);
subplot(3,1,1);
hold on;
plot(time,[base(3,:),-pi-base(4,:)]*180/pi*-1);
plot(time,[data1(3,:),-pi-data1(4,:)]*180/pi*-1);
plot(time,[data2(3,:),-pi-data2(4,:)]*180/pi*-1);
% plot(time,[data3(3,:),-pi-data3(4,:)]*180/pi*-1);
legend('base line','knee eff 5','knee eff 10');
title('Hip angle');
hold off;
subplot(3,1,2);
hold on;
plot(time,[base(9,:),-base(10,:)]*-1);
plot(time,[data1(9,:),-data1(10,:)]*-1);
plot(time,[data2(9,:),-data2(10,:)]*-1);
% plot(time,[data3(9,:),-data3(10,:)]*-1);
% legend('base line','knee eff 5','knee eff 10');
title('Hip torque');
hold off;
subplot(3,1,3);
work_base = [base(9,:),-base(10,:)].*[base_v(2,:),-base(5,:)];
work_data1 = [data1(9,:),-data1(10,:)].*[data1_v(3,:),-data1_v(4,:)];
work_data2 = [data2(9,:),-data2(10,:)].*[data2_v(3,:),-data2_v(4,:)];
% work_data3 = [data3(9,:),-data3(10,:)].*[data3_v(3,:),-data3_v(4,:)];
title('Hip Work');
hold on;
plot(time,work_base);
plot(time,work_data1);
plot(time,work_data2);
% plot(time,work_data3);
% legend('base line','knee eff 5','knee eff 10');

% GRF
time_grf = linspace(50,100,size(base,2));
figure(4);
subplot(2,1,1);
hold on;
plot(time_grf,base(14,:)/9.8/75);
plot(time_grf,data1(14,:)/9.8/75);
plot(time_grf,data2(14,:)/9.8/75);
% plot(time_grf,data3(14,:)/9.8/75);
xlim([0,100]);
title('Toe GRF: Fy');
legend('base line','knee eff 5','knee eff 10','knee, ankle eff 10');
hold off;
subplot(2,1,2);
hold on;
plot(time_grf,base(16,:)/9.8/75);
plot(time_grf,data1(16,:)/9.8/75);
plot(time_grf,data2(16,:)/9.8/75);
% plot(time_grf,data3(16,:)/9.8/75);
xlim([0,100]);
title('Heel GRF: Fy');
legend('base line','knee eff 5','knee eff 10','knee, ankle eff 10');
hold off;
