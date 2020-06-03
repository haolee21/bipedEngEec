%% data analyzing 

data1 = load('04201837.mat').result;
data2 = load('04201122.mat').result;

fig = figure(1);
% drawRobot_self(data1.x,fig);

% drawRobot_video(data1,'04201837');

figure(2);
subplot(2,2,1);
plot(data1.x(15,:));
title('with ankle-push, rear knee');
subplot(2,2,3);
plot(data1.x(14,:));
title('with ankle-push, rear hip');



% drawRobot_self(data2.x,fig);
% drawRobot_video(data2,'04201122');
subplot(2,2,2);
plot(data2.x(15,:));
title('without ankle-push, rear knee');
subplot(2,2,4);
plot(data2.x(14,:));
title('without ankle-push, rear hip');


figure(3);
subplot(2,2,1);
plot(data1.x(12,:));
title('with ankle-push, front knee');
subplot(2,2,3);
plot(data1.x(13,:));
title('with ankle-push, front hip');



% drawRobot_self(data2.x,fig);
% drawRobot_video(data2,'04201122');
subplot(2,2,2);
plot(data2.x(12,:));
title('without ankle-push, front knee');
subplot(2,2,4);
plot(data2.x(13,:));
title('without ankle-push, front hip');


figure(4);
subplot(2,1,1);
plot(data1.x(1,:));
subplot(2,1,2);
plot(data2.x(1,:));




