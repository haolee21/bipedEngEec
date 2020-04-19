%% Test the robot dynamic accuracy
addpath ../plotRobot/
ut=0:0.1:10;

u1=0.01*sin(ut+10/180*pi);
u2=0.01*sin(ut+20/180*pi);
u3=0.01*sin(ut+30/180*pi);
u4=0.01*sin(ut+40/180*pi);
u5=0.00001*sin(ut+50/180*pi);

% u1=zeros(size(u1));
% u2=zeros(size(u1));
% u3=zeros(size(u1));
% u4=zeros(size(u1));
% u5=zeros(size(u1));

u = [u1.',u2.',u3.',u4.',u5.'];

x0 =[70/180*pi,-10/180*pi,45/180*pi,-pi,20/180*pi,0,0,0,0,0];

opts = odeset('RelTol',1e-2,'AbsTol',1e-4);
tspan=[0,10];
sol = ode45(@(t,x)five_link_dyn(t,x,ut,u),tspan,x0,opts);

fig=figure(1);
drawRobot_self(sol.y,fig);