%% Test the robot dynamic accuracy
addpath ../plotRobot/
addpath dyn/
ut=0:0.1:10;

u1=0*sin(ut+10/180*pi);
u2=0*sin(ut+20/180*pi);
u3=0*sin(ut+30/180*pi);
u4=0*sin(ut+40/180*pi);
u5=0*sin(ut+50/180*pi);
u6=0*sin(ut+50/180*pi);
% u1=zeros(size(u1));
% u2=zeros(size(u1));
% u3=zeros(size(u1));
% u4=zeros(size(u1));
% u5=zeros(size(u1));

u = [u1.',u2.',u3.',u4.',u5.',u6.'];

% x0 =[70/180*pi,-10/180*pi,45/180*pi,-pi,20/180*pi,-125/180*pi,...
%     0,0,0,0,0,0];
x0 =[10/180*pi,0,45/180*pi,-90/180*pi,0,1/180*pi,...
    0,0,0,0,0,0];
fig=figure(1);
drawRobot_self(x0.',fig)

opts = odeset('RelTol',1e-2,'AbsTol',1e-4);
tspan=[0,10];
sol = ode45(@(t,x)six_link_dyn(t,x,ut,u),tspan,x0,opts);
% 
% fig=figure(1);
% drawRobot_self(sol.y,fig);