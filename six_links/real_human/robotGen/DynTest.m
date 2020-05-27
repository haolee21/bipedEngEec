%%
%  Test the accuracy of dynamic/grf equations
%
%   
%
%

addpath 'dyn/'
addpath 'grf/'

q1 = 70;
q2 = -5;
q3 =45;
q4 = -180;
q5 =2;
q6 = -180-q1-q2-q3-q4-q5;%+atan2d(param.heel_h,0.26);%0.26 is feet length

p.numJ=6;
p.toe_th =-6e-2;

tspan = [0 2];
ic=[q1/180*pi,q2/180*pi,q3/180*pi,q4/180*pi,q5/180*pi,q6/180*pi,zeros(1,p.numJ)];
opts = odeset('RelTol',1e-2,'AbsTol',1e-4);

ut = linspace(0,5,200);
u = zeros(p.numJ,length(ut));


fdyn =@(t,x)robotTestDyn(t,x,ut,u,p);

[t,y]=ode45(fdyn,tspan,ic,opts);











function dydt = robotTestDyn(t,x,ut,u,p)
numJ = p.numJ;
dydt = zeros(2*numJ,1);
u = interp1(ut,u.',t);
q = x(1:numJ).';
dq = x(numJ+1:2*numJ).';
M_ext = Mext_toe(q);
beta_grf = beta_grf_toe(q,p.toe_th);

M = six_M(q(2),q(3),q(4),q(5),q(6));
G = six_G(q(1),q(2),q(3),q(4),q(5),q(6));
V = six_V(q(2),q(3),q(4),q(5),q(6),dq(1),dq(2),dq(3),dq(4),dq(5),dq(6));

extF = M_ext*(u-G).';
if(extF(2)>0)
    extF=zeros(size(extF,1),size(extF,2));
    %M_ext = zeros(size(M_ext,1),size(M_ext,2));
end
dydt(1:numJ,1) = dq.';
dydt(numJ+1:2*numJ,1)= ((u-V-G-extF.'*beta_grf.')/M).';

end