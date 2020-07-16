function dydt = sim_dyn(t,x,ut,u,f,p)
u = interp1(ut,u.',t).';
f = interp1(ut,f.',t).';
M = six_M(x(2,1),x(3,1),x(4,1),x(5,1),x(6,1));
V = six_V(x(2,1),x(3,1),x(4,1),x(5,1),x(6,1),x(7,1),x(8,1),x(9,1),x(10,1),x(11,1),x(12,1));
G = six_G(x(1,1),x(2,1),x(3,1),x(4,1),x(5,1),x(6,1));
J = six_J(x(1,1),x(2,1),x(3,1),x(4,1),x(5,1),x(6,1));
J = J(1:2,:);

tend1 = [(pi/2-x(1,1))*p.ank_stiff,0,0,0,0,0];
tend2 = [0,0,0,0,0,(-pi/2-x(6,1))*p.ank_stiff];

ddq = (u.'-G-V+f.'*J+tend1+tend2)/M;
dydt = [x(7:12,1);ddq.'];

end