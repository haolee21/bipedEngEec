function dydt = five_link_dyn(t,x,ut,u)
% x is the state vector, [q1,q2,q3,q4,q5,dq1,dq2,dq3,dq4,dq5]
u = interp1(ut,u,t);
q = x(1:5);
dq=x(6:10);

M = five_M_simp(q(2),q(3),q(4),q(5));
V = (five_V(q(2),q(3),q(4),q(5),dq(1),dq(2),dq(3),dq(4),dq(5)))*dq;
G = five_G(q(1),q(2),q(3),q(4),q(5));


ddq = M\(u.'-V-G.');

dydt(1,1)=x(6);
dydt(2,1)=x(7);
dydt(3,1)=x(8);
dydt(4,1)=x(9);
dydt(5,1)=x(10);
dydt(6,1)=ddq(1);
dydt(7,1)=ddq(2);
dydt(8,1)=ddq(3);
dydt(9,1)=ddq(4);
dydt(10,1)=ddq(5);



end