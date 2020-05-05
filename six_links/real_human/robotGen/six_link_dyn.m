function dydt = six_link_dyn(t,x,ut,u)
% x is the state vector, [q1,q2,q3,q4,q5,dq1,dq2,dq3,dq4,dq5]
u = interp1(ut,u,t);
q = x(1:6);
dq=x(7:12);

M = six_M(q(2),q(3),q(4),q(5),q(6));
V = six_V(q(2),q(3),q(4),q(5),q(6),dq(1),dq(2),dq(3),dq(4),dq(5),dq(6)).';
G = six_G(q(1),q(2),q(3),q(4),q(5),q(6));
grf = beta_grf(q(1),q(2),q(3),q(4),q(5),q(6),0.01);

ddq = M\(u.'-V-G.');

dydt(1,1)=x(7);
dydt(2,1)=x(8);
dydt(3,1)=x(9);
dydt(4,1)=x(10);
dydt(5,1)=x(11);
dydt(6,1)=x(12);

dydt(7,1)=ddq(1);
dydt(8,1)=ddq(2);
dydt(9,1)=ddq(3);
dydt(10,1)=ddq(4);
dydt(11,1)=ddq(5);
dydt(12,1)=ddq(6);



end