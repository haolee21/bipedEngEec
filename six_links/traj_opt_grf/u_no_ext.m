function u = u_no_ext(x,ddq,p)

q = x(1:p.numJ,1).';
dq = x(p.numJ+1:2*p.numJ,1).';

M = six_M(q(2),q(3),q(4),q(5),q(6));
G = six_G(q(1),q(2),q(3),q(4),q(5),q(6));
V = six_V(q(2),q(3),q(4),q(5),q(6),dq(1),dq(2),dq(3),dq(4),dq(5),dq(6));

tau_tend_1 = p.ank_stiff*[pi/2-q(1),0,0,0,0,0];
tau_tend_2 = p.ank_stiff*[0,0,0,0,0,-pi/2-q(6)];

tau_kne_tend_1 = [0,-p.knee_stiff*q(2),0,0,0,0];
tau_kne_tend_2 = [0,0,0,0,-p.knee_stiff*q(5),0];


u = (ddq.'*M+V+G+dq*eye(p.numJ)*p.joint_fri-tau_tend_1-tau_tend_2-tau_kne_tend_1-tau_kne_tend_2).';

end