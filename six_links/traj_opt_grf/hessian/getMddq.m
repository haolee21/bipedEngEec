function Mddq_out = getMddq(x,p)


q = x(1:p.numJ,1).';
dq = x(p.numJ+1:2*p.numJ,1).';
u = x(p.numJ*2+1:3*p.numJ,1).';
fext_toe = [x(p.numJ*3+1:p.numJ*3+2,1).',0,0,0,0];
fext_heel= [x(p.numJ*3+3:p.numJ*3+4,1).',0,0,0,0];


G = six_G(q(1),q(2),q(3),q(4),q(5),q(6));
V = six_V(q(2),q(3),q(4),q(5),q(6),dq(1),dq(2),dq(3),dq(4),dq(5),dq(6));
J_toe = six_J(q(1),q(2),q(3),q(4),q(5),q(6));
J_heel = six_J2(q(1),q(2),q(3),q(4),q(5),q(6));


% 
tau_tend_1 = (0.5*tanh(200*(q(1)-pi/2))+0.5)*p.ank_stiff*[pi/2-q(1),0,0,0,0,0];

tau_tend_2 = (0.5*tanh(200*(-3*pi/2-q(6)))+0.5)*p.ank_stiff*[0,0,0,0,0,-3*pi/2-q(6)];
% knee tendon

tau_kne_tend_1 = [0,-p.knee_stiff*q(2),0,0,0,0];
tau_kne_tend_2 = [0,0,0,0,-p.knee_stiff*q(5),0];
Mddq_out = (u-V-G+fext_toe*J_toe+fext_heel*J_heel-dq*eye(p.numJ)*p.joint_fri+tau_tend_1+tau_tend_2+tau_kne_tend_1+tau_kne_tend_2);


end