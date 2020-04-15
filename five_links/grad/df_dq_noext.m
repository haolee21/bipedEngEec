function grad = df_dq_noext(x)
% u is the joint torque with dim=5
% external force is excluded since we will calculate it later
q = x(1:5);
dq = x(6:10);
u = x(11:15);


syms q1 q2 q3 q4 q5 real

% calculating equations with 5 symbols is too difficult, we will take
% derivative of each variable individually

% calculate df/dq1

M = five_M(q(2),q(3),q(4),q(5));
G = five_G(q1,q(2),q(3),q(4),q(5));
V = five_V(q(2),q(3),q(4),q(5),dq(1),dq(2),dq(3),dq(4),dq(5));

ddq = M\(u.'-V*dq.'-G.'); %no external force since it is generated with torque in previous time, is not function of joint angles now
grad_q1 = double(subs(diff(ddq,q1),q1,q(1))).';

% df/dq2
M = five_M(q2,q(3),q(4),q(5));
G = five_G(q(1),q2,q(3),q(4),q(5));
V = five_V(q2,q(3),q(4),q(5),dq(1),dq(2),dq(3),dq(4),dq(5));

ddq = M\(u.'-V*dq.'-G.'); 
grad_q2 = double(subs(diff(ddq,q2),q2,q(2))).';

% df/dq3
M = five_M(q(2),q3,q(4),q(5));
G = five_G(q(1),q(2),q3,q(4),q(5));
V = five_V(q(2),q3,q(4),q(5),dq(1),dq(2),dq(3),dq(4),dq(5));

ddq = M\(u.'-V*dq.'-G.'); 
grad_q3 = double(subs(diff(ddq,q3),q3,q(3))).';

% df/dq4
M = five_M(q(2),q(3),q4,q(5));
G = five_G(q(1),q(2),q(3),q4,q(5));
V = five_V(q(2),q(3),q4,q(5),dq(1),dq(2),dq(3),dq(4),dq(5));

ddq = M\(u.'-V*dq.'-G.'); 
grad_q4 = double(subs(diff(ddq,q4),q4,q(4))).';


% df/dq5
M = five_M(q(2),q(3),q(4),q5);
G = five_G(q(1),q(2),q(3),q(4),q5);
V = five_V(q(2),q(3),q(4),q5,dq(1),dq(2),dq(3),dq(4),dq(5));

ddq = M\(u.'-V*dq.'-G.'); 
grad_q5 = double(subs(diff(ddq,q5),q5,q(5))).';

% df/ddq

syms dq1 dq2 dq3 dq4 dq5 real
M = five_M(q(2),q(3),q(4),q(5));
G = five_G(q(1),q(2),q(3),q(4),q(5));
V = five_V(q(2),q(3),q(4),q(5),dq1,dq2,dq3,dq4,dq5);

ddq = M\(u.'-V*dq.'-G.'); 
grad_dq1 = double(subs(diff(ddq,dq1),dq1,dq(1))).';
grad_dq2 = double(subs(diff(ddq,dq2),dq2,dq(2))).';
grad_dq3 = double(subs(diff(ddq,dq3),dq3,dq(3))).';
grad_dq4 = double(subs(diff(ddq,dq4),dq4,dq(4))).';
grad_dq5 = double(subs(diff(ddq,dq5),dq5,dq(5))).';


% combine and add [0 1] row in the matrix since the states are [q;dq]
grad = [grad_q1;grad_q2;grad_q3;grad_q4;grad_q5;grad_dq1;grad_dq2;grad_dq3;grad_dq4;grad_dq5;M\eye(5)];
end