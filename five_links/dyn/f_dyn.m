function f = f_dyn(x)
%% calculate [dq;ddq]
% here z is 1 x n vector, represent the state at one specific of time

q = x(1:5);
dq = x(6:10);
u = x(11:15)

M = five_M(q(2),q(3),q(4),q(5));
V = five_V(q(2),q(3),q(4),q(5),dq(1),dq(2),dq(3),dq(4),dq(5));
G = five_G(q(1),q(2),q(3),q(4),q(5));

%% need to tranpose the result since we want a row vector
f = [dq,(M\(u.'-V*dq.'-G.')).'];


end