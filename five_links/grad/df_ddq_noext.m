function grad = df_ddq_noext(x)
q = x(1:5);
dq = x(6:10);
u = x(11:15);
syms dq1 dq2 dq3 dq4 dq5


M = five_M(q(2),q(3),q(4),q(5));
G = five_G(q(1),q(2),q(3),q(4),q(5));
V = five_V(q(2),q(3),q(4),q(5),dq1,dq2,dq3,dq4,dq5);

ddq = M\(u.'-V*dq.'-G.'); 
grad_dq1 = double(subs(diff(ddq,dq1),dq1,dq(1)));
grad_dq1 = double(subs(diff(ddq,dq1),dq1,dq(1)));
grad_dq1 = double(subs(diff(ddq,dq1),dq1,dq(1)));
grad_dq1 = double(subs(diff(ddq,dq1),dq1,dq(1)));
grad_dq1 = double(subs(diff(ddq,dq1),dq1,dq(1)));
grad = grad_dq1;
end