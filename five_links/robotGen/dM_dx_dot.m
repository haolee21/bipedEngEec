function out = dM_dx_dot(q,beta)


syms q2 q3 q4 q5

out(1,:) = zeros(1,5);
out(2,:) = beta*double(subs(diff(eye(5)\five_M(q2,q(3),q(4),q(5)),q2),q2,q(2)));
out(3,:) = beta*double(subs(diff(eye(5)\five_M(q(2),q3,q(4),q(5)),q3),q3,q(3)));
out(4,:) = beta*double(subs(diff(eye(5)\five_M(q(2),q(3),q4,q(5)),q4),q4,q(4)));
out(5,:) = beta*double(subs(diff(eye(5)\five_M(q(2),q(3),q(4),q5),q5),q5,q(5)));

out(6:15,:) = zeros(10,5);


end