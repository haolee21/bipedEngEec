function [f,dfu] = f_forward(u,q,dq,p)



M = six_M(q(2),q(3),q(4),q(5),q(6));
G = six_G(q(1),q(2),q(3),q(4),q(5),q(6));
V = six_V(q(2),q(3),q(4),q(5),q(6),dq(1),dq(2),dq(3),dq(4),dq(5),dq(6));

% [tau_ext_toe,dTau_ext_toe,~,~]=toe_grf(x.',p);
% [tau_ext_heel,dTau_ext_heel,~,~]=heel_grf(x.',p);


% ankle tendon
% 


tau_tend_1 = p.ank_stiff*[pi/2-q(1),0,0,0,0,0];



tau_tend_2 = p.ank_stiff*[0,0,0,0,0,-3*pi/2-q(6)];

% knee tendon


tau_kne_tend_1 = [0,-p.knee_stiff*q(2),0,0,0,0];
tau_kne_tend_2 = [0,0,0,0,-p.knee_stiff*q(5),0];




ddq = (u.'-V-G-dq.'*eye(p.numJ)*p.joint_fri+tau_tend_1+tau_tend_2+tau_kne_tend_1+tau_kne_tend_2)/M;

f = [dq.',ddq];









% grad =(dTaudx-dVdx-dGdx-dMdx+df_ext_toe1+df_ext_toe2+df_ext_heel1+df_ext_heel2-[zeros(p.numJ);p.joint_fri*eye(p.numJ);zeros(p.numJ);zeros(4,p.numJ)]+dtau_tend_1+dtau_tend_2+dtau_kne_tend_1+dtau_kne_tend_2)/M;

f = f.'; % we later wanna us it as column vector

% dfs = [[zeros(p.numJ);eye(p.numJ)],grad(1:p.numJ*2,:)]; %zeros(p.numJ) is because the output of f is [dq;ddq], 
dfu = [zeros(p.numJ),eye(p.numJ)/M];

 

end