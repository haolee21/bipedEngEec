function [f,dfx]=f_dyn(x,p,f_toe,f_heel,u)

q = x(1:p.numJ,1).';
dq = x(p.numJ+1:2*p.numJ,1).';

M = six_M(q(2),q(3),q(4),q(5),q(6));
G = six_G(q(1),q(2),q(3),q(4),q(5),q(6));
V = six_V(q(2),q(3),q(4),q(5),q(6),dq(1),dq(2),dq(3),dq(4),dq(5),dq(6));
J_toe = six_J(q(1),q(2),q(3),q(4),q(5),q(6));
J_heel = six_J2(q(1),q(2),q(3),q(4),q(5),q(6));


dtau_tend_1 = zeros(p.numJ*2,p.numJ);
dtau_tend_2 = zeros(p.numJ*2,p.numJ);
dtau_kne_tend_1 = zeros(p.numJ*2,p.numJ);
dtau_kne_tend_2 = zeros(p.numJ*2,p.numJ);
% ankle tendon
tau_tend_1 = p.ank_stiff*[pi/2-q(1),0,0,0,0,0];
dtau_tend_1(1,1) = -p.ank_stiff;
tau_tend_2 = p.ank_stiff*[0,0,0,0,0,-pi/2-q(6)];
dtau_tend_2(6,6)=-p.ank_stiff;
% knee tendon
tau_kne_tend_1 = [0,-p.knee_stiff*q(2),0,0,0,0];
tau_kne_tend_2 = [0,0,0,0,-p.knee_stiff*q(5),0];
dtau_kne_tend_1(2,2) = -p.knee_stiff;
dtau_kne_tend_2(5,5) = -p.knee_stiff;



ddq = (u.'-V-G+[f_toe.',0,0,0,0]*J_toe+[f_heel.',0,0,0,0]*J_heel-dq*eye(p.numJ)*p.joint_fri+tau_tend_1+tau_tend_2+tau_kne_tend_1+tau_kne_tend_2)/M;
d_joint_fri = [zeros(p.numJ);-eye(p.numJ)*p.joint_fri];
f = [dq,ddq].';
dVdx =  dV_dx(q(2),q(3),q(4),q(5),q(6),dq(1),dq(2),dq(3),dq(4),dq(5),dq(6));
dGdx =  [dG_dx(q(1),q(2),q(3),q(4),q(5),q(6));zeros(p.numJ)];

out_t_diag = repmat({ddq},1,p.numJ);
dMdx = blkdiag(out_t_diag{:})*[zeros(p.numJ);
                               dM_dx2(q(2),q(3),q(4),q(5),q(6));
                               dM_dx3(q(2),q(3),q(4),q(5),q(6));
                               dM_dx4(q(2),q(3),q(4),q(5),q(6));
                               dM_dx5(q(2),q(3),q(4),q(5),q(6));
                               dM_dx6(q(2),q(3),q(4),q(5),q(6))];
dMdx = [dMdx;zeros(p.numJ)];




df_ext_toe_diag = repmat({[f_toe.',0,0,0,0]},1,p.numJ);
df_ext_toe_dx = blkdiag(df_ext_toe_diag{:})*[dJ_q1(q);
                                           dJ_q2(q);
                                           dJ_q3(q);
                                           dJ_q4(q);
                                           dJ_q5(q);
                                           dJ_q6(q)]; %dJ^T/dx*Fext
df_ext_toe_dx=[df_ext_toe_dx;zeros(6)];
df_ext_heel_diag = repmat({[f_heel.',0,0,0,0]},1,p.numJ);
df_ext_heel_dx = blkdiag(df_ext_heel_diag{:})*[dJ2_q1(q);
                                            dJ2_q2(q);
                                            dJ2_q3(q);
                                            dJ2_q4(q);
                                            dJ2_q5(q);
                                            dJ2_q6(q)]; %dJ^T/dx*Fext

df_ext_heel_dx=[df_ext_heel_dx;zeros(6)];  


dfx = [[zeros(p.numJ);eye(p.numJ)],...
       (-dVdx-dGdx-dMdx+dtau_tend_1+dtau_tend_2+dtau_kne_tend_1+dtau_kne_tend_2+d_joint_fri+df_ext_toe_dx+df_ext_heel_dx)/M];

end