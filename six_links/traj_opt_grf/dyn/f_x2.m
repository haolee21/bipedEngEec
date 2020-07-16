function [f,dfs,dfu,df_fext] = f_x2(x,p)


q = x(1:p.numJ,1).';
dq = x(p.numJ+1:2*p.numJ,1).';
u = x(p.numJ*2+1:3*p.numJ,1).';
fext_toe = [x(p.numJ*3+1:p.numJ*3+2,1).',0,0,0,0];
fext_heel= [x(p.numJ*3+3:p.numJ*3+4,1).',0,0,0,0];

M = six_M(q(2),q(3),q(4),q(5),q(6));
G = six_G(q(1),q(2),q(3),q(4),q(5),q(6));
V = six_V(q(2),q(3),q(4),q(5),q(6),dq(1),dq(2),dq(3),dq(4),dq(5),dq(6));
J_toe = six_J(q(1),q(2),q(3),q(4),q(5),q(6));
J_heel = six_J2(q(1),q(2),q(3),q(4),q(5),q(6));
% [tau_ext_toe,dTau_ext_toe,~,~]=toe_grf(x.',p);
% [tau_ext_heel,dTau_ext_heel,~,~]=heel_grf(x.',p);


% ankle tendon
% tau_tend_1 = [0,0,0,0,0,0];
% tau_tend_2 = [0,0,0,0,0,0];
if(~isnumeric(x))
    dtau_tend_1 = sym(zeros(p.numJ*2,p.numJ));
    dtau_tend_2 = sym(zeros(p.numJ*2,p.numJ));
    dtau_kne_tend_1 = sym(zeros(p.numJ*2,p.numJ));
    dtau_kne_tend_2 = sym(zeros(p.numJ*2,p.numJ));
else
    dtau_tend_1 = zeros(p.numJ*2,p.numJ);
    dtau_tend_2 = zeros(p.numJ*2,p.numJ);
    dtau_kne_tend_1 = zeros(p.numJ*2,p.numJ);
    dtau_kne_tend_2 = zeros(p.numJ*2,p.numJ);
end
% 
% tau_tend_1 = (0.5*tanh(200*(q(1)-pi/2))+0.5)*p.ank_stiff*[pi/2-q(1),0,0,0,0,0];
% dtau_tend_1(1,1) = p.ank_stiff*(100*tanh(200*q(1) - 100*pi)^2 - 100)*(q(1) - pi/2) - p.ank_stiff*(tanh(200*q(1) - 100*pi)/2 + 1/2);

tau_tend_1 = p.ank_stiff*[pi/2-q(1),0,0,0,0,0];
dtau_tend_1(1,1) = -p.ank_stiff;

% tau_tend_2 = (0.5*tanh(200*(-3*pi/2-q(6)))+0.5)*p.ank_stiff*[0,0,0,0,0,-3*pi/2-q(6)];
% dtau_tend_2(6,6)=p.ank_stiff*(tanh(200*q(6) + 300*pi)/2 - 1/2) - p.ank_stiff*(100*tanh(200*q(6) + 300*pi)^2 - 100)*(q(6) + (3*pi)/2);
tau_tend_2 = p.ank_stiff*[0,0,0,0,0,-3*pi/2-q(6)];
dtau_tend_2(6,6)=-p.ank_stiff;
% knee tendon


tau_kne_tend_1 = [0,-p.knee_stiff*q(2),0,0,0,0];
tau_kne_tend_2 = [0,0,0,0,-p.knee_stiff*q(5),0];
dtau_kne_tend_1(2,2) = -p.knee_stiff;
dtau_kne_tend_2(5,5) = -p.knee_stiff;



ddq = (u-V-G+fext_toe*J_toe+fext_heel*J_heel-dq*eye(p.numJ)*p.joint_fri+tau_tend_1+tau_tend_2+tau_kne_tend_1+tau_kne_tend_2)/M;
d_joint_fri =[zeros(p.numJ);-eye(p.numJ)*p.joint_fri];

f = [dq,ddq];


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
df_ext_toe_diag = repmat({fext_toe},1,p.numJ);
df_ext_toe_dx = blkdiag(df_ext_toe_diag{:})*[dJ_q1(q);
                                           dJ_q2(q);
                                           dJ_q3(q);
                                           dJ_q4(q);
                                           dJ_q5(q);
                                           dJ_q6(q)]; %dJ^T/dx*Fext
df_ext_toe_dx=[df_ext_toe_dx;zeros(6)];
df_ext_heel_diag = repmat({fext_heel},1,p.numJ);
df_ext_heel_dx = blkdiag(df_ext_heel_diag{:})*[dJ2_q1(q);
                                            dJ2_q2(q);
                                            dJ2_q3(q);
                                            dJ2_q4(q);
                                            dJ2_q5(q);
                                            dJ2_q6(q)]; %dJ^T/dx*Fext
df_ext_heel_dx=[df_ext_heel_dx;zeros(6)];                                        
df_ext_toe_df = [zeros(p.numJ),J_toe/M];
df_ext_heel_df= [zeros(p.numJ),J_heel/M];
                                        
dfs = [[zeros(p.numJ);eye(p.numJ)],...
       (-dVdx-dGdx-dMdx+dtau_tend_1+dtau_tend_2+dtau_kne_tend_1+dtau_kne_tend_2+df_ext_toe_dx+df_ext_heel_dx+d_joint_fri )/M];

% grad =(dTaudx-dVdx-dGdx-dMdx+df_ext_toe1+df_ext_toe2+df_ext_heel1+df_ext_heel2-[zeros(p.numJ);p.joint_fri*eye(p.numJ);zeros(p.numJ);zeros(4,p.numJ)]+dtau_tend_1+dtau_tend_2+dtau_kne_tend_1+dtau_kne_tend_2)/M;

f = f.'; % we later wanna us it as column vector

% dfs = [[zeros(p.numJ);eye(p.numJ)],grad(1:p.numJ*2,:)]; %zeros(p.numJ) is because the output of f is [dq;ddq], 
dfu = [zeros(p.numJ),eye(p.numJ)/M];
df_fext = [df_ext_toe_df(1:2,:);df_ext_heel_df(1:2,:)];
 

end