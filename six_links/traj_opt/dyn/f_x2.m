function [f,dfs,dfu] = f_x2(x,p)


q = x(1:p.numJ,1).';
dq = x(p.numJ+1:2*p.numJ,1).';
u = x(p.numJ*2+1:3*p.numJ,1).';


M = six_M(q(2),q(3),q(4),q(5),q(6));
G = six_G(q(1),q(2),q(3),q(4),q(5),q(6));
V = six_V(q(2),q(3),q(4),q(5),q(6),dq(1),dq(2),dq(3),dq(4),dq(5),dq(6));

[tau_ext_toe,dTau_ext_toe,~,~]=toe_grf(x.',p);
[tau_ext_heel,dTau_ext_heel,~,~]=heel_grf(x.',p);


% ankle tendon
tau_tend_1 = [0,0,0,0,0,0];
tau_tend_2 = [0,0,0,0,0,0];
dtau_tend_1 = zeros(p.numJ*3,p.numJ);
dtau_tend_2 = zeros(p.numJ*3,p.numJ);

if(q(1)>pi/2)
    tau_tend_1 =  [p.ank_stiff*(pi/2-q(1)),0,0,0,0,0];
    dtau_tend_1(1,1) = -p.ank_stiff;
end
if(q(6)<-3*pi/2)
    tau_tend_2 = [0,0,0,0,0,p.ank_stiff*(-3*pi/2-q(6))];
    dtau_tend_2(6,6) = -p.ank_stiff;
end
    
% knee tendon

dtau_kne_tend_1 = zeros(p.numJ*3,p.numJ);
dtau_kne_tend_2 = zeros(p.numJ*3,p.numJ);
tau_kne_tend_1 = [0,-p.knee_stiff*q(2),0,0,0,0];
tau_kne_tend_2 = [0,0,0,0,-p.knee_stiff*q(5),0];
dtau_kne_tend_1(2,2) = -p.knee_stiff;
dtau_kne_tend_2(5,5) = -p.knee_stiff;



ddq = (u-V-G+tau_ext_toe+tau_ext_heel-dq*eye(p.numJ)*p.joint_fri+tau_tend_1+tau_tend_2+tau_kne_tend_1+tau_kne_tend_2)/M;

f = [dq,ddq];

dTaudx = [zeros(p.numJ);zeros(p.numJ);eye(p.numJ)];

dVdx = dV_dx(dq(1),dq(2),dq(3),dq(4),dq(5),dq(6),q(2),q(3),q(4),q(5),q(6));
dGdx = dG_dx(q(1),q(2),q(3),q(4),q(5),q(6));

out_t_diag = repmat({ddq},1,3*p.numJ);
dMdx = blkdiag(out_t_diag{:})*[zeros(p.numJ);
                               dM_dx2(q(2),q(3),q(4),q(5),q(6));
                               dM_dx3(q(2),q(3),q(4),q(5),q(6));
                               dM_dx4(q(2),q(3),q(4),q(5),q(6));
                               dM_dx5(q(2),q(3),q(4),q(5),q(6));
                               dM_dx6(q(2),q(3),q(4),q(5),q(6));
                               zeros(2*p.numJ*p.numJ,p.numJ)];



grad_ddq =(dTaudx-dVdx-dGdx-dMdx+dTau_ext_toe+dTau_ext_heel-[zeros(p.numJ);p.joint_fri*eye(p.numJ);zeros(p.numJ)]+dtau_tend_1+dtau_tend_2+dtau_kne_tend_1+dtau_kne_tend_2)/M;

f = f.'; % we later wanna us it as column vector

dfs = [[zeros(p.numJ);eye(p.numJ)],grad_ddq(1:p.numJ*2,:)];
dfu = [zeros(p.numJ),grad_ddq(p.numJ*2+1:end,:)];

 end