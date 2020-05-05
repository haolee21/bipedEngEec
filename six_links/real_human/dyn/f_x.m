function [out_t,grad_t,grad_t_k] = f_x(x,p)
%% this inputs has the dimension numJoints*3 X 2
%  only the current and delay ones are included
%  th is the ypos height to have grf
%  ypos is pre-calculated since it is easy to use gpuarray to achieve 

numJ = p.numJ;
if(length(x)>numJ*3)
    x = [x(1:numJ*3);x(numJ*3+1:end)]; %later comes first, [second,first]
end

q1 = x(1,1:numJ);
dq1 = x(1,numJ+1:2*numJ);
u1 = x(1,numJ*2+1:3*numJ);

M = six_M(q1(2),q1(3),q1(4),q1(5),q1(6));
G = six_G(q1(1),q1(2),q1(3),q1(4),q1(5),q1(6));
V = six_V(q1(2),q1(3),q1(4),q1(5),q1(6),dq1(1),dq1(2),dq1(3),dq1(4),dq1(5),dq1(6));
beta_out = beta_grf(q1(1),q1(2),q1(3),q1(4),q1(5),q1(6),p.toe_th);

out_t = (M\(u1.'-V.'-G.'-beta_out*(u1-G).')).';


dTaudx = [zeros(numJ);zeros(numJ);eye(numJ)];
% dTaudx = gpuArray(dTaudx);

dVdx = dV_dx(dq1(1),dq1(2),dq1(3),dq1(4),dq1(5),dq1(6),q1(2),q1(3),q1(4),q1(5),q1(6));
dGdx = dG_dx(q1(1),q1(2),q1(3),q1(4),q1(5),q1(6));
grad_t_k=zeros(numJ*3,numJ);

dBeta_dx = zeros(numJ*3,numJ);
dBeta_dx(1,:) = (G-u1)*dbeta_dx1(q1(1),q1(2),q1(3),q1(4),q1(5),q1(6),p.toe_th).';
dBeta_dx(2,:) = (G-u1)*dbeta_dx2(q1(1),q1(2),q1(3),q1(4),q1(5),q1(6),p.toe_th).';
dBeta_dx(3,:) = (G-u1)*dbeta_dx3(q1(1),q1(2),q1(3),q1(4),q1(5),q1(6),p.toe_th).';
dBeta_dx(4,:) = (G-u1)*dbeta_dx4(q1(1),q1(2),q1(3),q1(4),q1(5),q1(6),p.toe_th).';
dBeta_dx(5,:) = (G-u1)*dbeta_dx5(q1(1),q1(2),q1(3),q1(4),q1(5),q1(6),p.toe_th).';
dBeta_dx(6,:) = (G-u1)*dbeta_dx6(q1(1),q1(2),q1(3),q1(4),q1(5),q1(6),p.toe_th).';


dMdx=zeros(numJ*3,numJ);
% dMdxBeta=gpuArray(dMdxBeta);

dMdx(2,:) = out_t*dM_dx2(q1(2),q1(3),q1(4),q1(5),q1(6)); % can just use row vec since M is symmetric
dMdx(3,:) = out_t*dM_dx3(q1(2),q1(3),q1(4),q1(5),q1(6));
dMdx(4,:) = out_t*dM_dx4(q1(2),q1(3),q1(4),q1(5),q1(6));
dMdx(5,:) = out_t*dM_dx5(q1(2),q1(3),q1(4),q1(5),q1(6));
dMdx(6,:) = out_t*dM_dx6(q1(2),q1(3),q1(4),q1(5),q1(6));

grad_t =(dTaudx-dVdx-dGdx-dMdx+dBeta_dx+dGdx*beta_out.'-dTaudx*beta_out.')/M; %we already add sigma term in grad_t

end
