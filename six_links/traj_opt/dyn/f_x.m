function [out_t,grad_t,grad_t_k] = f_x(x,p,i)
%% this inputs has the dimension numJoints*3 X 2
%  only the current and delay ones are included
%  th is the ypos height to have grf
%  ypos is pre-calculated since it is easy to use gpuarray to achieve 

%  i is the index,less than half the grf is on the toe, more than half the
%  grf is on ankle
numJ = p.numJ;
if(length(x)>numJ*3)
    x = [x(1:numJ*3);x(numJ*3+1:end)]; %later comes first, [second,first]
end

q1 = x(1,1:numJ);
dq1 = x(1,numJ+1:2*numJ);
u1 = x(1,numJ*2+1:3*numJ);

M = six_M(q1(2),q1(3),q1(4),q1(5),q1(6));
G1 = six_G(q1(1),q1(2),q1(3),q1(4),q1(5),q1(6));
V = six_V(q1(2),q1(3),q1(4),q1(5),q1(6),dq1(1),dq1(2),dq1(3),dq1(4),dq1(5),dq1(6));


% if(i<floor(p.gaitT/p.sampT/2))
% %     grf act on toe
%     Mext = @Mext_toe;
%     beta_grf = @(q)beta_grf_toe(q,p.toe_th);
%     
%     dMext_dx1=@dMext_toe_dx1;
%     dMext_dx2=@dMext_toe_dx2;
%     dMext_dx3=@dMext_toe_dx3;
%     dMext_dx4=@dMext_toe_dx4;
%     dMext_dx5=@dMext_toe_dx5;
%     dMext_dx6=@dMext_toe_dx6;
%     
%     dbeta_dx1=@(q)dbeta_toe_dx1(q,p.toe_th);
%     dbeta_dx2=@(q)dbeta_toe_dx2(q,p.toe_th);
%     dbeta_dx3=@(q)dbeta_toe_dx3(q,p.toe_th);
%     dbeta_dx4=@(q)dbeta_toe_dx4(q,p.toe_th);
%     dbeta_dx5=@(q)dbeta_toe_dx5(q,p.toe_th);
%     dbeta_dx6=@(q)dbeta_toe_dx6(q,p.toe_th);
% else
% %     grf act on ank
%     Mext = @Mext_heel;
%     beta_grf = @(q)beta_grf_heel(q,p.toe_th);
%     
%     dMext_dx1=@dMext_heel_dx1;
%     dMext_dx2=@dMext_heel_dx2;
%     dMext_dx3=@dMext_heel_dx3;
%     dMext_dx4=@dMext_heel_dx4;
%     dMext_dx5=@dMext_heel_dx5;
%     dMext_dx6=@dMext_heel_dx6;
%     
%     dbeta_dx1=@(q)dbeta_heel_dx1(q,p.toe_th);
%     dbeta_dx2=@(q)dbeta_heel_dx2(q,p.toe_th);
%     dbeta_dx3=@(q)dbeta_heel_dx3(q,p.toe_th);
%     dbeta_dx4=@(q)dbeta_heel_dx4(q,p.toe_th);
%     dbeta_dx5=@(q)dbeta_heel_dx5(q,p.toe_th);
%     dbeta_dx6=@(q)dbeta_heel_dx6(q,p.toe_th);
% end


% M_ext = p.floor_stiff*Mext(q1(1:numJ));
% beta1 = beta_grf(q1(1:numJ));

% calculate dMext_dx first since later we need to mod the value based on
% extF

% dMextdx1 = p.floor_stiff*dMext_dx1(q1(1:numJ));
% dMextdx2 = p.floor_stiff*dMext_dx2(q1(1:numJ));
% dMextdx3 = p.floor_stiff*dMext_dx3(q1(1:numJ));
% dMextdx4 = p.floor_stiff*dMext_dx4(q1(1:numJ));
% dMextdx5 = p.floor_stiff*dMext_dx5(q1(1:numJ));
% dMextdx6 = p.floor_stiff*dMext_dx6(q1(1:numJ));
% 
% 
% extF = M_ext*(u1-G1).';
% if(extF(2)>0)
%     extF=zeros(size(extF,1),size(extF,2));
%     M_ext = zeros(size(extF,1),numJ);
%     dMextdx1 = zeros(size(extF,1),numJ);
%     dMextdx2 = zeros(size(extF,1),numJ);
%     dMextdx3= zeros(size(extF,1),numJ);
%     dMextdx4= zeros(size(extF,1),numJ);
%     dMextdx5 = zeros(size(extF,1),numJ);
%     dMextdx6= zeros(size(extF,1),numJ);
% end
% dMext_dx = zeros(numJ*3,numJ);
% dMext_dx(1,:) = (u1-G1)*dMextdx1.'*beta1.';
% dMext_dx(2,:) = (u1-G1)*dMextdx2.'*beta1.';
% dMext_dx(3,:) = (u1-G1)*dMextdx3.'*beta1.';
% dMext_dx(4,:) = (u1-G1)*dMextdx4.'*beta1.';
% dMext_dx(5,:) = (u1-G1)*dMextdx5.'*beta1.';
% dMext_dx(6,:) = (u1-G1)*dMextdx6.'*beta1.';





grad_t_k=zeros(numJ*3,numJ);

% out_t = (u1-V-G1-extF.'*beta1.')/M;

% if(i<floor(p.gaitT/p.sampT/2))
%     [tau_ext,dTau_ext]=toe_grf(x,p);
%     
% else
%     % force act on the heel
%     [tau_ext,dTau_ext]=heel_grf(x,p);
% end

[tau_ext_toe,dTau_ext_toe,~,~]=toe_grf(x,p);
[tau_ext_heel,dTau_ext_heel,~,~]=heel_grf(x,p);
out_t = (u1-V-G1+tau_ext_toe+tau_ext_heel)/M;






dTaudx = [zeros(numJ);zeros(numJ);eye(numJ)];
% dTaudx = gpuArray(dTaudx);

dVdx = dV_dx(dq1(1),dq1(2),dq1(3),dq1(4),dq1(5),dq1(6),q1(2),q1(3),q1(4),q1(5),q1(6));
dGdx = dG_dx(q1(1),q1(2),q1(3),q1(4),q1(5),q1(6));


% dBeta1_dx = zeros(numJ*3,numJ);
% dBeta1_dx(1,:) = (u1-G1)*M_ext.'*dbeta_dx1(q1(1:numJ)).';
% dBeta1_dx(2,:) = (u1-G1)*M_ext.'*dbeta_dx2(q1(1:numJ)).';
% dBeta1_dx(3,:) = (u1-G1)*M_ext.'*dbeta_dx3(q1(1:numJ)).';
% dBeta1_dx(4,:) = (u1-G1)*M_ext.'*dbeta_dx4(q1(1:numJ)).';
% dBeta1_dx(5,:) = (u1-G1)*M_ext.'*dbeta_dx5(q1(1:numJ)).';
% dBeta1_dx(6,:) = (u1-G1)*M_ext.'*dbeta_dx6(q1(1:numJ)).';

% dFextdx = dFext([q1,dq1],p.toe_th);






% 
% if(size(x,1)>1)
%     q2 = x(2,1:numJ);
%     u2 = x(2,numJ*2+1:3*numJ);
%     G2 = six_G(q2(1),q2(2),q2(3),q2(4),q2(5),q2(6));
%     beta2 = beta_grf(q2(1),q2(2),q2(3),q2(4),q2(5),q2(6),p.toe_th);
%     out_t = out_t-(u2-G2)*beta2.'/M;
%     dGdx2 = dG_dx(q2(1),q2(2),q2(3),q2(4),q2(5),q2(6));
%     
%     dBeta2_dx = zeros(numJ*3,numJ);
%     dBeta2_dx(1,:) = (G2-u2)*dbeta_dx1(q1(1),q1(2),q1(3),q1(4),q1(5),q1(6),p.toe_th).';
%     dBeta2_dx(2,:) = (G2-u2)*dbeta_dx2(q1(1),q1(2),q1(3),q1(4),q1(5),q1(6),p.toe_th).';
%     dBeta2_dx(3,:) = (G2-u2)*dbeta_dx3(q1(1),q1(2),q1(3),q1(4),q1(5),q1(6),p.toe_th).';
%     dBeta2_dx(4,:) = (G2-u2)*dbeta_dx4(q1(1),q1(2),q1(3),q1(4),q1(5),q1(6),p.toe_th).';
%     dBeta2_dx(5,:) = (G2-u2)*dbeta_dx5(q1(1),q1(2),q1(3),q1(4),q1(5),q1(6),p.toe_th).';
%     dBeta2_dx(6,:) = (G2-u2)*dbeta_dx6(q1(1),q1(2),q1(3),q1(4),q1(5),q1(6),p.toe_th).';
% 
%     
%     
%     
%     grad_t_k = (-dBeta2_dx-dTaudx*beta2.'+dGdx2*beta2.')/M;
%     
% end



% dMdxBeta=gpuArray(dMdxBeta);


% dMdx2=zeros(numJ*3,numJ);
% dMdx2(2,:) = out_t*dM_dx2(q1(2),q1(3),q1(4),q1(5),q1(6)); % can just use row vec since M is symmetric
% dMdx2(3,:) = out_t*dM_dx3(q1(2),q1(3),q1(4),q1(5),q1(6));
% dMdx2(4,:) = out_t*dM_dx4(q1(2),q1(3),q1(4),q1(5),q1(6));
% dMdx2(5,:) = out_t*dM_dx5(q1(2),q1(3),q1(4),q1(5),q1(6));
% dMdx2(6,:) = out_t*dM_dx6(q1(2),q1(3),q1(4),q1(5),q1(6));
out_t_diag = repmat({out_t},1,3*p.numJ);
dMdx = blkdiag(out_t_diag{:})*[zeros(p.numJ);
                               dM_dx2(q1(2),q1(3),q1(4),q1(5),q1(6));
                               dM_dx3(q1(2),q1(3),q1(4),q1(5),q1(6));
                               dM_dx4(q1(2),q1(3),q1(4),q1(5),q1(6));
                               dM_dx5(q1(2),q1(3),q1(4),q1(5),q1(6));
                               dM_dx6(q1(2),q1(3),q1(4),q1(5),q1(6));
                               zeros(2*p.numJ*p.numJ,p.numJ)];



grad_t =(dTaudx-dVdx-dGdx-dMdx+dTau_ext_toe+dTau_ext_heel)/M;
% grad_t =(dTaudx-dVdx-dGdx-dMdx-dBeta1_dx-dMext_dx-(dTaudx-dGdx)*M_ext.'*beta1.')/M; %we already add sigma term in grad_t



% if(~isnan(x))
%     if(cond(M_ext.'*beta1.')>1e5)
%         test =2;
%     end
% end

%% add stiffness to the knee
% when q2<0.5 deg, spring activate
% when q5<0.5 deg, spring activate
out_t(1,2)=out_t(1,2)+sigma_knee(q1(2))*p.knee_stiff*(q1(2)); %-0.0087 = -0.5 deg
out_t(1,5)=out_t(1,5)+sigma_knee(q1(5))*p.knee_stiff*(q1(5));
grad_t(2,2) = grad_t(2,2)+dsigma_knee(q1(2))*p.knee_stiff*(q1(2))+sigma_knee(q1(2))*p.knee_stiff;
grad_t(5,5) = grad_t(5,5)+dsigma_knee(q1(5))*p.knee_stiff*(q1(5))+sigma_knee(q1(5))*p.knee_stiff;

%% add ankle tendon
out_t(1,1)=out_t(1,1)+p.ank_stiff*(pi/2-q1(1));
out_t(1,6)=out_t(1,6)-p.ank_stiff*(3*pi/2+q1(6));

grad_t(1,1)=grad_t(1,1)-p.ank_stiff;
grad_t(6,6)=grad_t(6,6)-p.ank_stiff;

end
