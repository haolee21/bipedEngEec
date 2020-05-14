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

if(i<floor(p.gaitT/p.sampT/2))
    % grf act on toe
    Mext = @Mext_toe;
    beta_grf = @beta_grf_toe;
    
    dMext_dx1=@dMext_toe_dx1;
    dMext_dx2=@dMext_toe_dx2;
    dMext_dx3=@dMext_toe_dx3;
    dMext_dx4=@dMext_toe_dx4;
    dMext_dx5=@dMext_toe_dx5;
    dMext_dx6=@dMext_toe_dx6;
    
    dbeta_dx1=@dbeta_toe_dx1;
    dbeta_dx2=@dbeta_toe_dx2;
    dbeta_dx3=@dbeta_toe_dx3;
    dbeta_dx4=@dbeta_toe_dx4;
    dbeta_dx5=@dbeta_toe_dx5;
    dbeta_dx6=@dbeta_toe_dx6;
else
    % grf act on ank
    Mext = @Mext_ank;
    beta_grf = @beta_grf_ank;
    
    dMext_dx1=@dMext_ank_dx1;
    dMext_dx2=@dMext_ank_dx2;
    dMext_dx3=@dMext_ank_dx3;
    dMext_dx4=@dMext_ank_dx4;
    dMext_dx5=@dMext_ank_dx5;
    dMext_dx6=@dMext_ank_dx6;
    
    dbeta_dx1=@dbeta_ank_dx1;
    dbeta_dx2=@dbeta_ank_dx2;
    dbeta_dx3=@dbeta_ank_dx3;
    dbeta_dx4=@dbeta_ank_dx4;
    dbeta_dx5=@dbeta_ank_dx5;
    dbeta_dx6=@dbeta_ank_dx6;
end


M_ext = Mext(q1(1:numJ));
beta1 = beta_grf(q1(1:numJ),p.toe_th);

% calculate dMext_dx first since later we need to mod the value based on
% extF

dMextdx1 = dMext_dx1(q1(1:numJ));
dMextdx2 = dMext_dx2(q1(1:numJ));
dMextdx3 = dMext_dx3(q1(1:numJ));
dMextdx4 = dMext_dx4(q1(1:numJ));
dMextdx5 = dMext_dx5(q1(1:numJ));
dMextdx6 = dMext_dx6(q1(1:numJ));


extF = M_ext*(u1-G1).';
if(extF(2)>0)
    extF=[0;0];
    M_ext = zeros(2,numJ);
    dMextdx1 = zeros(2,numJ);
    dMextdx2 = zeros(2,numJ);
    dMextdx3= zeros(2,numJ);
    dMextdx4= zeros(2,numJ);
    dMextdx5 = zeros(2,numJ);
    dMextdx6= zeros(2,numJ);
end
dMext_dx = zeros(numJ*3,numJ);
dMext_dx(1,:) = (u1-G1)*dMextdx1.'*beta1.';
dMext_dx(2,:) = (u1-G1)*dMextdx2.'*beta1.';
dMext_dx(3,:) = (u1-G1)*dMextdx3.'*beta1.';
dMext_dx(4,:) = (u1-G1)*dMextdx4.'*beta1.';
dMext_dx(5,:) = (u1-G1)*dMextdx5.'*beta1.';
dMext_dx(6,:) = (u1-G1)*dMextdx6.'*beta1.';





grad_t_k=zeros(numJ*3,numJ);

out_t = (u1-V-G1-extF.'*beta1.')/M;



dTaudx = [zeros(numJ);zeros(numJ);eye(numJ)];
% dTaudx = gpuArray(dTaudx);

dVdx = dV_dx(dq1(1),dq1(2),dq1(3),dq1(4),dq1(5),dq1(6),q1(2),q1(3),q1(4),q1(5),q1(6));
dGdx = dG_dx(q1(1),q1(2),q1(3),q1(4),q1(5),q1(6));


dBeta1_dx = zeros(numJ*3,numJ);
dBeta1_dx(1,:) = (u1-G1)*M_ext.'*dbeta_dx1(q1(1:numJ),p.toe_th).';
dBeta1_dx(2,:) = (u1-G1)*M_ext.'*dbeta_dx2(q1(1:numJ),p.toe_th).';
dBeta1_dx(3,:) = (u1-G1)*M_ext.'*dbeta_dx3(q1(1:numJ),p.toe_th).';
dBeta1_dx(4,:) = (u1-G1)*M_ext.'*dbeta_dx4(q1(1:numJ),p.toe_th).';
dBeta1_dx(5,:) = (u1-G1)*M_ext.'*dbeta_dx5(q1(1:numJ),p.toe_th).';
dBeta1_dx(6,:) = (u1-G1)*M_ext.'*dbeta_dx6(q1(1:numJ),p.toe_th).';








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


dMdx=zeros(numJ*3,numJ);
dMdx(2,:) = out_t*dM_dx2(q1(2),q1(3),q1(4),q1(5),q1(6)); % can just use row vec since M is symmetric
dMdx(3,:) = out_t*dM_dx3(q1(2),q1(3),q1(4),q1(5),q1(6));
dMdx(4,:) = out_t*dM_dx4(q1(2),q1(3),q1(4),q1(5),q1(6));
dMdx(5,:) = out_t*dM_dx5(q1(2),q1(3),q1(4),q1(5),q1(6));
dMdx(6,:) = out_t*dM_dx6(q1(2),q1(3),q1(4),q1(5),q1(6));

grad_t =(dTaudx-dVdx-dGdx-dMdx-dBeta1_dx-dMext_dx-(dTaudx-dGdx)*M_ext.'*beta1.')/M; %we already add sigma term in grad_t


end
