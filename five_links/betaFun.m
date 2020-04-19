function [out,grad]=betaFun(x)
%% the input of the function is only one row (one time step)
% the num of joints are hard coded since if we change it, everything need
% to be modify

numJ = 5;

q = x(1,1:numJ);
dq = x(1,numJ+1:numJ*2);
u = x(1,numJ*2+1:numJ*3);


M = five_M(q(2),q(3),q(4),q(5));
G = five_G(q(1),q(2),q(3),q(4),q(5));
V = five_V(q(2),q(3),q(4),q(5),dq(1),dq(2),dq(3),dq(4),dq(5));

out = (M\(u.'-V*dq.'-G.')).';

% calculate gradient

dTaudx = [zeros(numJ);zeros(numJ);eye(numJ)];
% dTaudx = gpuArray(dTaudx);

dVdx = dV_dx(dq(1),dq(2),dq(3),dq(4),dq(5),q(2),q(3),q(4),q(5));
dGdx = dG_dx(q(1),q(2),q(3),q(4),q(5)).';

dMdxBeta=zeros(numJ*3,numJ);
% dMdxBeta=gpuArray(dMdxBeta);

dMdxBeta(2,:) = out*dM_dx2(q(2),q(3),q(4),q(5)); % can just use row vec since M is symmetric
dMdxBeta(3,:) = out*dM_dx3(q(2),q(3),q(4),q(5));
dMdxBeta(4,:) = out*dM_dx4(q(2),q(3),q(4),q(5));
dMdxBeta(5,:) = out*dM_dx5(q(2),q(3),q(4),q(5));




grad = (dTaudx-dVdx-dGdx-dMdxBeta)/M;





end