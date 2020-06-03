function [out_t,grad_t,grad_t_k] = f_x(x,th)
%% this inputs has the dimension numJoints*3 X 2
%  only the current and delay ones are included
%  th is the ypos height to have grf
%  ypos is pre-calculated since it is easy to use gpuarray to achieve 

numJ = 5;
if(length(x)>numJ*3)
    x = [x(1:numJ*3);x(numJ*3+1:end)]; %later comes first, [second,first]
end

q1 = x(1,1:numJ);
dq1 = x(1,numJ+1:2*numJ);
u1 = x(1,numJ*2+1:3*numJ);

M = five_M(q1(2),q1(3),q1(4),q1(5));
G = five_G(q1(1),q1(2),q1(3),q1(4),q1(5));
V = five_V(q1(2),q1(3),q1(4),q1(5),dq1(1),dq1(2),dq1(3),dq1(4),dq1(5));


out_t = (M\(u1.'-V*dq1.'-G.')).';


dTaudx = [zeros(numJ);zeros(numJ);eye(numJ)];
% dTaudx = gpuArray(dTaudx);

dVdx = dV_dx(dq1(1),dq1(2),dq1(3),dq1(4),dq1(5),q1(2),q1(3),q1(4),q1(5));
dGdx = dG_dx(q1(1),q1(2),q1(3),q1(4),q1(5)).';
grad_t_k=[];
grad_t = zeros(numJ*3,numJ);
if(size(x,1)>1 )
    u2 = x(2,numJ*2+1:3*numJ);
    sigma = sigma_out(x(1,1),x(1,2),x(1,3),x(1,4),x(1,5),th);
    dsigma = [dsigma_dq(x(1,1),x(1,2),x(1,3),x(1,4),x(1,5),th);zeros(10,1)]; % I should use param here.....
    
    out_t = out_t - u2*sigma/M;
    grad_t = -dsigma*u2/M;
    
    % find the gradient
    
    grad_t_k = -[zeros(numJ*2,numJ);eye(numJ)]*sigma/M;
    
end

dMdxBeta=zeros(numJ*3,numJ);
% dMdxBeta=gpuArray(dMdxBeta);

dMdxBeta(2,:) = out_t*dM_dx2(q1(2),q1(3),q1(4),q1(5)); % can just use row vec since M is symmetric
dMdxBeta(3,:) = out_t*dM_dx3(q1(2),q1(3),q1(4),q1(5));
dMdxBeta(4,:) = out_t*dM_dx4(q1(2),q1(3),q1(4),q1(5));
dMdxBeta(5,:) = out_t*dM_dx5(q1(2),q1(3),q1(4),q1(5));

grad_t = grad_t + (dTaudx-dVdx-dGdx-dMdxBeta)/M; %we already add sigma term in grad_t
end
