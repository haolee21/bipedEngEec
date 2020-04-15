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

M = five_M_simp(q1(2),q1(3),q1(4),q1(5));



[out_beta_t,grad_beta_t] = betaFun(x(1,:));

if(size(x,1)<2 )
    out_t = out_beta_t;
    grad_t = grad_beta_t;
    grad_t_k =[];
else
    sigma = sigma_out(x(1,1),x(1,2),x(1,3),x(1,4),x(1,5),th);
    dsigma = [dsigma_dq(x(1,1),x(1,2),x(1,3),x(1,4),x(1,5),th);zeros(10,1)]; % I should use param here.....
    [out_beta_t_k,grad_beta_t_k] = betaFun(x(2,:));
    out_t = out_beta_t - out_beta_t_k/M*sigma;
    
    
    % we need to solve the term dM/dx (M\beta/M), yet dM/dx is an tensor,
    % we name this term q1, and M\beta = q2, we cannot /M since the
    % dimension cannot work if we did not multiply tensor
    q2 = M\out_beta_t.'*sigma;
    q1 = gpuArray(zeros(3*numJ,numJ));
    q1(2,:) = dM_dx2(q1(2),q1(3),q1(4),q1(5))*q2;
    q1(3,:) = dM_dx3(q1(2),q1(3),q1(4),q1(5))*q2;
    q1(4,:) = dM_dx4(q1(2),q1(3),q1(4),q1(5))*q2;
    q1(5,:) = dM_dx5(q1(2),q1(3),q1(4),q1(5))*q2;
    
    grad_t = grad_beta_t +(q1-dsigma*out_beta_t_k)/M;
    
    % find the gradient
    
    grad_t_k = -grad_beta_t_k*sigma/M;
    
end

end
