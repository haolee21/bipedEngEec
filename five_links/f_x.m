function [out_t,grad_t,grad_t_k] = f_x(x,th)
%% this inputs has the dimension numJoints*3 X 2
%  only the current and delay ones are included
%  th is the ypos height to have grf
%  ypos is pre-calculated since it is easy to use gpuarray to achieve 

numJ = 5;
if(length(x)>numJ*3)
    x = [x(1:numJ*3);x(numJ*3+1:end)];
end

q1 = x(1,1:numJ);

M = five_M_simp(q1(2),q1(3),q1(4),q1(5));

ypos = end_y_pos([q1(1),q1(2),q1(3),q1(4),q1(5)]);

[out_beta_t,grad_beta_t] = betaFun(x(1,:));

if(ypos>=th || size(x,1)<2 )
    out_t = out_beta_t;
    grad_t = grad_beta_t;
    grad_t_k =gpuArray(zeros(size(x,2),numJ));
else
    M_t = five_M(q1(2),q1(3),q1(4),q1(5));
    [out_beta_t_k,grad_beta_t_k] = betaFun(x(2,:));
    out_t = out_beta_t - out_beta_t_k/M_t;
    
    
    % we need to solve the term dM/dx (M\beta/M), yet dM/dx is an tensor,
    % we name this term q1, and M\beta = q2, we cannot /M since the
    % dimension cannot work if we did not multiply tensor
    q2 = (M\out_beta_t.');
    q1 = gpuArray(zeros(3*numJ,numJ));
    q1(2,:) = dM_dx2(q1(2),q1(3),q1(4),q1(5))*q2;
    q1(3,:) = dM_dx3(q1(2),q1(3),q1(4),q1(5))*q2;
    q1(4,:) = dM_dx4(q1(2),q1(3),q1(4),q1(5))*q2;
    q1(5,:) = dM_dx5(q1(2),q1(3),q1(4),q1(5))*q2;
    
    grad_t = grad_beta_t +q1/M;
    
    % find the gradient
    
    grad_t_k = -grad_beta_t_k/M;
    
end

end
