function [u_out,du_out_du] = newton_interp(u,ut)

% the input dimension of u will be [m,n], m is the dimension of each u, n
% is the total time

% u_out will be in dimension [m,1]
% du_out_du will be [m,n,m]
u_out = zeros(size(u,1),1);
du_out_du = zeros(size(u,1),size(u,2),size(u,1));
for i=1:size(u,2)
    coeff =1;
    for t_now =1:size(u,2)
        if t_now~=i
            
            coeff = coeff*(ut-t_now)/(i-t_now);
      
        end
    end
    u_out = u_out + coeff*u(:,i);
    du_out_du(:,i,:) = coeff*eye(size(u,1));
end



end