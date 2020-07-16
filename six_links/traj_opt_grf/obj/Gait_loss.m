function loss = Gait_loss(x,p)

loss = (hip_x_pos(x(:,1).')-hip_x_pos(x(:,end).')-p.hipLen)^2;
end