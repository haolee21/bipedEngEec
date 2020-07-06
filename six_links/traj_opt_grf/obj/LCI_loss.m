function loss = LCI_loss(x,p)
loss = 0;
dq = gradient(x(1:p.numJ,:))/p.sampT;
for i=1:size(x,2)
    q = x(1:p.numJ,i);
    c_toe = x(p.numJ+1,i);
    c_heel = x(p.numJ+2,i);
    toe_pos = end_y_pos(q.');
    heel_pos = heel_y_pos(q.');
    
    v_toe = six_J(q(1),q(2),q(3),q(4),q(5),q(6))*dq(:,i);
    v_heel = six_J2(q(1),q(2),q(3),q(4),q(5),q(6))*dq(:,i);
    
    toe_x_vel = v_toe(1);
    heel_x_vel = v_heel(1);
    
    loss = loss + c_toe*((toe_pos-p.floor_h)^2+toe_x_vel^2)*1e-3;
    loss = loss + c_heel*((heel_pos-p.floor_h)^2+heel_x_vel^2)*1e-3;
    
end

end