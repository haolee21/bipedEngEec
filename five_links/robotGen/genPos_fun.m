function genPos_fun()
syms q1 q2 q3 q4 q5
end_pos = endPos(q1,q2,q3,q4,q5);
ypos = end_pos(2);
matlabFunction(ypos,'file','end_y_pos');
ypos_grad=[diff(ypos,q1),diff(ypos,q2),diff(ypos,q3),diff(ypos,q4),diff(ypos,q5)];
matlabFunction(ypos_grad,'file','end_y_poss_grad');
% xpos = end_pos(1);
% matlabFunction(xpos,'file','end_x_pos');
% xpos_grad=[diff(xpos,q1),diff(xpos,q2),diff(xpos,q3),diff(xpos,q4),diff(xpos,q5)];
% matlabFunction(xpos_grad,'file','end_x_poss_grad');

end