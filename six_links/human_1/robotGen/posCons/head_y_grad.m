function head_y_grad = head_y_grad(in1)
%HEAD_Y_GRAD
%    HEAD_Y_GRAD = HEAD_Y_GRAD(IN1)

%    This function was generated by the Symbolic Math Toolbox version 8.4.
%    27-May-2020 00:26:10

q1 = in1(:,1);
q2 = in1(:,2);
q3 = in1(:,3);
t2 = q1+q2;
t3 = cos(t2);
t4 = q3+t2;
t5 = cos(t4);
t7 = t3.*4.38012e-1;
t6 = t5.*7.6014e-1;
head_y_grad = [t6+t7+cos(q1).*4.5252e-1,t6+t7,t6];
