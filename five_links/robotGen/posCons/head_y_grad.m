function head_y_grad = head_y_grad(in1)
%HEAD_Y_GRAD
%    HEAD_Y_GRAD = HEAD_Y_GRAD(IN1)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    21-Apr-2020 22:42:55

q1 = in1(:,1);
q2 = in1(:,2);
q3 = in1(:,3);
t2 = q1+q2;
t3 = cos(t2);
t4 = q3+t2;
t5 = cos(t4);
t6 = t3.*(9.0./2.0e+1);
t7 = t5.*(9.0./1.0e+1);
head_y_grad = [t6+t7+cos(q1).*(9.0./2.0e+1),t6+t7,0.0];
