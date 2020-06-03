function end_y_grad = end_y_grad(in1)
%END_Y_GRAD
%    END_Y_GRAD = END_Y_GRAD(IN1)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    24-Apr-2020 12:32:12

q1 = in1(:,1);
q2 = in1(:,2);
q3 = in1(:,3);
q4 = in1(:,4);
q5 = in1(:,5);
t2 = q1+q2;
t3 = cos(t2);
t4 = q3+q4+t2;
t5 = cos(t4);
t6 = q5+t4;
t8 = t3.*(2.0./5.0);
t7 = cos(t6);
t9 = t5.*(2.0./5.0);
t10 = t7.*(4.3e+1./1.0e+2);
t11 = t9+t10;
end_y_grad = [t8+t11+cos(q1).*(4.3e+1./1.0e+2),t8+t11,t11,t11,t10];
