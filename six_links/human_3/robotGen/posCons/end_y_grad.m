function end_y_grad = end_y_grad(in1)
%END_Y_GRAD
%    END_Y_GRAD = END_Y_GRAD(IN1)

%    This function was generated by the Symbolic Math Toolbox version 8.4.
%    24-Jun-2020 17:48:29

q1 = in1(:,1);
q2 = in1(:,2);
q3 = in1(:,3);
q4 = in1(:,4);
q5 = in1(:,5);
q6 = in1(:,6);
t2 = q1+q2;
t8 = atan(3.472941176470588);
t9 = 1.535969075209524e+3;
t3 = cos(t2);
t4 = q3+q4+t2;
t10 = -t8;
t5 = cos(t4);
t6 = q5+t4;
t14 = t3.*4.38012e-1;
t7 = cos(t6);
t12 = q6+t6+t10;
t16 = t5.*4.38012e-1;
t11 = t7.*4.5252e-1;
t13 = sin(t12);
t15 = t9.*t13.*1.8e-4;
t17 = -t15;
t18 = t11+t16+t17;
end_y_grad = [t14+t18+cos(q1).*4.5252e-1,t14+t18,t18,t18,t11+t17,t17];
