function dFs_toe_c1_2 = dFs_toe_c1_2(in1,in2,us,ud)
%DFS_TOE_C1_2
%    DFS_TOE_C1_2 = DFS_TOE_C1_2(IN1,IN2,S_TOE,US,UD)

%    This function was generated by the Symbolic Math Toolbox version 8.4.
%    24-Jun-2020 13:53:47

q1 = in1(:,1);
q2 = in1(:,2);
q3 = in1(:,3);
q4 = in1(:,4);
q5 = in1(:,5);
q6 = in1(:,6);
qd1 = in1(:,7);
qd2 = in1(:,8);
qd3 = in1(:,9);
qd4 = in1(:,10);
qd5 = in1(:,11);
qd6 = in1(:,12);
t2 = sin(q1);
t3 = q1+q2;
t6 = -us;
t4 = cos(t3);
t5 = sin(t3);
t7 = q3+q4+t3;
t11 = t6+ud;
t14 = t2.*4.5252e-1;
t8 = cos(t7);
t9 = q5+t7;
t10 = sin(t7);
t18 = t4.*4.38012e-1;
t19 = t5.*4.38012e-1;
t12 = sin(t9);
t13 = cos(t9);
t20 = q6+t9-1.290439793566535;
t21 = t8.*4.38012e-1;
t22 = t10.*4.38012e-1;
t15 = t13.*4.5252e-1;
t16 = t12.*4.5252e-1;
t23 = cos(t20);
t24 = sin(t20);
t17 = -t15;
t25 = qd6.*t23.*2.764744335377143e+1;
t26 = qd6.*t24.*2.764744335377143e+1;
t28 = t23.*2.764744335377143e-1;
t29 = t24.*2.764744335377143e-1;
t27 = -t26;
t30 = -t29;
t31 = t16+t28;
t32 = t17+t29;
t33 = qd5.*t31.*1.0e+2;
t34 = qd5.*(t15+t30).*-1.0e+2;
t35 = qd5.*(t15+t30).*1.0e+2;
t36 = t22+t31;
t37 = t15+t21+t30;
t38 = qd3.*t36.*1.0e+2;
t39 = qd4.*t36.*1.0e+2;
t40 = qd3.*t37.*1.0e+2;
t41 = qd4.*t37.*1.0e+2;
t42 = t19+t36;
t43 = t18+t37;
t44 = qd2.*t42.*1.0e+2;
t45 = qd1.*t43.*1.0e+2;
t46 = qd2.*t43.*1.0e+2;
t47 = t14+t42;
t48 = qd1.*t47.*1.0e+2;
t49 = t27+t35+t40+t41+t45+t46;
t50 = t25+t33+t38+t39+t44+t48+5.0;
t51 = t25+t33+t38+t39+t44+t48-5.0;
t52 = tanh(t50);
t53 = tanh(t51);
dFs_toe_c1_2 = t11.*((t49.*(t52.^2-1.0))./2.0-(t49.*(t53.^2-1.0))./2.0).*(us+t11.*(t52.*(-1.0./2.0)+t53./2.0+1.0)).*-2.0;
