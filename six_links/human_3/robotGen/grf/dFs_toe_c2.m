function dFs_toe_c2 = dFs_toe_c2(in1,Fx_toe)
%DFS_TOE_C2
%    DFS_TOE_C2 = DFS_TOE_C2(IN1,FX_TOE)

%    This function was generated by the Symbolic Math Toolbox version 8.4.
%    25-Jun-2020 12:39:53

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
t4 = cos(t3);
t5 = sin(t3);
t6 = q3+q4+t3;
t12 = t2.*4.5252e-1;
t7 = cos(t6);
t8 = q5+t6;
t9 = sin(t6);
t16 = t4.*4.38012e-1;
t17 = t5.*4.38012e-1;
t10 = sin(t8);
t11 = cos(t8);
t18 = q6+t8-1.290439793566535;
t19 = t7.*4.38012e-1;
t20 = t9.*4.38012e-1;
t13 = t11.*4.5252e-1;
t14 = t10.*4.5252e-1;
t21 = cos(t18);
t22 = sin(t18);
t15 = -t13;
t23 = t21.*2.764744335377143e-1;
t24 = t22.*2.764744335377143e-1;
t27 = qd6.*t22.*(-2.764744335377143e-1);
t25 = qd6.*t24;
t26 = -t24;
t28 = t14+t23;
t29 = t15+t24;
t30 = -qd5.*(t13+t26);
t31 = qd5.*(t13+t26);
t32 = t20+t28;
t33 = t13+t19+t26;
t34 = Fx_toe.*t32;
t35 = qd1.*t33;
t36 = qd2.*t33;
t37 = qd3.*t33;
t38 = qd4.*t33;
t40 = t17+t32;
t41 = t16+t33;
t39 = -t34;
t42 = qd2.*t41;
t43 = t12+t40;
t44 = t27+t31+t35+t36+t37+t38;
t45 = Fx_toe.*t44;
t46 = -t45;
dFs_toe_c2 = [-Fx_toe.*(t27+t31+t37+t38+t42+qd1.*(t41+cos(q1).*4.5252e-1));-Fx_toe.*(t27+t31+t37+t38+t42+qd1.*t41);t46;t46;-Fx_toe.*(t27+t31+qd1.*(t13+t26)+qd2.*(t13+t26)+qd3.*(t13+t26)+qd4.*(t13+t26));Fx_toe.*(t25+qd1.*t24+qd2.*t24+qd3.*t24+qd4.*t24+qd5.*t24);-Fx_toe.*t43;-Fx_toe.*t40;t39;t39;-Fx_toe.*t28;Fx_toe.*t21.*(-2.764744335377143e-1);0.0;0.0;0.0;0.0;0.0;0.0;qd6.*t21.*(-2.764744335377143e-1)-qd5.*t28-qd3.*t32-qd4.*t32-qd2.*t40-qd1.*t43;0.0;0.0;0.0;0.0;0.0];
