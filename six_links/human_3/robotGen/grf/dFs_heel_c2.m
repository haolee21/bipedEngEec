function dFs_heel_c2 = dFs_heel_c2(in1,Fx_heel)
%DFS_HEEL_C2
%    DFS_HEEL_C2 = DFS_HEEL_C2(IN1,FX_HEEL)

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
t15 = t2.*4.5252e-1;
t7 = cos(t6);
t8 = q5+t6;
t9 = sin(t6);
t24 = t4.*4.38012e-1;
t25 = t5.*4.38012e-1;
t10 = q6+t8;
t11 = sin(t8);
t12 = cos(t8);
t27 = t7.*4.38012e-1;
t28 = t9.*4.38012e-1;
t13 = cos(t10);
t14 = sin(t10);
t16 = t12.*4.5252e-1;
t17 = t11.*4.5252e-1;
t18 = -t16;
t19 = t13.*4.5252e-1;
t20 = t14.*4.5252e-1;
t23 = qd6.*t14.*(-4.5252e-1);
t21 = qd6.*t20;
t22 = -t20;
t26 = t17+t19;
t29 = t18+t20;
t30 = qd5.*t29;
t32 = t26+t28;
t33 = t16+t22+t27;
t31 = -t30;
t34 = Fx_heel.*t32;
t35 = qd1.*t33;
t36 = qd2.*t33;
t37 = qd3.*t33;
t38 = qd4.*t33;
t40 = t25+t32;
t41 = t24+t33;
t39 = -t34;
t42 = qd2.*t41;
t43 = t15+t40;
t44 = t23+t31+t35+t36+t37+t38;
t45 = Fx_heel.*t44;
t46 = -t45;
dFs_heel_c2 = [-Fx_heel.*(t23+t31+t37+t38+t42+qd1.*(t41+cos(q1).*4.5252e-1));-Fx_heel.*(t23+t31+t37+t38+t42+qd1.*t41);t46;t46;Fx_heel.*(t21+t30+qd1.*t29+qd2.*t29+qd3.*t29+qd4.*t29);Fx_heel.*(t21+qd1.*t20+qd2.*t20+qd3.*t20+qd4.*t20+qd5.*t20);-Fx_heel.*t43;-Fx_heel.*t40;t39;t39;-Fx_heel.*t26;Fx_heel.*t13.*(-4.5252e-1);0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;qd6.*t13.*(-4.5252e-1)-qd5.*t26-qd3.*t32-qd4.*t32-qd2.*t40-qd1.*t43;0.0;0.0;0.0];
