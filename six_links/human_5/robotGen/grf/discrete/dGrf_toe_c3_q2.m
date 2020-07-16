function dGrf_toe_c3_q2 = dGrf_toe_c3_q2(in1,in2,s)
%DGRF_TOE_C3_Q2
%    DGRF_TOE_C3_Q2 = DGRF_TOE_C3_Q2(IN1,IN2,S)

%    This function was generated by the Symbolic Math Toolbox version 8.4.
%    14-Jul-2020 22:57:53

q_t1 = in1(:,1);
q_t2 = in1(:,2);
q_t3 = in1(:,3);
q_t4 = in1(:,4);
q_t5 = in1(:,5);
q_t6 = in1(:,6);
t2 = cos(q_t1);
t3 = cos(q_t2);
t4 = cos(q_t3);
t5 = cos(q_t4);
t6 = cos(q_t5);
t7 = cos(q_t6);
t8 = sin(q_t1);
t9 = sin(q_t2);
t10 = sin(q_t3);
t11 = sin(q_t4);
t12 = sin(q_t5);
t13 = sin(q_t6);
t14 = t2.*t3;
t15 = t2.*t9;
t16 = t3.*t8;
t17 = t8.*t9;
t18 = -t17;
t19 = t15+t16;
t33 = t14.*4.453122e-1;
t34 = t17.*4.453122e-1;
t20 = t14+t18;
t21 = t4.*t19;
t22 = t10.*t19;
t35 = -t33;
t23 = t4.*t20;
t24 = t10.*t20;
t25 = -t22;
t26 = t21+t24;
t27 = t23+t25;
t30 = -t5.*(t22-t23);
t31 = -t11.*(t22-t23);
t37 = t5.*(t22-t23).*(-4.453122e-1);
t39 = t5.*(t22-t23).*4.453122e-1;
t28 = t5.*t26;
t29 = t11.*t26;
t32 = -t29;
t36 = t29.*4.453122e-1;
t38 = t28+t31;
t43 = -t6.*(t29+t5.*(t22-t23));
t44 = -t12.*(t29+t5.*(t22-t23));
t47 = t6.*(t29+t5.*(t22-t23)).*(-4.60062e-1);
t48 = t6.*(t29+t5.*(t22-t23)).*4.60062e-1;
t40 = t30+t32;
t41 = t6.*t38;
t42 = t12.*t38;
t45 = -t42;
t46 = t42.*4.60062e-1;
t49 = t41+t44;
t52 = t13.*(t42+t6.*(t29+t5.*(t22-t23))).*(-7.7775e-2);
t54 = t7.*(t42+t6.*(t29+t5.*(t22-t23))).*(-2.70108e-1);
t55 = t7.*(t42+t6.*(t29+t5.*(t22-t23))).*2.70108e-1;
t50 = t43+t45;
t51 = t7.*t49.*7.7775e-2;
t53 = t13.*t49.*2.70108e-1;
t56 = t36+t39+t46+t48+t51+t52+t53+t55;
t57 = (s.*t56)./2.0;
t58 = -t57;
dGrf_toe_c3_q2 = [s.*(t2.*(-4.60062e-1)+t34+t35+t56).*(-1.0./2.0);s.*(t34+t35+t56).*(-1.0./2.0);t58;t58;s.*(t46+t48+t51+t52+t53+t55).*(-1.0./2.0);s.*(t51+t52+t53+t55).*(-1.0./2.0)];
