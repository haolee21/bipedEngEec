function dGrf_toe_c5_q1 = dGrf_toe_c5_q1(in1,in2,s,sampT)
%DGRF_TOE_C5_Q1
%    DGRF_TOE_C5_Q1 = DGRF_TOE_C5_Q1(IN1,IN2,S,SAMPT)

%    This function was generated by the Symbolic Math Toolbox version 8.4.
%    14-Jul-2020 11:54:57

q_t1 = in1(:,1);
q_t2 = in1(:,2);
q_t3 = in1(:,3);
q_t4 = in1(:,4);
q_t5 = in1(:,5);
q_t6 = in1(:,6);
qd_t1 = in2(:,1);
qd_t2 = in2(:,2);
qd_t3 = in2(:,3);
qd_t4 = in2(:,4);
qd_t5 = in2(:,5);
qd_t6 = in2(:,6);
t2 = sin(q_t1);
t3 = q_t1+q_t2;
t6 = 1.0./sampT;
t13 = atan(3.472941176470588);
t14 = 1.535969075209524e+3;
t4 = cos(t3);
t5 = sin(t3);
t7 = q_t3+q_t4+t3;
t15 = -t13;
t16 = t2.*4.2738e-1;
t8 = cos(t7);
t9 = q_t5+t7;
t10 = sin(t7);
t22 = t4.*4.13678e-1;
t23 = t5.*4.13678e-1;
t11 = sin(t9);
t12 = cos(t9);
t19 = q_t6+t9+t15;
t24 = t8.*4.13678e-1;
t25 = t10.*4.13678e-1;
t17 = t12.*4.2738e-1;
t18 = t11.*4.2738e-1;
t20 = cos(t19);
t21 = sin(t19);
t26 = t14.*t20.*1.7e-4;
t27 = t14.*t21.*1.7e-4;
t31 = qd_t6.*t14.*t21.*(-1.7e-4);
t28 = qd_t6.*t26;
t29 = qd_t6.*t27;
t30 = -t27;
t32 = t18+t26;
t33 = t17+t30;
t34 = qd_t5.*t32;
t36 = t25+t32;
t35 = qd_t5.*t33;
t37 = t24+t33;
t38 = qd_t3.*t36;
t39 = qd_t4.*t36;
t44 = t23+t36;
t40 = qd_t1.*t37;
t41 = qd_t2.*t37;
t42 = qd_t3.*t37;
t43 = qd_t4.*t37;
t45 = t22+t37;
t46 = qd_t2.*t44;
t48 = t16+t44;
t47 = qd_t2.*t45;
t49 = qd_t1.*t48;
t50 = t31+t35+t40+t41+t42+t43;
t51 = t28+t34+t38+t39+t46+t49;
t52 = s.*t6.*t36.*t51.*2.0;
t54 = s.*t50.*t51;
t53 = -t52;
t55 = t53+t54;
dGrf_toe_c5_q1 = [s.*t51.*(t31+t35+t42+t43+t47+qd_t1.*(t45+cos(q_t1).*4.2738e-1))-s.*t6.*t48.*t51.*2.0;s.*t51.*(t31+t35+t42+t43+t47+qd_t1.*t45)-s.*t6.*t44.*t51.*2.0;t55;t55;s.*t51.*(t31+t35+qd_t1.*t33+qd_t2.*t33+qd_t3.*t33+qd_t4.*t33)-s.*t6.*t32.*t51.*2.0;-s.*t51.*(t29+qd_t1.*t27+qd_t2.*t27+qd_t3.*t27+qd_t4.*t27+qd_t5.*t27)-s.*t6.*t14.*t20.*t51.*3.4e-4];
