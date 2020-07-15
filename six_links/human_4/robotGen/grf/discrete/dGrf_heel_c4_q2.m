function dGrf_heel_c4_q2 = dGrf_heel_c4_q2(in1,in2,Fx,sampT)
%DGRF_HEEL_C4_Q2
%    DGRF_HEEL_C4_Q2 = DGRF_HEEL_C4_Q2(IN1,IN2,FX,SAMPT)

%    This function was generated by the Symbolic Math Toolbox version 8.4.
%    14-Jul-2020 11:54:42

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
t2 = q_t1+q_t2;
t5 = 1.0./sampT;
t3 = cos(t2);
t4 = sin(t2);
t6 = q_t3+q_t4+t2;
t7 = cos(t6);
t8 = q_t5+t6;
t9 = sin(t6);
t23 = t3.*4.13678e-1;
t24 = t4.*4.13678e-1;
t10 = q_t6+t8;
t11 = sin(t8);
t12 = cos(t8);
t25 = t7.*4.13678e-1;
t26 = t9.*4.13678e-1;
t13 = cos(t10);
t14 = sin(t10);
t15 = t12.*4.2738e-1;
t16 = t11.*4.2738e-1;
t17 = -t15;
t18 = t13.*4.2738e-1;
t19 = t14.*4.2738e-1;
t22 = qd_t6.*t14.*(-4.2738e-1);
t20 = qd_t6.*t19;
t21 = -t19;
t27 = t17+t19;
t30 = t16+t18+t26;
t28 = qd_t5.*t27;
t31 = t15+t21+t25;
t36 = Fx.*t5.*t30;
t29 = -t28;
t32 = qd_t1.*t31;
t33 = qd_t2.*t31;
t34 = qd_t3.*t31;
t35 = qd_t4.*t31;
t37 = -t36;
t38 = t23+t31;
t39 = qd_t2.*t38;
t40 = t22+t29+t32+t33+t34+t35;
t41 = (Fx.*t40)./2.0;
t42 = -t41;
t43 = t37+t42;
dGrf_heel_c4_q2 = [Fx.*(t22+t29+t34+t35+t39+qd_t1.*(t38+cos(q_t1).*4.2738e-1)).*(-1.0./2.0)-Fx.*t5.*(t24+t30+sin(q_t1).*4.2738e-1);Fx.*(t22+t29+t34+t35+t39+qd_t1.*t38).*(-1.0./2.0)-Fx.*t5.*(t24+t30);t43;t43;(Fx.*(t20+t28+qd_t1.*t27+qd_t2.*t27+qd_t3.*t27+qd_t4.*t27))./2.0-Fx.*t5.*(t16+t18);(Fx.*(t20+qd_t1.*t19+qd_t2.*t19+qd_t3.*t19+qd_t4.*t19+qd_t5.*t19))./2.0-Fx.*t5.*t13.*4.2738e-1];
