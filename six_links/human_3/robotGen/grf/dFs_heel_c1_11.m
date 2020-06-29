function dFs_heel_c1_11 = dFs_heel_c1_11(in1,in2,us,ud)
%DFS_HEEL_C1_11
%    DFS_HEEL_C1_11 = DFS_HEEL_C1_11(IN1,IN2,S_HEEL,US,UD)

%    This function was generated by the Symbolic Math Toolbox version 8.4.
%    24-Jun-2020 13:54:02

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
t5 = -us;
t4 = sin(t3);
t6 = q3+q4+t3;
t9 = t5+ud;
t13 = t2.*4.5252e-1;
t7 = q5+t6;
t8 = sin(t6);
t20 = t4.*4.38012e-1;
t10 = q6+t7;
t11 = sin(t7);
t22 = t8.*4.38012e-1;
t12 = cos(t10);
t14 = t11.*4.5252e+1;
t17 = t11.*4.5252e-1;
t15 = t12.*4.5252e+1;
t18 = t12.*4.5252e-1;
t16 = qd6.*t15;
t19 = t14+t15;
t21 = t17+t18;
t23 = qd5.*t21.*1.0e+2;
t24 = t21+t22;
t25 = qd3.*t24.*1.0e+2;
t26 = qd4.*t24.*1.0e+2;
t27 = t20+t24;
t28 = qd2.*t27.*1.0e+2;
t29 = t13+t27;
t30 = qd1.*t29.*1.0e+2;
t31 = t16+t23+t25+t26+t28+t30+5.0;
t33 = t16+t23+t25+t26+t28+t30-5.0;
t32 = tanh(t31);
t34 = tanh(t33);
dFs_heel_c1_11 = t9.*((t19.*(t32.^2-1.0))./2.0-(t19.*(t34.^2-1.0))./2.0).*(us+t9.*(t32.*(-1.0./2.0)+t34./2.0+1.0)).*-2.0;
