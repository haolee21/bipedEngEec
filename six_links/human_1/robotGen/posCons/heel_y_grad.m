function heel_y_grad = heel_y_grad(in1)
%HEEL_Y_GRAD
%    HEEL_Y_GRAD = HEEL_Y_GRAD(IN1)

%    This function was generated by the Symbolic Math Toolbox version 8.4.
%    27-May-2020 00:26:09

q1 = in1(:,1);
q2 = in1(:,2);
q3 = in1(:,3);
q4 = in1(:,4);
q5 = in1(:,5);
t2 = cos(q1);
t3 = cos(q2);
t4 = cos(q3);
t5 = cos(q4);
t6 = cos(q5);
t7 = sin(q1);
t8 = sin(q2);
t9 = sin(q3);
t10 = sin(q4);
t11 = sin(q5);
t12 = t2.*t3;
t13 = t2.*t8;
t14 = t3.*t7;
t15 = t7.*t8;
t16 = -t15;
t17 = t13+t14;
t31 = t12.*4.38012e-1;
t32 = t15.*4.38012e-1;
t18 = t12+t16;
t19 = t4.*t17;
t20 = t9.*t17;
t33 = -t32;
t21 = t4.*t18;
t22 = t9.*t18;
t23 = -t20;
t24 = t19+t22;
t25 = t21+t23;
t28 = -t5.*(t20-t21);
t29 = -t10.*(t20-t21);
t36 = t5.*(t20-t21).*(-4.38012e-1);
t26 = t5.*t24;
t27 = t10.*t24;
t30 = -t27;
t34 = t27.*4.38012e-1;
t37 = t26+t29;
t41 = t6.*(t27+t5.*(t20-t21)).*(-5.2902e-1);
t35 = -t34;
t38 = t28+t30;
t39 = t11.*t37.*5.2902e-1;
t40 = -t39;
t42 = t35+t36+t40+t41;
heel_y_grad = [t2.*4.5252e-1+t31+t33+t42,t31+t33+t42,t42,t42,t40+t41,0.0];
