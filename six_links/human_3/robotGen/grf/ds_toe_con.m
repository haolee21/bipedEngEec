function ds_toe_con = ds_toe_con(in1,s_toe,th)
%DS_TOE_CON
%    DS_TOE_CON = DS_TOE_CON(IN1,S_TOE,TH)

%    This function was generated by the Symbolic Math Toolbox version 8.4.
%    24-Jun-2020 17:51:36

q1 = in1(:,1);
q2 = in1(:,2);
q3 = in1(:,3);
q4 = in1(:,4);
q5 = in1(:,5);
q6 = in1(:,6);
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
t12 = q1+q2;
t14 = -th;
t24 = atan(3.472941176470588);
t25 = 1.535969075209524e+3;
t13 = cos(t12);
t15 = t2.*t3;
t16 = t2.*t8;
t17 = t3.*t7;
t18 = q3+q4+t12;
t19 = t7.*t8;
t27 = -t24;
t29 = t7.*4.5252e-1;
t20 = cos(t18);
t21 = q5+t18;
t22 = -t19;
t26 = t16+t17;
t45 = t13.*4.38012e-1;
t23 = cos(t21);
t28 = t15+t22;
t30 = t4.*t26;
t31 = t9.*t26;
t36 = q6+t21+t27;
t47 = t20.*4.38012e-1;
t32 = t4.*t28;
t33 = t9.*t28;
t34 = -t31;
t35 = t23.*4.5252e-1;
t37 = sin(t36);
t38 = t30+t33;
t39 = t32+t34;
t42 = -t5.*(t31-t32);
t43 = -t10.*(t31-t32);
t46 = t25.*t37.*1.8e-4;
t40 = t5.*t38;
t41 = t10.*t38;
t48 = -t46;
t44 = -t41;
t49 = t40+t43;
t51 = t35+t47+t48;
t50 = t42+t44;
t52 = s_toe.*t51;
ds_toe_con = [s_toe.*(t2.*4.5252e-1+t45+t51);s_toe.*(t45+t51);t52;t52;s_toe.*(t35+t48);s_toe.*t25.*t37.*(-1.8e-4);0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;t14+t29+sin(t12).*4.38012e-1+sin(t18).*4.38012e-1+sin(t21).*4.5252e-1+t25.*cos(t36).*1.8e-4;t14+t16.*4.38012e-1+t17.*4.38012e-1+t29+t40.*4.38012e-1+sin(q6).*(t11.*(t41+t5.*(t31-t32))-t6.*t49).*7.65e-2-t11.*(t41+t5.*(t31-t32)).*4.5252e-1+t6.*t49.*4.5252e-1-t10.*(t31-t32).*4.38012e-1-cos(q6).*(t6.*(t41+t5.*(t31-t32))+t11.*t49).*7.65e-2];
