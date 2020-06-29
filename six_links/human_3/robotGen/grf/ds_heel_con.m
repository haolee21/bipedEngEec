function ds_heel_con = ds_heel_con(in1,s_heel,th)
%DS_HEEL_CON
%    DS_HEEL_CON = DS_HEEL_CON(IN1,S_HEEL,TH)

%    This function was generated by the Symbolic Math Toolbox version 8.4.
%    24-Jun-2020 17:52:13

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
t7 = cos(q6);
t8 = sin(q1);
t9 = sin(q2);
t10 = sin(q3);
t11 = sin(q4);
t12 = sin(q5);
t13 = sin(q6);
t14 = -th;
t15 = t2.*t3;
t16 = t2.*t9;
t17 = t3.*t8;
t18 = t8.*t9;
t22 = t8.*4.5252e-1;
t19 = -t18;
t20 = t16+t17;
t35 = t15.*4.38012e-1;
t36 = t18.*4.38012e-1;
t21 = t15+t19;
t23 = t4.*t20;
t24 = t10.*t20;
t37 = -t36;
t25 = t4.*t21;
t26 = t10.*t21;
t27 = -t24;
t28 = t23+t26;
t29 = t25+t27;
t32 = -t5.*(t24-t25);
t33 = -t11.*(t24-t25);
t40 = t5.*(t24-t25).*(-4.38012e-1);
t41 = t5.*(t24-t25).*4.38012e-1;
t30 = t5.*t28;
t31 = t11.*t28;
t34 = -t31;
t38 = t31.*4.38012e-1;
t42 = t30+t33;
t46 = -t12.*(t31+t5.*(t24-t25));
t47 = -t6.*(t31+t5.*(t24-t25));
t48 = t6.*(t31+t5.*(t24-t25));
t39 = -t38;
t43 = t32+t34;
t44 = t12.*t42;
t45 = t6.*t42;
t51 = t48.*(-4.5252e-1);
t52 = t48.*4.5252e-1;
t49 = t44.*4.5252e-1;
t53 = t45+t46;
t54 = t44+t48;
t50 = -t49;
t55 = t7.*t53.*7.65e-2;
t57 = t13.*t54.*7.65e-2;
t56 = -t55;
t58 = -t57;
t59 = t38+t41+t49+t52+t55+t58;
t60 = s_heel.*t59;
t61 = -t60;
ds_heel_con = [s_heel.*(t2.*4.5252e-1+t35-t36-t59);-s_heel.*(-t35+t36+t59);t61;t61;-s_heel.*(t49+t52+t55+t58);-s_heel.*(t55+t58);0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;t14+t22+sin(q1+q2+q3+q4+q5).*4.5252e-1+sin(q1+q2+q3+q4).*4.38012e-1+sin(q1+q2).*4.38012e-1+1.535969075209524e+3.*cos(q1+q2+q3+q4+q5+q6-atan(3.472941176470588)).*1.8e-4;t14+t16.*4.38012e-1+t17.*4.38012e-1+t22+t30.*4.38012e-1+t45.*4.5252e-1-t12.*(t31+t5.*(t24-t25)).*4.5252e-1-t7.*t54.*7.65e-2-t13.*t53.*7.65e-2-t11.*(t24-t25).*4.38012e-1];
