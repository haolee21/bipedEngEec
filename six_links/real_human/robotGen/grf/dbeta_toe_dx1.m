function dbeta_toe_dx1 = dbeta_toe_dx1(in1,th)
%DBETA_TOE_DX1
%    DBETA_TOE_DX1 = DBETA_TOE_DX1(IN1,TH)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    17-May-2020 03:22:58

q1 = in1(:,1);
q2 = in1(:,2);
q3 = in1(:,3);
q4 = in1(:,4);
q5 = in1(:,5);
q6 = in1(:,6);
t2 = cos(q1);
t3 = sin(q1);
t4 = q1+q2;
t8 = atan(7.0./2.6e+1);
t9 = atan(2.6e+1./7.0);
t10 = th.*4.0e+2;
t13 = sqrt(2.9e+1);
t5 = cos(t4);
t6 = sin(t4);
t7 = q3+q4+t4;
t11 = t2.*1.72e+2;
t12 = t3.*1.72e+2;
t17 = -t10;
t24 = -t9;
t25 = t2.*(4.3e+1./1.0e+2);
t26 = t3.*(4.3e+1./1.0e+2);
t14 = cos(t7);
t15 = q5+t7;
t16 = sin(t7);
t20 = t5.*(2.0./5.0);
t21 = t5.*1.6e+2;
t22 = t6.*(2.0./5.0);
t23 = t6.*1.6e+2;
t18 = sin(t15);
t19 = cos(t15);
t27 = t14.*(2.0./5.0);
t28 = t14.*1.6e+2;
t29 = t16.*(2.0./5.0);
t30 = t16.*1.6e+2;
t33 = q6+t8+t15;
t38 = q6+t15+t24;
t31 = t19.*1.72e+2;
t32 = t18.*1.72e+2;
t34 = cos(t33);
t35 = sin(t33);
t36 = t19.*(4.3e+1./1.0e+2);
t37 = t18.*(4.3e+1./1.0e+2);
t39 = t13.*t34.*2.0e+1;
t40 = t13.*t35.*2.0e+1;
t41 = (t13.*t34)./2.0e+1;
t42 = (t13.*t35)./2.0e+1;
t43 = t36+t41;
t44 = t37+t42;
t49 = t11+t21+t28+t31+t39;
t50 = t12+t17+t23+t30+t32+t40;
t45 = t27+t43;
t46 = t29+t44;
t51 = tanh(t50);
t47 = t20+t45;
t48 = t22+t46;
t54 = t51.^2;
t55 = t51./2.0;
t52 = t25+t47;
t53 = t26+t48;
t56 = t54-1.0;
t57 = t55-1.0./2.0;
t58 = t45.*t57;
t59 = t46.*t57;
t60 = (t45.*t49.*t56)./2.0;
t61 = (t46.*t49.*t56)./2.0;
t62 = -t61;
t63 = t59+t60;
t64 = t58+t62;
dbeta_toe_dx1 = reshape([t52.*t57-(t49.*t53.*t56)./2.0,t47.*t57-(t48.*t49.*t56)./2.0,t64,t64,t43.*t57-(t44.*t49.*t56)./2.0,t13.*t57.*sin(t38).*(-1.0./2.0e+1)-(t13.*t49.*t56.*cos(t38))./4.0e+1,t53.*t57+(t49.*t52.*t56)./2.0,t48.*t57+(t47.*t49.*t56)./2.0,t63,t63,t44.*t57+(t43.*t49.*t56)./2.0,t42.*t57+(t13.*t34.*t49.*t56)./4.0e+1],[6,2]);
