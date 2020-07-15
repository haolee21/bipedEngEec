function out1 = dL2_5(in1,in2,sampT)
%DL2_5
%    OUT1 = DL2_5(IN1,IN2,SAMPT)

%    This function was generated by the Symbolic Math Toolbox version 8.4.
%    14-Jul-2020 11:51:04

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
t2 = cos(q_t1);
t3 = cos(q_t2);
t4 = cos(q_t3);
t5 = cos(q_t4);
t6 = cos(q_t5);
t7 = sin(q_t1);
t8 = sin(q_t2);
t9 = sin(q_t3);
t10 = sin(q_t4);
t11 = sin(q_t5);
t12 = q_t1+q_t2;
t13 = cos(t12);
t14 = sin(t12);
t15 = t2.*t3;
t16 = t2.*t8;
t17 = t3.*t7;
t18 = q_t3+q_t4+t12;
t19 = t7.*t8;
t36 = qd_t1.*t2.*4.811871021864006e+19;
t37 = qd_t1.*t7.*4.811871021864006e+19;
t38 = qd_t1.*t2.*9.623742043728013e+19;
t39 = qd_t1.*t7.*9.623742043728013e+19;
t20 = cos(t18);
t21 = q_t5+t18;
t22 = sin(t18);
t25 = -t19;
t29 = t16+t17;
t40 = qd_t1.*t13.*4.65760021662843e+19;
t41 = qd_t2.*t13.*4.65760021662843e+19;
t42 = qd_t1.*t14.*4.65760021662843e+19;
t43 = qd_t2.*t14.*4.65760021662843e+19;
t44 = qd_t1.*t13.*9.31520043325686e+19;
t45 = qd_t2.*t13.*9.31520043325686e+19;
t46 = qd_t1.*t14.*9.31520043325686e+19;
t47 = qd_t2.*t14.*9.31520043325686e+19;
t23 = q_t6+t21;
t24 = sin(t21);
t26 = cos(t21);
t30 = t15+t25;
t31 = t4.*t29;
t32 = t9.*t29;
t48 = qd_t1.*t20.*4.65760021662843e+19;
t49 = qd_t2.*t20.*4.65760021662843e+19;
t50 = qd_t3.*t20.*4.65760021662843e+19;
t51 = qd_t4.*t20.*4.65760021662843e+19;
t52 = qd_t1.*t22.*4.65760021662843e+19;
t53 = qd_t2.*t22.*4.65760021662843e+19;
t54 = qd_t3.*t22.*4.65760021662843e+19;
t55 = qd_t4.*t22.*4.65760021662843e+19;
t66 = qd_t1.*t20.*9.31520043325686e+19;
t67 = qd_t2.*t20.*9.31520043325686e+19;
t68 = qd_t3.*t20.*9.31520043325686e+19;
t69 = qd_t4.*t20.*9.31520043325686e+19;
t70 = qd_t1.*t22.*9.31520043325686e+19;
t71 = qd_t2.*t22.*9.31520043325686e+19;
t72 = qd_t3.*t22.*9.31520043325686e+19;
t73 = qd_t4.*t22.*9.31520043325686e+19;
t27 = cos(t23);
t28 = sin(t23);
t33 = t4.*t30;
t34 = t9.*t30;
t35 = -t32;
t56 = qd_t1.*t26.*2.083540152467115e+19;
t57 = qd_t2.*t26.*2.083540152467115e+19;
t58 = qd_t3.*t26.*2.083540152467115e+19;
t59 = qd_t4.*t26.*2.083540152467115e+19;
t60 = qd_t5.*t26.*2.083540152467115e+19;
t61 = qd_t1.*t24.*2.083540152467115e+19;
t62 = qd_t2.*t24.*2.083540152467115e+19;
t63 = qd_t3.*t24.*2.083540152467115e+19;
t64 = qd_t4.*t24.*2.083540152467115e+19;
t65 = qd_t5.*t24.*2.083540152467115e+19;
t86 = qd_t1.*t26.*9.623742043728013e+19;
t87 = qd_t2.*t26.*9.623742043728013e+19;
t88 = qd_t3.*t26.*9.623742043728013e+19;
t89 = qd_t4.*t26.*9.623742043728013e+19;
t90 = qd_t5.*t26.*9.623742043728013e+19;
t91 = qd_t1.*t24.*9.623742043728013e+19;
t92 = qd_t2.*t24.*9.623742043728013e+19;
t93 = qd_t3.*t24.*9.623742043728013e+19;
t94 = qd_t4.*t24.*9.623742043728013e+19;
t95 = qd_t5.*t24.*9.623742043728013e+19;
t74 = qd_t1.*t27.*2.49457040483832e+19;
t75 = qd_t2.*t27.*2.49457040483832e+19;
t76 = qd_t3.*t27.*2.49457040483832e+19;
t77 = qd_t4.*t27.*2.49457040483832e+19;
t78 = qd_t5.*t27.*2.49457040483832e+19;
t79 = qd_t6.*t27.*2.49457040483832e+19;
t80 = qd_t1.*t28.*2.49457040483832e+19;
t81 = qd_t2.*t28.*2.49457040483832e+19;
t82 = qd_t3.*t28.*2.49457040483832e+19;
t83 = qd_t4.*t28.*2.49457040483832e+19;
t84 = qd_t5.*t28.*2.49457040483832e+19;
t85 = qd_t6.*t28.*2.49457040483832e+19;
t96 = t31+t34;
t97 = t33+t35;
t100 = -t5.*(t32-t33);
t101 = -t10.*(t32-t33);
t105 = t37+t42+t43+t52+t53+t54+t55+t61+t62+t63+t64+t65;
t106 = t36+t40+t41+t48+t49+t50+t51+t56+t57+t58+t59+t60;
t98 = t5.*t96;
t99 = t10.*t96;
t107 = t38+t44+t45+t66+t67+t68+t69+t74+t75+t76+t77+t78+t79+t86+t87+t88+t89+t90;
t108 = t39+t46+t47+t70+t71+t72+t73+t80+t81+t82+t83+t84+t85+t91+t92+t93+t94+t95;
t102 = -t99;
t103 = t98+t101;
t104 = t100+t102;
out1 = sin(q_t6).*(t11.*(t99+t5.*(t32-t33))-t6.*t103).*(-4.615484057022602e-1)+t6.*(t99+t5.*(t32-t33)).*4.21740078932502+t108.*(t74+t75+t76+t77+t78+t79+t86+t87+t88+t89+t90).*8.375730661184093e-42-t107.*(t80+t81+t82+t83+t84+t85+t91+t92+t93+t94+t95).*8.375730661184093e-42+t11.*t103.*4.21740078932502+t105.*(t56+t57+t58+t59+t60).*1.058887993078164e-40-t106.*(t61+t62+t63+t64+t65).*1.058887993078164e-40+(qd_t1.*4.472206631138927e-2+qd_t2.*4.472206631138927e-2+qd_t3.*4.472206631138927e-2+qd_t4.*4.472206631138927e-2+qd_t5.*4.472206631138927e-2+t24.*t105.*4.41247130108735e-21+t26.*t106.*4.41247130108735e-21+t107.*(t26.*9.623742043728013e+19+t27.*2.49457040483832e+19).*1.675146132236819e-41+t108.*(t24.*9.623742043728013e+19+t28.*2.49457040483832e+19).*1.675146132236819e-41)./sampT+cos(q_t6).*(t6.*(t99+t5.*(t32-t33))+t11.*t103).*4.615484057022602e-1;
