function dFn_heel_c2 = dFn_heel_c2(in1,in2,s_heel,th,k,cmax,dmax)
%DFN_HEEL_C2
%    DFN_HEEL_C2 = DFN_HEEL_C2(IN1,IN2,S_HEEL,TH,K,CMAX,DMAX)

%    This function was generated by the Symbolic Math Toolbox version 8.4.
%    24-Jun-2020 14:56:40

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
t14 = q1+q2;
t15 = th.*2.0;
t18 = 1.0./dmax;
t19 = -th;
t16 = cos(t14);
t17 = sin(t14);
t20 = -t15;
t21 = t2.*t3;
t22 = t2.*t9;
t23 = t3.*t8;
t24 = q3+q4+t14;
t25 = t8.*t9;
t37 = t2.*4.5252e-1;
t38 = t8.*9.0504e-1;
t39 = t8.*4.5252e-1;
t26 = cos(t24);
t27 = q5+t24;
t28 = sin(t24);
t31 = -t25;
t35 = t22+t23;
t61 = t16.*4.38012e-1;
t62 = t21.*4.38012e-1;
t63 = t21.*8.76024e-1;
t64 = t17.*4.38012e-1;
t65 = t22.*4.38012e-1;
t66 = t22.*8.76024e-1;
t67 = t23.*4.38012e-1;
t68 = t23.*8.76024e-1;
t69 = t25.*4.38012e-1;
t70 = t25.*8.76024e-1;
t29 = q6+t27;
t30 = sin(t27);
t32 = cos(t27);
t36 = t21+t31;
t40 = t4.*t35;
t41 = t10.*t35;
t72 = -t70;
t73 = -t69;
t74 = t26.*4.38012e-1;
t75 = t28.*4.38012e-1;
t33 = cos(t29);
t34 = sin(t29);
t42 = t4.*t36;
t43 = t10.*t36;
t44 = -t41;
t45 = t32.*4.5252e-1;
t46 = t30.*4.5252e-1;
t47 = -t45;
t48 = t33.*4.5252e-1;
t49 = t34.*4.5252e-1;
t53 = qd6.*t34.*(-4.5252e-1);
t54 = t40+t43;
t55 = t42+t44;
t58 = -t5.*(t41-t42);
t59 = -t11.*(t41-t42);
t86 = t5.*(t41-t42).*(-4.38012e-1);
t87 = t5.*(t41-t42).*(-8.76024e-1);
t88 = t11.*(t41-t42).*(-4.38012e-1);
t89 = t11.*(t41-t42).*(-8.76024e-1);
t90 = t5.*(t41-t42).*8.76024e-1;
t91 = t5.*(t41-t42).*4.38012e-1;
t50 = qd6.*t48;
t51 = qd6.*t49;
t52 = -t49;
t56 = t5.*t54;
t57 = t11.*t54;
t71 = t46+t48;
t76 = t47+t49;
t60 = -t57;
t77 = qd5.*t71;
t78 = qd5.*t76;
t80 = t56.*4.38012e-1;
t81 = t56.*8.76024e-1;
t82 = t57.*4.38012e-1;
t83 = t57.*8.76024e-1;
t92 = t56+t59;
t96 = -t12.*(t57+t5.*(t41-t42));
t97 = -t6.*(t57+t5.*(t41-t42));
t98 = t6.*(t57+t5.*(t41-t42));
t99 = t71+t75;
t100 = t45+t52+t74;
t113 = t12.*(t57+t5.*(t41-t42)).*(-9.0504e-1);
t114 = t12.*(t57+t5.*(t41-t42)).*(-4.5252e-1);
t79 = -t78;
t84 = -t83;
t85 = -t82;
t93 = t58+t60;
t94 = t12.*t92;
t95 = t6.*t92;
t101 = qd1.*t99;
t102 = qd2.*t99;
t103 = qd3.*t99;
t104 = qd4.*t99;
t105 = qd3.*t100;
t106 = qd4.*t100;
t115 = t98.*(-9.0504e-1);
t116 = t98.*(-4.5252e-1);
t117 = t98.*9.0504e-1;
t118 = t98.*4.5252e-1;
t119 = t64+t99;
t120 = t61+t100;
t107 = t94.*9.0504e-1;
t108 = t94.*4.5252e-1;
t109 = t95.*9.0504e-1;
t110 = t95.*4.5252e-1;
t121 = qd2.*t119;
t122 = qd2.*t120;
t123 = t37+t120;
t125 = t95+t96;
t126 = t94+t98;
t143 = t50+t77+t101+t102+t103+t104;
t111 = -t107;
t112 = -t108;
t124 = qd1.*t123;
t127 = t7.*t125.*(1.53e+2./1.0e+3);
t128 = t7.*t125.*7.65e-2;
t129 = t13.*t125.*(1.53e+2./1.0e+3);
t130 = t13.*t125.*7.65e-2;
t135 = t7.*t126.*(1.53e+2./1.0e+3);
t136 = t7.*t126.*7.65e-2;
t137 = t13.*t126.*(1.53e+2./1.0e+3);
t138 = t13.*t126.*7.65e-2;
t131 = -t127;
t132 = -t128;
t133 = -t129;
t134 = -t130;
t139 = -t135;
t140 = -t136;
t141 = -t137;
t142 = -t138;
t144 = t53+t79+t105+t106+t122+t124;
t145 = t83+t90+t107+t117+t127+t141;
t146 = t82+t91+t108+t118+t128+t142;
t147 = t20+t38+t66+t68+t81+t89+t109+t113+t133+t139;
t148 = t19+t39+t65+t67+t80+t88+t110+t114+t134+t140;
t149 = t18.*t147;
t158 = k.*t146.*t148.*2.0;
t150 = tanh(t149);
t151 = t150.^2;
t152 = t150./2.0;
t153 = t151-1.0;
t154 = t152-1.0./2.0;
t155 = cmax.*s_heel.*t100.*t154;
t157 = cmax.*t143.*t154;
t159 = (cmax.*t18.*t144.*t145.*t153)./2.0;
t156 = -t155;
t160 = -t159;
t161 = t157+t158+t160;
t162 = s_heel.*t161;
dFn_heel_c2 = [s_heel.*(k.*t148.*(t37+t62-t69-t146).*-2.0+cmax.*t154.*(t50+t77+t103+t104+t121+qd1.*(t39+t119))+(cmax.*t18.*t144.*t153.*(t2.*9.0504e-1+t63-t70-t145))./2.0);s_heel.*(cmax.*t154.*(t50+t77+t103+t104+t121+qd1.*t119)+k.*t148.*(-t62+t69+t146).*2.0-(cmax.*t18.*t144.*t153.*(-t63+t70+t145))./2.0);t162;t162;s_heel.*(k.*t148.*(t108+t118+t128+t142).*2.0+cmax.*t154.*(t50+t77+qd1.*t71+qd2.*t71+qd3.*t71+qd4.*t71)-(cmax.*t18.*t144.*t153.*(t107+t117+t127+t141))./2.0);s_heel.*(cmax.*t154.*(t50+qd1.*t48+qd2.*t48+qd3.*t48+qd4.*t48+qd5.*t48)+k.*t148.*(t128+t142).*2.0-(cmax.*t18.*t144.*t153.*(t127+t141))./2.0);-cmax.*s_heel.*t123.*t154;-cmax.*s_heel.*t120.*t154;t156;t156;cmax.*s_heel.*t76.*t154;cmax.*s_heel.*t49.*t154;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;1.0;0.0;-k.*t148.^2-cmax.*t144.*t154];
