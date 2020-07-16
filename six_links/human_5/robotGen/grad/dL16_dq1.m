function dL16_dq1 = dL16_dq1(in1,in2,sampT)
%DL16_DQ1
%    DL16_DQ1 = DL16_DQ1(IN1,IN2,SAMPT)

%    This function was generated by the Symbolic Math Toolbox version 8.4.
%    14-Jul-2020 22:50:35

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
t7 = cos(q_t6);
t8 = sin(q_t1);
t9 = sin(q_t2);
t10 = sin(q_t3);
t11 = sin(q_t4);
t12 = sin(q_t5);
t13 = sin(q_t6);
t14 = q_t1+q_t2;
t17 = 1.0./sampT;
t15 = cos(t14);
t16 = sin(t14);
t18 = t2.*t3;
t19 = t2.*t9;
t20 = t3.*t8;
t21 = q_t3+q_t4+t14;
t22 = t8.*t9;
t39 = t2.*3.31509608282772e+16;
t40 = t8.*3.31509608282772e+16;
t23 = cos(t21);
t24 = q_t5+t21;
t25 = sin(t21);
t28 = -t22;
t32 = t19+t20;
t41 = qd_t1.*t39;
t42 = qd_t1.*t40;
t43 = t16.*3.208812572773657e+16;
t44 = t15.*3.208812572773657e+16;
t26 = q_t6+t24;
t27 = sin(t24);
t29 = cos(t24);
t33 = t18+t28;
t34 = t4.*t32;
t35 = t10.*t32;
t45 = qd_t1.*t44;
t46 = qd_t2.*t44;
t47 = qd_t1.*t43;
t48 = qd_t2.*t43;
t49 = t23.*3.208812572773657e+16;
t50 = t25.*3.208812572773657e+16;
t30 = cos(t26);
t31 = sin(t26);
t36 = t4.*t33;
t37 = t10.*t33;
t38 = -t35;
t51 = qd_t1.*t49;
t52 = qd_t2.*t49;
t53 = qd_t3.*t49;
t54 = qd_t4.*t49;
t55 = qd_t1.*t50;
t56 = qd_t2.*t50;
t57 = qd_t3.*t50;
t58 = qd_t4.*t50;
t61 = t29.*3.31509608282772e+16;
t62 = t27.*3.31509608282772e+16;
t59 = t30.*8.593061347490117e+15;
t60 = t31.*8.593061347490117e+15;
t75 = qd_t1.*t61;
t76 = qd_t2.*t61;
t77 = qd_t3.*t61;
t78 = qd_t4.*t61;
t79 = qd_t5.*t61;
t80 = qd_t1.*t62;
t81 = qd_t2.*t62;
t82 = qd_t3.*t62;
t83 = qd_t4.*t62;
t84 = qd_t5.*t62;
t85 = t34+t37;
t86 = t36+t38;
t89 = -t5.*(t35-t36);
t90 = -t11.*(t35-t36);
t63 = qd_t1.*t59;
t64 = qd_t2.*t59;
t65 = qd_t3.*t59;
t66 = qd_t4.*t59;
t67 = qd_t5.*t59;
t68 = qd_t6.*t59;
t69 = qd_t1.*t60;
t70 = qd_t2.*t60;
t71 = qd_t3.*t60;
t72 = qd_t4.*t60;
t73 = qd_t5.*t60;
t74 = qd_t6.*t60;
t87 = t5.*t85;
t88 = t11.*t85;
t92 = t59+t61;
t93 = t60+t62;
t91 = -t88;
t94 = t87+t90;
t98 = t50+t93;
t99 = t49+t92;
t100 = -t6.*(t88+t5.*(t35-t36));
t101 = -t12.*(t88+t5.*(t35-t36));
t111 = t69+t70+t71+t72+t73+t74;
t112 = t63+t64+t65+t66+t67+t68;
t95 = t89+t91;
t96 = t6.*t94;
t97 = t12.*t94;
t103 = t43+t98;
t104 = t44+t99;
t107 = t30.*t99.*1.700474910257266e-18;
t108 = t31.*t98.*1.700474910257266e-18;
t115 = t98.*t112.*9.894465089289424e-35;
t116 = t99.*t111.*9.894465089289424e-35;
t119 = t80+t81+t82+t83+t84+t111;
t120 = t75+t76+t77+t78+t79+t112;
t102 = -t97;
t105 = t39+t104;
t106 = t40+t103;
t109 = t96+t101;
t114 = t13.*(t97+t6.*(t88+t5.*(t35-t36))).*(-3.005100506766376e-1);
t117 = -t115;
t118 = t107+t108;
t122 = t55+t56+t57+t58+t119;
t123 = t51+t52+t53+t54+t120;
t110 = t100+t102;
t113 = t7.*t109.*3.005100506766376e-1;
t121 = t17.*t118;
t124 = t45+t46+t123;
t125 = t47+t48+t122;
t126 = t30.*t122.*1.700474910257266e-18;
t127 = t31.*t123.*1.700474910257266e-18;
t138 = t111.*t122.*4.947232544644712e-35;
t139 = t112.*t123.*4.947232544644712e-35;
t128 = -t127;
t129 = t41+t124;
t130 = t42+t125;
t131 = t31.*t129.*8.50237455128633e-19;
t132 = t31.*t129.*1.700474910257266e-18;
t133 = t30.*t130.*8.50237455128633e-19;
t134 = t30.*t130.*1.700474910257266e-18;
t140 = t111.*t130.*4.947232544644712e-35;
t141 = t112.*t129.*4.947232544644712e-35;
t135 = -t132;
t136 = -t134;
t137 = -t133;
t142 = -t140;
t143 = -t141;
t144 = t116+t117+t121+t131+t137;
t146 = t126+t128+t132+t136;
t145 = t17.*t144;
t147 = (t17.*t146)./2.0;
t148 = t113+t114+t138+t139+t142+t143+t145+t147;
dL16_dq1 = [t113+t114+t17.*(t131+t137+t105.*t111.*9.894465089289424e-35-t106.*t112.*9.894465089289424e-35+t17.*(t30.*t105.*1.700474910257266e-18+t31.*t106.*1.700474910257266e-18));t113+t114+t142+t143+t17.*(t131+t137-t103.*t112.*9.894465089289424e-35+t104.*t111.*9.894465089289424e-35+t17.*(t30.*t104.*1.700474910257266e-18+t31.*t103.*1.700474910257266e-18))+t111.*t125.*4.947232544644712e-35+t112.*t124.*4.947232544644712e-35+(t17.*(t132+t136+t30.*t125.*1.700474910257266e-18-t31.*t124.*1.700474910257266e-18))./2.0;t148;t148;t113+t114+t142+t143+t17.*(t131+t137+t92.*t111.*9.894465089289424e-35-t93.*t112.*9.894465089289424e-35+t17.*(t30.*t92.*1.700474910257266e-18+t31.*t93.*1.700474910257266e-18))+t111.*t119.*4.947232544644712e-35+t112.*t120.*4.947232544644712e-35+(t17.*(t132+t136+t30.*t119.*1.700474910257266e-18-t31.*t120.*1.700474910257266e-18))./2.0;t113+t114+t142+t143+t17.*(t131+t137+t30.*t111.*8.50237455128633e-19-t31.*t112.*8.50237455128633e-19+t17.*(t30.^2.*1.461228522370844e-2+t31.^2.*1.461228522370844e-2))+(t17.*(t132+t136+t30.*t111.*1.700474910257266e-18-t31.*t112.*1.700474910257266e-18))./2.0+t111.^2.*4.947232544644712e-35+t112.^2.*4.947232544644712e-35];
