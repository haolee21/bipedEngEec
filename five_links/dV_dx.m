function dV_dx = dV_dx(dq1,dq2,dq3,dq4,dq5,q2,q3,q4,q5)
%DV_DX
%    DV_DX = DV_DX(DQ1,DQ2,DQ3,DQ4,DQ5,Q2,Q3,Q4,Q5)

%    This function was generated by the Symbolic Math Toolbox version 8.4.
%    12-Apr-2020 23:51:52

t2 = cos(q2);
t3 = cos(q3);
t4 = cos(q5);
t5 = sin(q2);
t6 = sin(q3);
t7 = sin(q5);
t8 = q2+q3;
t9 = q3+q4;
t10 = dq1.^2;
t11 = dq2.^2;
t12 = dq3.^2;
t13 = dq4.^2;
t14 = dq5.^2;
t26 = dq1+dq2+dq3+dq4;
t15 = cos(t8);
t16 = cos(t9);
t17 = q4+t8;
t18 = q5+t9;
t19 = sin(t8);
t20 = sin(t9);
t27 = dq5+t26;
t33 = dq1.*t5.*2.9838375e+1;
t34 = dq2.*t5.*2.9838375e+1;
t35 = dq1.*t6.*2.362122e+1;
t36 = dq2.*t6.*2.362122e+1;
t37 = dq3.*t6.*2.362122e+1;
t38 = dq1.*dq3.*t3.*2.362122e+1;
t39 = dq2.*dq3.*t3.*2.362122e+1;
t44 = dq1.*t7.*8.76015e-1;
t45 = dq2.*t7.*8.76015e-1;
t46 = dq3.*t7.*8.76015e-1;
t47 = dq4.*t7.*8.76015e-1;
t48 = dq5.*t7.*8.76015e-1;
t49 = t3.*t12.*1.181061e+1;
t52 = dq1.*dq5.*t4.*8.76015e-1;
t53 = dq2.*dq5.*t4.*8.76015e-1;
t54 = dq3.*dq5.*t4.*8.76015e-1;
t55 = dq4.*dq5.*t4.*8.76015e-1;
t84 = t4.*t14.*4.380075e-1;
t111 = t7.*t26.*8.76015e-1;
t21 = cos(t17);
t22 = cos(t18);
t23 = q5+t17;
t24 = sin(t17);
t25 = sin(t18);
t32 = t15.*1.063805175900365e+17;
t40 = -t34;
t41 = -t35;
t42 = -t36;
t43 = -t37;
t50 = -t38;
t51 = -t39;
t56 = dq1.*t19.*2.362122e+1;
t57 = dq2.*t19.*2.362122e+1;
t58 = dq3.*t19.*2.362122e+1;
t59 = -t48;
t60 = -t49;
t61 = dq1.*t20.*3.240405;
t62 = dq2.*t20.*3.240405;
t63 = dq3.*t20.*3.240405;
t64 = dq4.*t20.*3.240405;
t65 = -t52;
t66 = -t53;
t67 = -t54;
t68 = -t55;
t72 = dq1.*dq2.*t15.*2.362122e+1;
t73 = dq1.*dq3.*t15.*2.362122e+1;
t74 = dq2.*dq3.*t15.*2.362122e+1;
t85 = dq1.*dq2.*t16.*3.240405;
t86 = dq1.*dq3.*t16.*3.240405;
t87 = dq1.*dq4.*t16.*3.240405;
t88 = dq2.*dq3.*t16.*3.240405;
t89 = dq2.*dq4.*t16.*3.240405;
t90 = dq3.*dq4.*t16.*3.240405;
t95 = t10.*t15.*1.181061e+1;
t96 = t11.*t15.*1.181061e+1;
t97 = t12.*t15.*1.181061e+1;
t106 = -t84;
t107 = t10.*t16.*1.6202025;
t108 = t11.*t16.*1.6202025;
t109 = t12.*t16.*1.6202025;
t110 = t13.*t16.*1.6202025;
t160 = t7.*t27.*8.76015e-1;
t28 = cos(t23);
t29 = sin(t23);
t31 = t21.*1.459348675052949e+16;
t69 = -t56;
t70 = -t57;
t71 = -t58;
t75 = -t61;
t76 = -t62;
t77 = -t63;
t78 = -t64;
t79 = dq1.*t25.*8.76015e-1;
t80 = dq2.*t25.*8.76015e-1;
t81 = dq3.*t25.*8.76015e-1;
t82 = dq4.*t25.*8.76015e-1;
t83 = dq5.*t25.*8.76015e-1;
t91 = dq1.*t24.*3.240405;
t92 = dq2.*t24.*3.240405;
t93 = dq3.*t24.*3.240405;
t94 = dq4.*t24.*3.240405;
t98 = -t72;
t99 = -t73;
t100 = -t74;
t112 = -t86;
t113 = -t87;
t114 = -t88;
t115 = -t89;
t116 = -t90;
t117 = dq1.*dq2.*t22.*8.76015e-1;
t118 = dq1.*dq3.*t22.*8.76015e-1;
t119 = dq1.*dq4.*t22.*8.76015e-1;
t120 = dq2.*dq3.*t22.*8.76015e-1;
t121 = dq1.*dq5.*t22.*8.76015e-1;
t122 = dq2.*dq4.*t22.*8.76015e-1;
t123 = dq2.*dq5.*t22.*8.76015e-1;
t124 = dq3.*dq4.*t22.*8.76015e-1;
t125 = dq3.*dq5.*t22.*8.76015e-1;
t126 = dq4.*dq5.*t22.*8.76015e-1;
t136 = -t96;
t137 = -t97;
t138 = dq1.*dq2.*t21.*3.240405;
t139 = dq1.*dq3.*t21.*3.240405;
t140 = dq1.*dq4.*t21.*3.240405;
t141 = dq2.*dq3.*t21.*3.240405;
t142 = dq2.*dq4.*t21.*3.240405;
t143 = dq3.*dq4.*t21.*3.240405;
t144 = -t109;
t145 = -t110;
t161 = t10.*t21.*1.6202025;
t162 = t11.*t21.*1.6202025;
t163 = t12.*t21.*1.6202025;
t164 = t13.*t21.*1.6202025;
t181 = t10.*t22.*4.380075e-1;
t182 = t11.*t22.*4.380075e-1;
t183 = t12.*t22.*4.380075e-1;
t184 = t13.*t22.*4.380075e-1;
t185 = t14.*t22.*4.380075e-1;
t186 = -t160;
t30 = t28.*3.945220827570965e+15;
t101 = -t79;
t102 = -t80;
t103 = -t81;
t104 = -t82;
t105 = -t83;
t127 = -t91;
t128 = -t92;
t129 = -t93;
t130 = -t94;
t131 = dq1.*t29.*8.76015e-1;
t132 = dq2.*t29.*8.76015e-1;
t133 = dq3.*t29.*8.76015e-1;
t134 = dq4.*t29.*8.76015e-1;
t135 = dq5.*t29.*8.76015e-1;
t146 = -t118;
t147 = -t119;
t148 = -t120;
t149 = -t121;
t150 = -t122;
t151 = -t123;
t152 = -t124;
t153 = -t125;
t154 = -t126;
t165 = -t138;
t166 = -t139;
t167 = -t140;
t168 = -t141;
t169 = -t142;
t170 = -t143;
t171 = dq1.*dq2.*t28.*8.76015e-1;
t172 = dq1.*dq3.*t28.*8.76015e-1;
t173 = dq1.*dq4.*t28.*8.76015e-1;
t174 = dq2.*dq3.*t28.*8.76015e-1;
t175 = dq1.*dq5.*t28.*8.76015e-1;
t176 = dq2.*dq4.*t28.*8.76015e-1;
t177 = dq2.*dq5.*t28.*8.76015e-1;
t178 = dq3.*dq4.*t28.*8.76015e-1;
t179 = dq3.*dq5.*t28.*8.76015e-1;
t180 = dq4.*dq5.*t28.*8.76015e-1;
t187 = -t162;
t188 = -t163;
t189 = -t164;
t200 = -t183;
t201 = -t184;
t202 = -t185;
t203 = t10.*t28.*4.380075e-1;
t204 = t11.*t28.*4.380075e-1;
t205 = t12.*t28.*4.380075e-1;
t206 = t13.*t28.*4.380075e-1;
t207 = t14.*t28.*4.380075e-1;
t155 = -t131;
t156 = -t132;
t157 = -t133;
t158 = -t134;
t159 = -t135;
t190 = -t171;
t191 = -t172;
t192 = -t173;
t193 = -t174;
t194 = -t175;
t195 = -t176;
t196 = -t177;
t197 = -t178;
t198 = -t179;
t199 = -t180;
t208 = -t204;
t209 = -t205;
t210 = -t206;
t211 = -t207;
t212 = t117+t181+t182+t203;
t213 = t85+t107+t108+t161+t212;
t214 = t65+t66+t67+t68+t106+t212;
dV_dx = reshape([0.0,t98+t99+t100+t136+t137+t165+t166+t167+t168+t169+t170+t187+t188+t189+t190+t191+t192+t193+t194+t195+t196+t197+t198+t199+t208+t209+t210+t211-t2.*t11.*1.49191875e+1-dq1.*dq2.*t2.*2.9838375e+1,t50+t51+t60+t98+t99+t100+t112+t113+t114+t115+t116+t136+t137+t144+t145+t146+t147+t148+t149+t150+t151+t152+t153+t154+t165+t166+t167+t168+t169+t170+t187+t188+t189+t190+t191+t192+t193+t194+t195+t196+t197+t198+t199+t200+t201+t202+t208+t209+t210+t211,t112+t113+t114+t115+t116+t144+t145+t146+t147+t148+t149+t150+t151+t152+t153+t154+t165+t166+t167+t168+t169+t170+t187+t188+t189+t190+t191+t192+t193+t194+t195+t196+t197+t198+t199+t200+t201+t202+t208+t209+t210+t211,t65+t66+t67+t68+t106+t146+t147+t148+t149+t150+t151+t152+t153+t154+t190+t191+t192+t193+t194+t195+t196+t197+t198+t199+t200+t201+t202+t208+t209+t210+t211,t40+t43+t59+t70+t71+t77+t78+t103+t104+t105+t128+t129+t130+t156+t157+t158+t159,-t33+t40+t43+t59+t69+t70+t71+t77+t78+t103+t104+t105+t127+t128+t129+t130+t155+t156+t157+t158+t159,t41+t42+t43+t59+t69+t70+t71+t75+t76+t77+t78+t101+t102+t103+t104+t105+t127+t128+t129+t130+t155+t156+t157+t158+t159,t59+t75+t76+t77+t78+t101+t102+t103+t104+t105+t127+t128+t129+t130+t155+t156+t157+t158+t159,t27.*(t7+t25+t29).*(-8.76015e-1),0.0,0.0,0.0,0.0,0.0,0.0,(t10.*(t2.*1.343800945313411e+17+t30+t31+t32))./9.007199254740992e+15,t50+t51+t60+t95+t112+t113+t114+t115+t116+t144+t145+t146+t147+t148+t149+t150+t151+t152+t153+t154+t161+t200+t201+t202+t203,t112+t113+t114+t115+t116+t144+t145+t146+t147+t148+t149+t150+t151+t152+t153+t154+t161+t200+t201+t202+t203,t65+t66+t67+t68+t106+t146+t147+t148+t149+t150+t151+t152+t153+t154+t200+t201+t202+t203,t33+t43+t56+t59+t77+t78+t91+t103+t104+t105+t131,t43+t59+t77+t78+t103+t104+t105,t41+t42+t43+t59+t75+t76+t77+t78+t101+t102+t103+t104+t105,t59+t75+t76+t77+t78+t101+t102+t103+t104+t105,t27.*(t7+t25).*(-8.76015e-1),0.0,0.0,0.0,0.0,0.0,0.0,(t10.*(t30+t31+t32))./9.007199254740992e+15,t95+t213+t3.*t10.*1.181061e+1+t3.*t11.*1.181061e+1+dq1.*dq2.*t3.*2.362122e+1,t213,t214,t35+t36+t56+t59+t61+t62+t79+t80+t91+t131,t35+t36+t59+t61+t62+t79+t80,t59,t59,t186,0.0,0.0,0.0,0.0,0.0,0.0,(t10.*(t30+t31))./9.007199254740992e+15,t213,t213,t214,t59+t61+t62+t79+t80+t91+t131,t59+t61+t62+t79+t80,t59,t59,t186,0.0,0.0,0.0,0.0,0.0,0.0,t203,t212,t212,t212+t4.*t10.*4.380075e-1+t4.*t11.*4.380075e-1+t4.*t12.*4.380075e-1+t4.*t13.*4.380075e-1+dq1.*dq2.*t4.*8.76015e-1+dq1.*dq3.*t4.*8.76015e-1+dq1.*dq4.*t4.*8.76015e-1+dq2.*dq3.*t4.*8.76015e-1+dq2.*dq4.*t4.*8.76015e-1+dq3.*dq4.*t4.*8.76015e-1,t44+t45+t46+t47+t79+t80+t131,t44+t45+t46+t47+t79+t80,t111,t111,0.0,0.0,0.0,0.0,0.0,0.0],[15,5]);
