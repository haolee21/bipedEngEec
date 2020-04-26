function dbeta_dx5 = dbeta_dx5(q1,q2,q3,q4,q5,th)
%DBETA_DX5
%    DBETA_DX5 = DBETA_DX5(Q1,Q2,Q3,Q4,Q5,TH)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    25-Apr-2020 22:24:42

t2 = cos(q5);
t3 = sin(q1);
t4 = sin(q5);
t5 = q1+q2;
t8 = q3+q4+q5;
t14 = th.*4.0e+2;
t6 = t2.^2;
t7 = t4.^2;
t9 = sin(t5);
t10 = cos(t8);
t11 = q3+q4+t5;
t12 = q2+t8;
t13 = sin(t8);
t15 = t2.*(2.0./5.0);
t16 = t4.*(2.0./5.0);
t17 = t3.*1.72e+2;
t19 = t5+t8;
t22 = -t14;
t18 = cos(t12);
t20 = sin(t11);
t21 = sin(t12);
t23 = sin(t19);
t24 = t10.^2;
t25 = t13.^2;
t26 = cos(t19);
t27 = t9.*1.6e+2;
t32 = t10.*(2.0./5.0);
t33 = t13.*(2.0./5.0);
t28 = t18.^2;
t29 = t21.^2;
t30 = t26.^2;
t31 = t23.^2;
t34 = t20.*1.6e+2;
t35 = t23.*1.72e+2;
t36 = t15.*t26;
t37 = t15.*t23;
t38 = t16.*t26;
t39 = t16.*t23;
t40 = t18.*(4.3e+1./1.0e+2);
t41 = t21.*(4.3e+1./1.0e+2);
t42 = t4.*t26.*(-2.0./5.0);
t43 = t15+t32;
t44 = t16+t33;
t45 = t26.*t43;
t46 = t23.*t43;
t47 = t26.*t44;
t48 = t23.*t44;
t49 = t7.*t24.*t30.*6.4e+3;
t50 = t6.*t25.*t30.*6.4e+3;
t51 = t7.*t24.*t31.*6.4e+3;
t52 = t6.*t25.*t31.*6.4e+3;
t54 = t7.*t28.*t30.*5.547e+3;
t55 = t6.*t29.*t30.*5.547e+3;
t56 = t7.*t28.*t31.*5.547e+3;
t57 = t6.*t29.*t31.*5.547e+3;
t58 = t36+t39;
t59 = t37+t42;
t60 = t2.*t4.*t10.*t21.*t30.*6.88e+3;
t61 = t2.*t4.*t13.*t18.*t30.*6.88e+3;
t62 = t41+t44;
t63 = t24.*t29.*t30.*1.849e+3;
t64 = t25.*t28.*t30.*1.849e+3;
t65 = t2.*t4.*t10.*t21.*t31.*6.88e+3;
t66 = t2.*t4.*t13.*t18.*t31.*6.88e+3;
t67 = t24.*t29.*t31.*1.849e+3;
t68 = t25.*t28.*t31.*1.849e+3;
t69 = t2.*t10.*t29.*t30.*3.698e+3;
t70 = t7.*t10.*t18.*t30.*6.88e+3;
t71 = t4.*t13.*t28.*t30.*3.698e+3;
t72 = t2.*t10.*t29.*t31.*3.698e+3;
t73 = t6.*t13.*t21.*t30.*6.88e+3;
t74 = t7.*t10.*t18.*t31.*6.88e+3;
t75 = t2.*t4.*t10.*t13.*t30.*1.28e+4;
t76 = t4.*t13.*t28.*t31.*3.698e+3;
t77 = t6.*t13.*t21.*t31.*6.88e+3;
t78 = t2.*t4.*t10.*t13.*t31.*1.28e+4;
t79 = t40+t43;
t86 = t2.*t4.*t18.*t21.*t30.*1.1094e+4;
t87 = t2.*t4.*t18.*t21.*t31.*1.1094e+4;
t90 = t4.*t10.*t18.*t21.*t30.*3.698e+3;
t91 = t2.*t13.*t18.*t21.*t30.*3.698e+3;
t92 = t4.*t10.*t18.*t21.*t31.*3.698e+3;
t93 = t2.*t13.*t18.*t21.*t31.*3.698e+3;
t98 = t10.*t13.*t18.*t21.*t30.*3.698e+3;
t99 = t10.*t13.*t18.*t21.*t31.*3.698e+3;
t107 = t17+t22+t27+t34+t35;
t53 = -t47;
t80 = -t60;
t81 = -t61;
t82 = -t65;
t83 = -t66;
t84 = -t75;
t85 = -t78;
t88 = -t86;
t89 = -t87;
t94 = -t90;
t95 = -t91;
t96 = -t92;
t97 = -t93;
t100 = -t98;
t101 = -t99;
t102 = t26.*t79;
t103 = t23.*t79;
t104 = t26.*t62;
t105 = t23.*t62;
t108 = tanh(t107);
t111 = t45+t48;
t106 = -t104;
t109 = t108.^2;
t112 = t46+t53;
t139 = t102+t105;
t277 = t49+t50+t51+t52+t54+t55+t56+t57+t63+t64+t67+t68+t69+t70+t71+t72+t73+t74+t76+t77+t80+t81+t82+t83+t84+t85+t88+t89+t94+t95+t96+t97+t100+t101;
t110 = t109-1.0;
t159 = t103+t106;
t278 = 1.0./t277;
t113 = t7.*t10.*t30.*t110.*1.376e+4;
t114 = t6.*t13.*t30.*t110.*1.376e+4;
t115 = t2.*t4.*t10.*t30.*t110.*1.376e+4;
t116 = t2.*t4.*t13.*t30.*t110.*1.376e+4;
t117 = t2.*t4.*t18.*t30.*t110.*2.2188e+4;
t118 = t2.*t4.*t21.*t30.*t110.*2.2188e+4;
t120 = t7.*t18.*t30.*t110.*2.2188e+4;
t121 = t6.*t21.*t30.*t110.*2.2188e+4;
t124 = t2.*t4.*t18.*t30.*t110.*2.9584e+5;
t125 = t2.*t4.*t21.*t30.*t110.*2.9584e+5;
t126 = t7.*t10.*t30.*t110.*5.504e+5;
t127 = t6.*t13.*t30.*t110.*5.504e+5;
t128 = t4.*t24.*t30.*t110.*5.504e+5;
t129 = t2.*t25.*t30.*t110.*5.504e+5;
t131 = t7.*t18.*t30.*t110.*2.9584e+5;
t132 = t6.*t21.*t30.*t110.*2.9584e+5;
t133 = t4.*t28.*t30.*t110.*3.18028e+5;
t134 = t2.*t29.*t30.*t110.*3.18028e+5;
t135 = t2.*t4.*t10.*t30.*t110.*5.504e+5;
t136 = t2.*t4.*t13.*t30.*t110.*5.504e+5;
t142 = t4.*t10.*t18.*t30.*t110.*7.396e+3;
t143 = t2.*t13.*t18.*t30.*t110.*7.396e+3;
t144 = t4.*t10.*t21.*t30.*t110.*7.396e+3;
t145 = t2.*t13.*t21.*t30.*t110.*7.396e+3;
t152 = t21.*t24.*t30.*t110.*7.396e+3;
t153 = t18.*t25.*t30.*t110.*7.396e+3;
t156 = t2.*t10.*t21.*t30.*t110.*1.4792e+4;
t157 = t4.*t13.*t18.*t30.*t110.*1.4792e+4;
t160 = t2.*t4.*t10.*t23.*t26.*t110.*1.376e+4;
t161 = t4.*t10.*t18.*t30.*t110.*2.9584e+5;
t162 = t2.*t10.*t21.*t30.*t110.*2.9584e+5;
t163 = t2.*t13.*t18.*t30.*t110.*2.9584e+5;
t164 = t2.*t4.*t13.*t23.*t26.*t110.*1.376e+4;
t166 = t4.*t10.*t21.*t30.*t110.*2.9584e+5;
t167 = t4.*t13.*t18.*t30.*t110.*2.9584e+5;
t168 = t2.*t13.*t21.*t30.*t110.*2.9584e+5;
t169 = t7.*t10.*t23.*t26.*t110.*1.376e+4;
t170 = t6.*t13.*t23.*t26.*t110.*1.376e+4;
t171 = t21.*t24.*t30.*t110.*2.9584e+5;
t172 = t18.*t25.*t30.*t110.*2.9584e+5;
t173 = t10.*t29.*t30.*t110.*3.18028e+5;
t174 = t13.*t28.*t30.*t110.*3.18028e+5;
t175 = t2.*t10.*t13.*t30.*t110.*5.504e+5;
t176 = t4.*t10.*t13.*t30.*t110.*5.504e+5;
t177 = t7.*t18.*t23.*t26.*t110.*2.2188e+4;
t178 = t6.*t21.*t23.*t26.*t110.*2.2188e+4;
t179 = t4.*t10.*t18.*t30.*t110.*5.9168e+5;
t180 = t2.*t10.*t21.*t30.*t110.*5.9168e+5;
t181 = t4.*t13.*t18.*t30.*t110.*5.9168e+5;
t182 = t2.*t13.*t21.*t30.*t110.*5.9168e+5;
t183 = t10.*t13.*t18.*t30.*t110.*7.396e+3;
t184 = t10.*t13.*t21.*t30.*t110.*7.396e+3;
t185 = t2.*t4.*t18.*t23.*t26.*t110.*2.2188e+4;
t186 = t2.*t18.*t21.*t30.*t110.*3.18028e+5;
t187 = t2.*t4.*t21.*t23.*t26.*t110.*2.2188e+4;
t190 = t4.*t18.*t21.*t30.*t110.*3.18028e+5;
t195 = t7.*t18.*t23.*t26.*t110.*2.9584e+5;
t196 = t6.*t21.*t23.*t26.*t110.*2.9584e+5;
t197 = t4.*t23.*t26.*t28.*t110.*3.18028e+5;
t198 = t2.*t23.*t26.*t29.*t110.*3.18028e+5;
t199 = t2.*t4.*t10.*t23.*t26.*t110.*5.504e+5;
t202 = t2.*t4.*t13.*t23.*t26.*t110.*5.504e+5;
t206 = t2.*t4.*t18.*t23.*t26.*t110.*2.9584e+5;
t208 = t2.*t4.*t21.*t23.*t26.*t110.*2.9584e+5;
t210 = t7.*t10.*t23.*t26.*t110.*5.504e+5;
t211 = t6.*t13.*t23.*t26.*t110.*5.504e+5;
t212 = t4.*t23.*t24.*t26.*t110.*5.504e+5;
t213 = t2.*t23.*t25.*t26.*t110.*5.504e+5;
t214 = t4.*t10.*t21.*t23.*t26.*t110.*7.396e+3;
t215 = t2.*t13.*t21.*t23.*t26.*t110.*7.396e+3;
t221 = t21.*t23.*t24.*t26.*t110.*7.396e+3;
t222 = t18.*t23.*t25.*t26.*t110.*7.396e+3;
t223 = t10.*t13.*t18.*t30.*t110.*2.9584e+5;
t224 = t10.*t13.*t21.*t30.*t110.*2.9584e+5;
t229 = t4.*t10.*t18.*t23.*t26.*t110.*7.396e+3;
t230 = t2.*t13.*t18.*t23.*t26.*t110.*7.396e+3;
t231 = t10.*t18.*t21.*t30.*t110.*3.18028e+5;
t232 = t4.*t13.*t18.*t23.*t26.*t110.*1.4792e+4;
t233 = t13.*t18.*t21.*t30.*t110.*3.18028e+5;
t236 = t2.*t10.*t21.*t23.*t26.*t110.*1.4792e+4;
t239 = t4.*t10.*t21.*t23.*t26.*t110.*2.9584e+5;
t240 = t4.*t13.*t18.*t23.*t26.*t110.*2.9584e+5;
t241 = t2.*t13.*t21.*t23.*t26.*t110.*2.9584e+5;
t243 = t21.*t23.*t24.*t26.*t110.*2.9584e+5;
t244 = t18.*t23.*t25.*t26.*t110.*2.9584e+5;
t245 = t10.*t23.*t26.*t29.*t110.*3.18028e+5;
t246 = t13.*t23.*t26.*t28.*t110.*3.18028e+5;
t247 = t2.*t10.*t13.*t23.*t26.*t110.*5.504e+5;
t248 = t4.*t10.*t13.*t23.*t26.*t110.*5.504e+5;
t250 = t4.*t10.*t18.*t23.*t26.*t110.*2.9584e+5;
t251 = t2.*t10.*t21.*t23.*t26.*t110.*2.9584e+5;
t252 = t2.*t13.*t18.*t23.*t26.*t110.*2.9584e+5;
t256 = t4.*t13.*t18.*t23.*t26.*t110.*5.9168e+5;
t257 = t2.*t13.*t21.*t23.*t26.*t110.*5.9168e+5;
t258 = t10.*t13.*t18.*t23.*t26.*t110.*7.396e+3;
t259 = t10.*t13.*t21.*t23.*t26.*t110.*7.396e+3;
t260 = t2.*t18.*t21.*t23.*t26.*t110.*3.18028e+5;
t261 = t4.*t18.*t21.*t23.*t26.*t110.*3.18028e+5;
t265 = t4.*t10.*t18.*t23.*t26.*t110.*5.9168e+5;
t266 = t2.*t10.*t21.*t23.*t26.*t110.*5.9168e+5;
t271 = t10.*t13.*t18.*t23.*t26.*t110.*2.9584e+5;
t272 = t10.*t13.*t21.*t23.*t26.*t110.*2.9584e+5;
t274 = t10.*t18.*t21.*t23.*t26.*t110.*3.18028e+5;
t275 = t13.*t18.*t21.*t23.*t26.*t110.*3.18028e+5;
t119 = -t113;
t122 = -t115;
t123 = -t117;
t130 = -t120;
t137 = -t124;
t138 = -t125;
t140 = -t128;
t141 = -t129;
t146 = -t131;
t147 = -t132;
t148 = -t133;
t149 = -t134;
t150 = -t135;
t151 = -t136;
t154 = -t142;
t155 = -t143;
t158 = -t153;
t165 = -t157;
t188 = -t169;
t189 = -t170;
t191 = -t171;
t192 = -t172;
t193 = -t177;
t194 = -t178;
t200 = -t179;
t201 = -t180;
t203 = -t181;
t204 = -t182;
t205 = -t183;
t207 = -t186;
t209 = -t190;
t216 = -t195;
t217 = -t196;
t218 = -t197;
t219 = -t198;
t220 = -t199;
t225 = -t206;
t226 = -t208;
t227 = -t210;
t228 = -t212;
t234 = -t221;
t235 = -t222;
t237 = -t231;
t238 = -t232;
t242 = -t233;
t249 = -t236;
t253 = -t239;
t254 = -t240;
t255 = -t241;
t262 = -t243;
t263 = -t245;
t264 = -t248;
t267 = -t260;
t268 = -t261;
t269 = -t265;
t270 = -t266;
t273 = -t272;
t276 = -t274;
t279 = t116+t118+t119+t130+t144+t145+t158+t160+t165+t184+t185+t189+t194+t229+t230+t234+t249+t258;
t280 = t114+t121+t122+t123+t152+t154+t155+t156+t164+t187+t188+t193+t205+t214+t215+t235+t238+t259;
t281 = t126+t131+t138+t141+t149+t151+t166+t167+t176+t190+t196+t204+t211+t218+t220+t225+t228+t247+t251+t252+t260+t269;
t282 = t127+t132+t137+t140+t148+t150+t162+t163+t175+t186+t198+t200+t202+t208+t213+t216+t227+t253+t254+t257+t264+t268;
t283 = t125+t126+t134+t146+t151+t166+t168+t173+t192+t197+t203+t206+t209+t211+t217+t220+t224+t242+t246+t250+t252+t262+t267+t270+t271+t276;
t284 = t124+t127+t133+t147+t150+t161+t163+t174+t191+t195+t201+t202+t207+t219+t223+t226+t227+t237+t244+t253+t255+t256+t261+t263+t273+t275;
t290 = t278.*(t38-t2.*t23.*(2.0./5.0)).*(t124-t127+t128+t133+t135+t147-t162-t163-t175+t179+t195-t202+t207+t210-t213+t219+t226+t239+t240+t248-t257+t261);
t294 = -t112.*t278.*(t124-t127+t128+t133+t135+t147-t162-t163-t175+t179+t195-t202+t207+t210-t213+t219+t226+t239+t240+t248-t257+t261);
t297 = -t159.*t278.*(t124-t127+t128+t133+t135+t147-t162-t163-t175+t179+t195-t202+t207+t210-t213+t219+t226+t239+t240+t248-t257+t261);
t285 = t58.*t278.*t279.*4.0e+1;
t287 = t278.*t280.*(t38-t2.*t23.*(2.0./5.0)).*-4.0e+1;
t288 = t278.*t280.*(t38-t2.*t23.*(2.0./5.0)).*4.0e+1;
t289 = t58.*t278.*t281;
t292 = t111.*t278.*t281;
t295 = t139.*t278.*t281;
t298 = t58.*t278.*t283;
t299 = -t278.*t284.*(t38-t2.*t23.*(2.0./5.0));
t300 = t278.*t284.*(t38-t2.*t23.*(2.0./5.0));
t286 = -t285;
t291 = -t289;
t293 = -t292;
t296 = -t295;
t305 = t298+t300;
t301 = t286+t288;
t302 = t290+t291;
t303 = t293+t294;
t304 = t296+t297;
dbeta_dx5 = reshape([t139.*t278.*t279.*-4.0e+1-t159.*t278.*t280.*4.0e+1,t139.*t278.*t283-t159.*t278.*t284,t304,t304,0.0,t111.*t278.*t279.*-4.0e+1-t112.*t278.*t280.*4.0e+1,t111.*t278.*t283-t112.*t278.*t284,t303,t303,0.0,t301,t305,t302,t302,0.0,t301,t305,t302,t302,0.0,0.0,0.0,0.0,0.0,0.0],[5,5]);
