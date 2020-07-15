function dL23_dq1 = dL23_dq1(in1,in2,sampT)
%DL23_DQ1
%    DL23_DQ1 = DL23_DQ1(IN1,IN2,SAMPT)

%    This function was generated by the Symbolic Math Toolbox version 8.4.
%    14-Jul-2020 11:58:56

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
t18 = 1.0./sampT;
t15 = cos(t14);
t16 = q_t3+t14;
t17 = sin(t14);
t19 = t2.*t3;
t21 = t2.*t9;
t22 = t3.*t8;
t25 = t8.*t9;
t44 = t2.*2.004946259110003e+18;
t45 = t8.*2.004946259110003e+18;
t46 = t2.*4.811871021864006e+19;
t47 = t8.*4.811871021864006e+19;
t51 = t2.*9.623742043728013e+19;
t53 = t8.*9.623742043728013e+19;
t20 = cos(t16);
t23 = q_t4+t16;
t24 = sin(t16);
t31 = -t25;
t37 = t21+t22;
t48 = qd_t1.*t44;
t49 = qd_t1.*t45;
t50 = t15.*1.940666756928512e+18;
t52 = t17.*1.940666756928512e+18;
t54 = qd_t1.*t46;
t55 = qd_t1.*t47;
t56 = t15.*4.65760021662843e+19;
t57 = t17.*4.65760021662843e+19;
t61 = qd_t1.*t51;
t65 = qd_t1.*t53;
t72 = t15.*9.31520043325686e+19;
t73 = t17.*9.31520043325686e+19;
t26 = cos(t23);
t27 = q_t5+t23;
t28 = sin(t23);
t38 = t19+t31;
t39 = t4.*t37;
t40 = t10.*t37;
t58 = qd_t1.*t50;
t59 = qd_t2.*t50;
t60 = t20.*2.026117081401678e+18;
t62 = qd_t1.*t52;
t63 = qd_t2.*t52;
t64 = t24.*2.026117081401678e+18;
t74 = qd_t1.*t56;
t75 = qd_t2.*t56;
t76 = qd_t1.*t57;
t77 = qd_t2.*t57;
t82 = qd_t1.*t72;
t83 = qd_t2.*t72;
t84 = qd_t1.*t73;
t85 = qd_t2.*t73;
t29 = q_t6+t27;
t30 = sin(t27);
t32 = cos(t27);
t35 = t26.^2;
t36 = t28.^2;
t41 = t4.*t38;
t42 = t10.*t38;
t43 = -t40;
t66 = qd_t1.*t60;
t67 = qd_t2.*t60;
t68 = qd_t3.*t60;
t69 = qd_t1.*t64;
t70 = qd_t2.*t64;
t71 = qd_t3.*t64;
t78 = t26.*4.65760021662843e+19;
t79 = t28.*4.65760021662843e+19;
t80 = t26.*2.01674089380011e+19;
t81 = t28.*2.01674089380011e+19;
t96 = t26.*9.31520043325686e+19;
t97 = t28.*9.31520043325686e+19;
t157 = t50+t60;
t158 = t52+t64;
t178 = t39.*3.30984181582429e+1;
t33 = cos(t29);
t34 = sin(t29);
t86 = t32.*2.083540152467115e+19;
t87 = t30.*2.083540152467115e+19;
t88 = qd_t1.*t80;
t89 = qd_t2.*t80;
t90 = qd_t3.*t80;
t91 = qd_t4.*t80;
t92 = qd_t1.*t81;
t93 = qd_t2.*t81;
t94 = qd_t3.*t81;
t95 = qd_t4.*t81;
t98 = qd_t1.*t78;
t99 = qd_t2.*t78;
t100 = qd_t3.*t78;
t101 = qd_t4.*t78;
t102 = qd_t1.*t79;
t103 = qd_t2.*t79;
t104 = qd_t3.*t79;
t105 = qd_t4.*t79;
t106 = t32.*9.623742043728013e+19;
t107 = t30.*9.623742043728013e+19;
t118 = qd_t1.*t96;
t119 = qd_t2.*t96;
t120 = qd_t3.*t96;
t121 = qd_t4.*t96;
t123 = qd_t1.*t97;
t124 = qd_t2.*t97;
t125 = qd_t3.*t97;
t126 = qd_t4.*t97;
t150 = t39+t42;
t151 = t41+t43;
t154 = -t5.*(t40-t41);
t155 = -t11.*(t40-t41);
t159 = t56+t80;
t160 = t57+t81;
t168 = t44+t157;
t169 = t45+t158;
t179 = t66+t67+t68;
t180 = t42.*3.30984181582429e+1;
t181 = t69+t70+t71;
t184 = t35.*2.816797234966514e-1;
t185 = t36.*2.816797234966514e-1;
t191 = t11.*(t40-t41).*(-7.442089063853112);
t108 = qd_t1.*t86;
t109 = qd_t2.*t86;
t110 = qd_t3.*t86;
t111 = qd_t4.*t86;
t112 = qd_t5.*t86;
t113 = qd_t1.*t87;
t114 = qd_t2.*t87;
t115 = qd_t3.*t87;
t116 = qd_t4.*t87;
t117 = qd_t5.*t87;
t122 = t33.*2.49457040483832e+19;
t127 = t34.*2.49457040483832e+19;
t140 = qd_t1.*t106;
t141 = qd_t2.*t106;
t142 = qd_t3.*t106;
t143 = qd_t4.*t106;
t144 = qd_t5.*t106;
t145 = qd_t1.*t107;
t146 = qd_t2.*t107;
t147 = qd_t3.*t107;
t148 = qd_t4.*t107;
t149 = qd_t5.*t107;
t152 = t5.*t150;
t153 = t11.*t150;
t161 = t78+t86;
t162 = t79+t87;
t176 = t46+t159;
t177 = t47+t160;
t197 = t88+t89+t90+t91;
t199 = t92+t93+t94+t95;
t202 = t62+t63+t181;
t203 = t58+t59+t179;
t128 = qd_t1.*t122;
t129 = qd_t2.*t122;
t130 = qd_t3.*t122;
t131 = qd_t4.*t122;
t132 = qd_t5.*t122;
t133 = qd_t6.*t122;
t134 = qd_t1.*t127;
t135 = qd_t2.*t127;
t136 = qd_t3.*t127;
t137 = qd_t4.*t127;
t138 = qd_t5.*t127;
t139 = qd_t6.*t127;
t156 = -t153;
t163 = t162.^2;
t164 = t161.^2;
t165 = t106+t122;
t166 = t107+t127;
t167 = t152+t155;
t173 = -t6.*(t153+t5.*(t40-t41));
t174 = -t12.*(t153+t5.*(t40-t41));
t182 = t57+t162;
t183 = t56+t161;
t190 = t152.*7.442089063853112;
t200 = t197.^2;
t201 = t199.^2;
t209 = t108+t109+t110+t111+t112;
t210 = t113+t114+t115+t116+t117;
t211 = t12.*(t153+t5.*(t40-t41)).*(-2.10870039466251);
t213 = t48+t203;
t214 = t49+t202;
t217 = t74+t75+t197;
t218 = t76+t77+t199;
t221 = t28.*t197.*6.983537755459679e-21;
t222 = t28.*t197.*1.396707551091936e-20;
t223 = t26.*t199.*6.983537755459679e-21;
t224 = t26.*t199.*1.396707551091936e-20;
t170 = t154+t156;
t171 = t6.*t167;
t172 = t12.*t167;
t186 = t96+t165;
t187 = t97+t166;
t192 = t46+t183;
t193 = t47+t182;
t194 = t164.*2.117775986156328e-40;
t195 = t163.*2.117775986156328e-40;
t216 = t201.*1.731391914779077e-40;
t219 = t200.*1.731391914779077e-40;
t220 = t134+t135+t136+t137+t138+t139;
t225 = t128+t129+t130+t131+t132+t133;
t226 = -t224;
t227 = -t223;
t228 = t54+t217;
t229 = t55+t218;
t230 = t20.*t214.*2.876800071864181e-18;
t231 = t24.*t213.*2.876800071864181e-18;
t232 = t20.*t214.*1.438400035932091e-18;
t234 = t24.*t213.*1.438400035932091e-18;
t244 = t98+t99+t100+t101+t209;
t245 = t102+t103+t104+t105+t210;
t248 = t181.*t214.*3.549646881553849e-37;
t249 = t179.*t213.*3.549646881553849e-37;
t175 = -t172;
t188 = t186.^2;
t189 = t187.^2;
t196 = t72+t186;
t198 = t73+t187;
t208 = t171.*2.10870039466251;
t212 = t171+t174;
t233 = -t231;
t235 = -t234;
t237 = t13.*(t172+t6.*(t153+t5.*(t40-t41))).*(-2.307742028511301e-1);
t238 = t28.*t228.*6.983537755459679e-21;
t239 = t28.*t228.*1.396707551091936e-20;
t240 = t26.*t229.*6.983537755459679e-21;
t241 = t26.*t229.*1.396707551091936e-20;
t246 = t245.^2;
t247 = t244.^2;
t250 = -t248;
t251 = -t249;
t252 = t76+t77+t245;
t253 = t74+t75+t244;
t256 = t145+t146+t147+t148+t149+t220;
t257 = t140+t141+t142+t143+t144+t225;
t260 = t199.*t229.*1.731391914779077e-40;
t261 = t197.*t228.*1.731391914779077e-40;
t264 = t162.*t244.*2.117775986156328e-40;
t265 = t161.*t245.*2.117775986156328e-40;
t266 = t162.*t244.*1.058887993078164e-40;
t268 = t161.*t245.*1.058887993078164e-40;
t204 = t51+t196;
t205 = t53+t198;
t206 = t188.*1.675146132236819e-41;
t207 = t189.*1.675146132236819e-41;
t215 = t173+t175;
t236 = t7.*t212.*2.307742028511301e-1;
t242 = -t239;
t243 = -t238;
t254 = t247.*5.294439965390821e-41;
t255 = t246.*5.294439965390821e-41;
t258 = t55+t252;
t259 = t54+t253;
t262 = -t260;
t263 = -t261;
t267 = -t265;
t269 = -t268;
t276 = t123+t124+t125+t126+t256;
t277 = t118+t119+t120+t121+t257;
t270 = t161.*t258.*2.117775986156328e-40;
t271 = t162.*t259.*2.117775986156328e-40;
t272 = t161.*t258.*1.058887993078164e-40;
t274 = t162.*t259.*1.058887993078164e-40;
t278 = t277.^2;
t279 = t276.^2;
t281 = t82+t83+t277;
t283 = t84+t85+t276;
t286 = t186.*t276.*8.375730661184093e-42;
t287 = t186.*t276.*1.675146132236819e-41;
t288 = t187.*t277.*8.375730661184093e-42;
t289 = t187.*t277.*1.675146132236819e-41;
t292 = t244.*t259.*5.294439965390821e-41;
t293 = t245.*t258.*5.294439965390821e-41;
t273 = -t271;
t275 = -t274;
t280 = t279.*4.187865330592046e-42;
t282 = t278.*4.187865330592046e-42;
t284 = t61+t281;
t285 = t65+t283;
t290 = -t287;
t291 = -t286;
t294 = -t292;
t295 = -t293;
t296 = t186.*t285.*8.375730661184093e-42;
t297 = t186.*t285.*1.675146132236819e-41;
t298 = t187.*t284.*8.375730661184093e-42;
t299 = t187.*t284.*1.675146132236819e-41;
t302 = t277.*t284.*4.187865330592046e-42;
t303 = t276.*t285.*4.187865330592046e-42;
t300 = -t299;
t301 = -t298;
t304 = -t302;
t305 = -t303;
dL23_dq1 = [t178+t180+t190+t191+t208+t211+t236+t237-t18.*(t232+t235+t240+t243+t272+t275+t296+t301+t169.*t179.*7.099293763107699e-37-t168.*t181.*7.099293763107699e-37+t177.*t197.*3.462783829558153e-40-t176.*t199.*3.462783829558153e-40-t192.*t245.*1.058887993078164e-40+t193.*t244.*1.058887993078164e-40-t204.*t276.*8.375730661184093e-42+t205.*t277.*8.375730661184093e-42+t18.*(t20.*t168.*2.876800071864181e-18+t24.*t169.*2.876800071864181e-18+t26.*t176.*1.396707551091936e-20+t28.*t177.*1.396707551091936e-20+t161.*t192.*2.117775986156328e-40+t162.*t193.*2.117775986156328e-40+t186.*t204.*1.675146132236819e-41+t187.*t205.*1.675146132236819e-41+4.163565255585724));t178+t180+t190+t191+t208+t211+t236+t237+t250+t251+t262+t263+t294+t295+t304+t305+t179.*t203.*3.549646881553849e-37+t181.*t202.*3.549646881553849e-37+t197.*t217.*1.731391914779077e-40+t199.*t218.*1.731391914779077e-40+t244.*t253.*5.294439965390821e-41+t245.*t252.*5.294439965390821e-41+t277.*t281.*4.187865330592046e-42+t276.*t283.*4.187865330592046e-42+(t18.*(t230+t233+t241+t242+t270+t273+t297+t300-t20.*t202.*2.876800071864181e-18+t24.*t203.*2.876800071864181e-18-t26.*t218.*1.396707551091936e-20+t28.*t217.*1.396707551091936e-20-t161.*t252.*2.117775986156328e-40+t162.*t253.*2.117775986156328e-40+t187.*t281.*1.675146132236819e-41-t186.*t283.*1.675146132236819e-41))./2.0-t18.*(t232+t235+t240+t243+t272+t275+t296+t301+t158.*t179.*7.099293763107699e-37-t157.*t181.*7.099293763107699e-37+t160.*t197.*3.462783829558153e-40-t159.*t199.*3.462783829558153e-40+t182.*t244.*1.058887993078164e-40-t183.*t245.*1.058887993078164e-40-t196.*t276.*8.375730661184093e-42+t198.*t277.*8.375730661184093e-42+t18.*(t20.*t157.*2.876800071864181e-18+t24.*t158.*2.876800071864181e-18+t26.*t159.*1.396707551091936e-20+t28.*t160.*1.396707551091936e-20+t161.*t183.*2.117775986156328e-40+t162.*t182.*2.117775986156328e-40+t186.*t196.*1.675146132236819e-41+t187.*t198.*1.675146132236819e-41+4.163565255585724));t178+t180+t190+t191+t208+t211+t216+t219+t236+t237+t250+t251+t254+t255+t262+t263+t280+t282+t294+t295+t304+t305-t18.*(t221+t227+t232+t235+t240+t243+t266+t269+t272+t275+t288+t291+t296+t301+t18.*(t184+t185+t194+t195+t206+t207+t20.^2.*5.828733765381593+t24.^2.*5.828733765381593+4.163565255585724)-t20.*t181.*1.438400035932091e-18+t24.*t179.*1.438400035932091e-18)+(t18.*(t222+t226+t230+t233+t241+t242+t264+t267+t270+t273+t289+t290+t297+t300-t20.*t181.*2.876800071864181e-18+t24.*t179.*2.876800071864181e-18))./2.0+t179.^2.*3.549646881553849e-37+t181.^2.*3.549646881553849e-37;t190+t191+t208+t211+t216+t219+t236+t237+t254+t255+t262+t263+t280+t282+t294+t295+t304+t305-t18.*(t221+t227+t240+t243+t266+t269+t272+t275+t288+t291+t296+t301+t18.*(t184+t185+t194+t195+t206+t207+2.014638691514606e-1))+(t18.*(t222+t226+t241+t242+t264+t267+t270+t273+t289+t290+t297+t300))./2.0;t208+t211+t236+t237-(t18.*(t30.*t259.*4.41247130108735e-21-t32.*t258.*4.41247130108735e-21+t161.*t210.*2.117775986156328e-40-t162.*t209.*2.117775986156328e-40+t186.*t256.*1.675146132236819e-41-t187.*t257.*1.675146132236819e-41-t165.*t285.*1.675146132236819e-41+t166.*t284.*1.675146132236819e-41))./2.0+t209.*t244.*5.294439965390821e-41+t210.*t245.*5.294439965390821e-41-t209.*t259.*5.294439965390821e-41-t210.*t258.*5.294439965390821e-41+t256.*t276.*4.187865330592046e-42+t257.*t277.*4.187865330592046e-42-t256.*t285.*4.187865330592046e-42-t257.*t284.*4.187865330592046e-42-t18.*(t30.*t244.*2.206235650543675e-21-t32.*t245.*2.206235650543675e-21-t30.*t259.*2.206235650543675e-21+t32.*t258.*2.206235650543675e-21-t165.*t276.*8.375730661184093e-42+t166.*t277.*8.375730661184093e-42+t165.*t285.*8.375730661184093e-42-t166.*t284.*8.375730661184093e-42+t18.*(t30.*t162.*4.41247130108735e-21+t32.*t161.*4.41247130108735e-21+t165.*t186.*1.675146132236819e-41+t166.*t187.*1.675146132236819e-41+4.472206631138927e-2));t236+t237+(t18.*(t33.*t285.*4.178769965257346e-22-t34.*t284.*4.178769965257346e-22-t186.*t220.*1.675146132236819e-41+t187.*t225.*1.675146132236819e-41))./2.0+t220.*t276.*4.187865330592046e-42+t225.*t277.*4.187865330592046e-42-t220.*t285.*4.187865330592046e-42-t225.*t284.*4.187865330592046e-42-t18.*(t33.*t276.*(-2.089384982628673e-22)+t34.*t277.*2.089384982628673e-22+t33.*t285.*2.089384982628673e-22-t34.*t284.*2.089384982628673e-22+t18.*(t33.*t186.*4.178769965257346e-22+t34.*t187.*4.178769965257346e-22))];
