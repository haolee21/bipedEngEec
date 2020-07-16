function dL14_dq1 = dL14_dq1(in1,in2,sampT)
%DL14_DQ1
%    DL14_DQ1 = DL14_DQ1(IN1,IN2,SAMPT)

%    This function was generated by the Symbolic Math Toolbox version 8.4.
%    14-Jul-2020 22:54:16

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
t41 = t2.*8.287740207069301e+15;
t42 = t8.*8.287740207069301e+15;
t43 = t8.*1.65754804141386e+16;
t46 = t2.*3.31509608282772e+16;
t47 = t8.*3.31509608282772e+16;
t48 = t2.*1.65754804141386e+16;
t23 = cos(t21);
t24 = q_t5+t21;
t25 = sin(t21);
t28 = -t22;
t34 = t19+t20;
t44 = qd_t1.*t41;
t45 = qd_t1.*t42;
t49 = qd_t1.*t46;
t50 = qd_t1.*t47;
t51 = qd_t1.*t48;
t52 = qd_t1.*t43;
t53 = t15.*8.022031431934143e+15;
t54 = t16.*8.022031431934143e+15;
t55 = t16.*3.208812572773657e+16;
t60 = t15.*1.604406286386829e+16;
t61 = t16.*1.604406286386829e+16;
t62 = t15.*3.208812572773657e+16;
t26 = q_t6+t24;
t27 = sin(t24);
t29 = cos(t24);
t32 = t23.^2;
t33 = t25.^2;
t35 = t18+t28;
t36 = t4.*t34;
t37 = t10.*t34;
t56 = qd_t1.*t53;
t57 = qd_t2.*t53;
t58 = qd_t1.*t54;
t59 = qd_t2.*t54;
t63 = qd_t1.*t61;
t64 = qd_t2.*t61;
t65 = t23.*3.473539610027484e+15;
t66 = qd_t1.*t62;
t67 = qd_t2.*t62;
t68 = t25.*3.473539610027484e+15;
t69 = qd_t1.*t55;
t70 = qd_t2.*t55;
t71 = qd_t1.*t60;
t72 = qd_t2.*t60;
t73 = t23.*1.604406286386829e+16;
t74 = t25.*1.604406286386829e+16;
t81 = t23.*3.208812572773657e+16;
t86 = t25.*3.208812572773657e+16;
t30 = cos(t26);
t31 = sin(t26);
t38 = t4.*t35;
t39 = t10.*t35;
t40 = -t37;
t75 = t29.*7.177183019322015e+15;
t76 = t27.*7.177183019322015e+15;
t77 = qd_t1.*t65;
t78 = qd_t2.*t65;
t79 = qd_t3.*t65;
t80 = qd_t4.*t65;
t82 = qd_t1.*t68;
t83 = qd_t2.*t68;
t84 = qd_t3.*t68;
t85 = qd_t4.*t68;
t87 = qd_t1.*t81;
t88 = qd_t2.*t81;
t89 = qd_t3.*t81;
t90 = qd_t4.*t81;
t91 = qd_t1.*t86;
t92 = qd_t2.*t86;
t93 = qd_t3.*t86;
t94 = qd_t4.*t86;
t96 = qd_t1.*t73;
t97 = qd_t2.*t73;
t98 = qd_t3.*t73;
t99 = qd_t4.*t73;
t101 = qd_t1.*t74;
t102 = qd_t2.*t74;
t103 = qd_t3.*t74;
t104 = qd_t4.*t74;
t105 = t29.*3.31509608282772e+16;
t111 = t27.*3.31509608282772e+16;
t146 = t53+t65;
t147 = t54+t68;
t167 = t32.*3.948475943260698e-1;
t168 = t33.*3.948475943260698e-1;
t95 = t30.*8.593061347490117e+15;
t100 = t31.*8.593061347490117e+15;
t106 = qd_t1.*t75;
t107 = qd_t2.*t75;
t108 = qd_t3.*t75;
t109 = qd_t4.*t75;
t110 = qd_t5.*t75;
t112 = qd_t1.*t76;
t113 = qd_t2.*t76;
t114 = qd_t3.*t76;
t115 = qd_t4.*t76;
t116 = qd_t5.*t76;
t129 = qd_t1.*t105;
t130 = qd_t2.*t105;
t131 = qd_t3.*t105;
t132 = qd_t4.*t105;
t133 = qd_t5.*t105;
t134 = qd_t1.*t111;
t135 = qd_t2.*t111;
t136 = qd_t3.*t111;
t137 = qd_t4.*t111;
t138 = qd_t5.*t111;
t139 = t36+t39;
t140 = t38+t40;
t143 = -t5.*(t37-t38);
t144 = -t11.*(t37-t38);
t148 = t73+t75;
t149 = t74+t76;
t154 = t41+t146;
t155 = t42+t147;
t176 = t77+t78+t79+t80;
t177 = t82+t83+t84+t85;
t182 = t11.*(t37-t38).*(-9.690955635804931);
t117 = qd_t1.*t95;
t118 = qd_t2.*t95;
t119 = qd_t3.*t95;
t120 = qd_t4.*t95;
t121 = qd_t5.*t95;
t122 = qd_t6.*t95;
t123 = qd_t1.*t100;
t124 = qd_t2.*t100;
t125 = qd_t3.*t100;
t126 = qd_t4.*t100;
t127 = qd_t5.*t100;
t128 = qd_t6.*t100;
t141 = t5.*t139;
t142 = t11.*t139;
t150 = t95+t105;
t151 = t100+t111;
t152 = t148.^2;
t153 = t149.^2;
t156 = t60+t148;
t158 = t61+t149;
t179 = t177.^2;
t180 = t176.^2;
t187 = t106+t107+t108+t109+t110;
t188 = t112+t113+t114+t115+t116;
t192 = t56+t57+t176;
t193 = t58+t59+t177;
t198 = t25.*t176.*1.136729787638568e-16;
t199 = t23.*t177.*1.136729787638568e-16;
t202 = t25.*t176.*5.683648938192842e-17;
t203 = t23.*t177.*5.683648938192842e-17;
t145 = -t142;
t157 = t141+t144;
t162 = t86+t151;
t163 = t81+t150;
t164 = -t6.*(t142+t5.*(t37-t38));
t165 = -t12.*(t142+t5.*(t37-t38));
t171 = t48+t156;
t172 = t43+t158;
t173 = t153.*2.501782997759019e-33;
t174 = t152.*2.501782997759019e-33;
t181 = t141.*9.690955635804931;
t190 = t180.*8.181350403756978e-33;
t191 = t179.*8.181350403756978e-33;
t194 = t12.*(t142+t5.*(t37-t38)).*(-2.745912041436713);
t197 = t123+t124+t125+t126+t127+t128;
t200 = t117+t118+t119+t120+t121+t122;
t201 = -t198;
t204 = -t202;
t205 = t44+t192;
t206 = t45+t193;
t215 = t96+t97+t98+t99+t187;
t216 = t101+t102+t103+t104+t188;
t159 = t143+t145;
t160 = t6.*t157;
t161 = t12.*t157;
t169 = t163.^2;
t170 = t162.^2;
t175 = t55+t162;
t178 = t62+t163;
t208 = t25.*t205.*1.136729787638568e-16;
t209 = t23.*t206.*1.136729787638568e-16;
t211 = t25.*t205.*5.683648938192842e-17;
t213 = t23.*t206.*5.683648938192842e-17;
t217 = t216.^2;
t218 = t215.^2;
t221 = t71+t72+t215;
t222 = t63+t64+t216;
t223 = t134+t135+t136+t137+t138+t197;
t224 = t129+t130+t131+t132+t133+t200;
t227 = t177.*t206.*8.181350403756978e-33;
t228 = t176.*t205.*8.181350403756978e-33;
t231 = t149.*t215.*2.501782997759019e-33;
t232 = t148.*t216.*2.501782997759019e-33;
t234 = t149.*t215.*1.25089149887951e-33;
t235 = t148.*t216.*1.25089149887951e-33;
t166 = -t161;
t183 = t46+t178;
t184 = t47+t175;
t185 = t169.*1.978893017857885e-34;
t186 = t170.*1.978893017857885e-34;
t189 = t160.*2.745912041436713;
t195 = t160+t165;
t210 = t13.*(t161+t6.*(t142+t5.*(t37-t38))).*(-3.005100506766376e-1);
t212 = -t209;
t214 = -t213;
t219 = t218.*6.254457494397548e-34;
t220 = t217.*6.254457494397548e-34;
t225 = t52+t222;
t226 = t51+t221;
t229 = -t227;
t230 = -t228;
t233 = -t231;
t236 = -t234;
t243 = t91+t92+t93+t94+t223;
t244 = t87+t88+t89+t90+t224;
t196 = t164+t166;
t207 = t7.*t195.*3.005100506766376e-1;
t237 = t148.*t225.*2.501782997759019e-33;
t238 = t149.*t226.*2.501782997759019e-33;
t240 = t148.*t225.*1.25089149887951e-33;
t241 = t149.*t226.*1.25089149887951e-33;
t245 = t243.^2;
t246 = t244.^2;
t249 = t66+t67+t244;
t250 = t69+t70+t243;
t253 = t163.*t243.*9.894465089289424e-35;
t254 = t163.*t243.*1.978893017857885e-34;
t255 = t162.*t244.*9.894465089289424e-35;
t256 = t162.*t244.*1.978893017857885e-34;
t259 = t215.*t226.*6.254457494397548e-34;
t260 = t216.*t225.*6.254457494397548e-34;
t269 = t167+t168+t173+t174+t185+t186+2.824041542309381e-1;
t239 = -t237;
t242 = -t240;
t247 = t245.*4.947232544644712e-35;
t248 = t246.*4.947232544644712e-35;
t251 = t49+t249;
t252 = t50+t250;
t257 = -t256;
t258 = -t255;
t261 = -t259;
t262 = -t260;
t270 = t17.*t269;
t263 = t163.*t252.*9.894465089289424e-35;
t264 = t163.*t252.*1.978893017857885e-34;
t265 = t162.*t251.*9.894465089289424e-35;
t266 = t162.*t251.*1.978893017857885e-34;
t271 = t243.*t252.*4.947232544644712e-35;
t272 = t244.*t251.*4.947232544644712e-35;
t267 = -t264;
t268 = -t263;
t273 = -t271;
t274 = -t272;
t276 = t17.*(t198-t199-t208+t209+t231-t232+t237-t238-t254+t256+t264-t266).*(-1.0./2.0);
t275 = t199+t201+t208+t212+t232+t233+t238+t239+t254+t257+t266+t267;
t277 = t203+t204+t211+t214+t235+t236+t241+t242+t253+t258+t265+t268+t270;
t278 = t17.*t277;
t279 = t181+t182+t189+t190+t191+t194+t207+t210+t219+t220+t229+t230+t247+t248+t261+t262+t273+t274+t276+t278;
dL14_dq1 = [t181+t182+t189+t194+t207+t210+t17.*(t211+t214+t241+t242+t265+t268+t154.*t177.*1.636270080751396e-32-t155.*t176.*1.636270080751396e-32+t171.*t216.*1.25089149887951e-33-t172.*t215.*1.25089149887951e-33+t183.*t243.*9.894465089289424e-35-t184.*t244.*9.894465089289424e-35+t17.*(t23.*t154.*1.136729787638568e-16+t25.*t155.*1.136729787638568e-16+t148.*t171.*2.501782997759019e-33+t149.*t172.*2.501782997759019e-33+t162.*t184.*1.978893017857885e-34+t163.*t183.*1.978893017857885e-34+2.824041542309381e-1));t181+t182+t189+t194+t207+t210+t229+t230+t261+t262+t273+t274+t17.*(t211+t214+t241+t242+t265+t268+t146.*t177.*1.636270080751396e-32-t147.*t176.*1.636270080751396e-32+t156.*t216.*1.25089149887951e-33-t158.*t215.*1.25089149887951e-33-t175.*t244.*9.894465089289424e-35+t178.*t243.*9.894465089289424e-35+t17.*(t23.*t146.*1.136729787638568e-16+t25.*t147.*1.136729787638568e-16+t148.*t156.*2.501782997759019e-33+t149.*t158.*2.501782997759019e-33+t162.*t175.*1.978893017857885e-34+t163.*t178.*1.978893017857885e-34+2.824041542309381e-1))+t176.*t192.*8.181350403756978e-33+t177.*t193.*8.181350403756978e-33+t215.*t221.*6.254457494397548e-34+t216.*t222.*6.254457494397548e-34+t243.*t250.*4.947232544644712e-35+t244.*t249.*4.947232544644712e-35+(t17.*(t208+t212+t238+t239+t266+t267+t23.*t193.*1.136729787638568e-16-t25.*t192.*1.136729787638568e-16+t148.*t222.*2.501782997759019e-33-t149.*t221.*2.501782997759019e-33-t162.*t249.*1.978893017857885e-34+t163.*t250.*1.978893017857885e-34))./2.0;t279;t279;t189+t194+t207+t210+(t17.*(t27.*t226.*1.795575444954456e-17-t29.*t225.*1.795575444954456e-17+t148.*t188.*2.501782997759019e-33-t149.*t187.*2.501782997759019e-33-t162.*t224.*1.978893017857885e-34+t163.*t223.*1.978893017857885e-34-t150.*t252.*1.978893017857885e-34+t151.*t251.*1.978893017857885e-34))./2.0+t187.*t215.*6.254457494397548e-34+t188.*t216.*6.254457494397548e-34-t187.*t226.*6.254457494397548e-34-t188.*t225.*6.254457494397548e-34+t223.*t243.*4.947232544644712e-35+t224.*t244.*4.947232544644712e-35-t223.*t252.*4.947232544644712e-35-t224.*t251.*4.947232544644712e-35+t17.*(t27.*t215.*(-8.97787722477228e-18)+t29.*t216.*8.97787722477228e-18+t27.*t226.*8.97787722477228e-18-t29.*t225.*8.97787722477228e-18+t150.*t243.*9.894465089289424e-35-t151.*t244.*9.894465089289424e-35-t150.*t252.*9.894465089289424e-35+t151.*t251.*9.894465089289424e-35+t17.*(t27.*t149.*1.795575444954456e-17+t29.*t148.*1.795575444954456e-17+t150.*t163.*1.978893017857885e-34+t151.*t162.*1.978893017857885e-34+6.26896394143647e-2));t207+t210-(t17.*(t30.*t252.*1.700474910257266e-18-t31.*t251.*1.700474910257266e-18-t163.*t197.*1.978893017857885e-34+t162.*t200.*1.978893017857885e-34))./2.0+t17.*(t30.*t243.*8.50237455128633e-19-t31.*t244.*8.50237455128633e-19-t30.*t252.*8.50237455128633e-19+t31.*t251.*8.50237455128633e-19+t17.*(t30.*t163.*1.700474910257266e-18+t31.*t162.*1.700474910257266e-18))+t197.*t243.*4.947232544644712e-35+t200.*t244.*4.947232544644712e-35-t197.*t252.*4.947232544644712e-35-t200.*t251.*4.947232544644712e-35];
