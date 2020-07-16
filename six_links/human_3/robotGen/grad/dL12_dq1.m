function dL12_dq1 = dL12_dq1(in1,in2,sampT)
%DL12_DQ1
%    DL12_DQ1 = DL12_DQ1(IN1,IN2,SAMPT)

%    This function was generated by the Symbolic Math Toolbox version 8.4.
%    05-Jul-2020 23:38:58

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
t15 = q_t2+q_t3;
t16 = qd_t1.^2;
t21 = 1.0./sampT;
t17 = cos(t14);
t18 = cos(t15);
t19 = sin(t14);
t20 = sin(t15);
t22 = t2.*t3;
t23 = t2.*t9;
t24 = t3.*t8;
t25 = q_t3+q_t4+t14;
t26 = t8.*t9;
t43 = t2.*2.547461129222121e+19;
t44 = t2.*5.094922258444242e+19;
t45 = t8.*2.547461129222121e+19;
t46 = t8.*5.094922258444242e+19;
t47 = t2.*1.018984451688848e+20;
t48 = t8.*1.018984451688848e+20;
t145 = qd_t1.*t9.*8.68583195846977;
t151 = t4.*1.514284177462349e+1;
t27 = cos(t25);
t28 = q_t5+t25;
t29 = sin(t25);
t32 = -t26;
t36 = t23+t24;
t49 = qd_t1.*t43;
t50 = qd_t1.*t44;
t51 = qd_t1.*t45;
t52 = qd_t1.*t46;
t53 = t17.*2.465788349979757e+19;
t54 = t19.*2.465788349979757e+19;
t55 = qd_t1.*t47;
t56 = qd_t1.*t48;
t57 = t17.*4.931576699959514e+19;
t58 = t19.*4.931576699959514e+19;
t63 = t17.*9.863153399919029e+19;
t64 = t19.*9.863153399919029e+19;
t149 = t23.*6.30747417130722e+1;
t150 = t24.*6.30747417130722e+1;
t152 = qd_t1.*t20.*7.822204368661843;
t153 = qd_t1.*t20.*3.911102184330921;
t154 = qd_t1.*qd_t2.*t18.*1.95555109216546;
t155 = qd_t1.*qd_t3.*t18.*1.95555109216546;
t156 = t16.*t18.*1.95555109216546;
t30 = q_t6+t28;
t31 = sin(t28);
t33 = cos(t28);
t37 = t22+t32;
t38 = t4.*t36;
t39 = t10.*t36;
t59 = qd_t1.*t57;
t60 = qd_t2.*t57;
t61 = qd_t1.*t58;
t62 = qd_t2.*t58;
t65 = qd_t1.*t53;
t66 = qd_t2.*t53;
t67 = qd_t1.*t54;
t68 = qd_t2.*t54;
t69 = qd_t1.*t63;
t70 = qd_t2.*t63;
t71 = qd_t1.*t64;
t72 = qd_t2.*t64;
t73 = t27.*2.465788349979757e+19;
t74 = t29.*2.465788349979757e+19;
t75 = t27.*2.13537271108247e+19;
t76 = t29.*2.13537271108247e+19;
t81 = t27.*9.863153399919029e+19;
t82 = t29.*9.863153399919029e+19;
t157 = -t154;
t158 = -t155;
t159 = -t156;
t34 = cos(t30);
t35 = sin(t30);
t40 = t4.*t37;
t41 = t10.*t37;
t42 = -t39;
t77 = qd_t1.*t76;
t78 = qd_t2.*t76;
t79 = qd_t3.*t76;
t80 = qd_t4.*t76;
t83 = t33.*1.103050668953178e+19;
t84 = qd_t1.*t73;
t85 = qd_t2.*t73;
t86 = qd_t3.*t73;
t87 = qd_t4.*t73;
t88 = t31.*1.103050668953178e+19;
t89 = qd_t1.*t74;
t90 = qd_t2.*t74;
t91 = qd_t3.*t74;
t92 = qd_t4.*t74;
t93 = qd_t1.*t75;
t94 = qd_t2.*t75;
t95 = qd_t3.*t75;
t96 = qd_t4.*t75;
t97 = qd_t1.*t82;
t98 = qd_t2.*t82;
t99 = qd_t3.*t82;
t100 = qd_t4.*t82;
t113 = t33.*1.018984451688848e+20;
t114 = t31.*1.018984451688848e+20;
t115 = qd_t1.*t81;
t116 = qd_t2.*t81;
t117 = qd_t3.*t81;
t118 = qd_t4.*t81;
t160 = t38.*4.23936095954724e+1;
t162 = t57+t75;
t163 = t58+t76;
t101 = t34.*2.641309840417044e+19;
t102 = qd_t1.*t83;
t103 = qd_t2.*t83;
t104 = qd_t3.*t83;
t105 = qd_t4.*t83;
t106 = qd_t5.*t83;
t107 = t35.*2.641309840417044e+19;
t108 = qd_t1.*t88;
t109 = qd_t2.*t88;
t110 = qd_t3.*t88;
t111 = qd_t4.*t88;
t112 = qd_t5.*t88;
t131 = qd_t1.*t113;
t132 = qd_t2.*t113;
t133 = qd_t3.*t113;
t134 = qd_t4.*t113;
t135 = qd_t5.*t113;
t136 = qd_t1.*t114;
t137 = qd_t2.*t114;
t138 = qd_t3.*t114;
t139 = qd_t4.*t114;
t140 = qd_t5.*t114;
t141 = t38+t41;
t142 = t40+t42;
t146 = -t5.*(t39-t40);
t147 = -t11.*(t39-t40);
t161 = t41.*4.23936095954724e+1;
t164 = t73+t83;
t165 = t74+t88;
t169 = t11.*(t39-t40).*(-9.5320875106278);
t177 = t46+t163;
t178 = t44+t162;
t188 = t93+t94+t95+t96;
t189 = t77+t78+t79+t80;
t193 = t27.*t162.*1.788951796939387e-20;
t194 = t29.*t163.*1.788951796939387e-20;
t119 = qd_t1.*t101;
t120 = qd_t2.*t101;
t121 = qd_t3.*t101;
t122 = qd_t4.*t101;
t123 = qd_t5.*t101;
t124 = qd_t6.*t101;
t125 = qd_t1.*t107;
t126 = qd_t2.*t107;
t127 = qd_t3.*t107;
t128 = qd_t4.*t107;
t129 = qd_t5.*t107;
t130 = qd_t6.*t107;
t143 = t5.*t141;
t144 = t11.*t141;
t166 = t101+t113;
t167 = t107+t114;
t179 = t53+t164;
t180 = t54+t165;
t195 = t102+t103+t104+t105+t106;
t196 = t108+t109+t110+t111+t112;
t199 = t61+t62+t189;
t200 = t59+t60+t188;
t209 = t163.*t188.*8.377702813447146e-40;
t210 = t162.*t189.*8.377702813447146e-40;
t148 = -t144;
t168 = t143.*9.5320875106278;
t170 = t143+t147;
t174 = -t6.*(t144+t5.*(t39-t40));
t175 = -t12.*(t144+t5.*(t39-t40));
t181 = t81+t166;
t182 = t82+t167;
t183 = t43+t179;
t184 = t45+t180;
t186 = t12.*(t144+t5.*(t39-t40)).*(-2.70089708993775);
t201 = t125+t126+t127+t128+t129+t130;
t202 = t119+t120+t121+t122+t123+t124;
t204 = t50+t200;
t205 = t52+t199;
t207 = t164.*t179.*1.024730315882094e-39;
t208 = t165.*t180.*1.024730315882094e-39;
t211 = t29.*t200.*8.944758984696935e-21;
t212 = t27.*t199.*8.944758984696935e-21;
t213 = t84+t85+t86+t87+t195;
t214 = t89+t90+t91+t92+t196;
t232 = t188.*t200.*2.094425703361786e-40;
t233 = t189.*t199.*2.094425703361786e-40;
t171 = t146+t148;
t172 = t6.*t170;
t173 = t12.*t170;
t187 = t64+t182;
t190 = t63+t181;
t215 = -t211;
t218 = t29.*t204.*8.944758984696935e-21;
t219 = t29.*t204.*1.788951796939387e-20;
t220 = t27.*t205.*8.944758984696935e-21;
t221 = t27.*t205.*1.788951796939387e-20;
t223 = t67+t68+t214;
t224 = t65+t66+t213;
t225 = t163.*t204.*4.188851406723573e-40;
t226 = t162.*t205.*4.188851406723573e-40;
t228 = t136+t137+t138+t139+t140+t201;
t229 = t131+t132+t133+t134+t135+t202;
t234 = t189.*t205.*2.094425703361786e-40;
t235 = t188.*t204.*2.094425703361786e-40;
t238 = t180.*t213.*1.024730315882094e-39;
t239 = t179.*t214.*1.024730315882094e-39;
t176 = -t173;
t185 = t172.*2.70089708993775;
t191 = t47+t190;
t192 = t48+t187;
t197 = t172+t175;
t206 = t13.*(t173+t6.*(t144+t5.*(t39-t40))).*(-2.955836564032501e-1);
t216 = t181.*t190.*2.026386450286474e-41;
t217 = t182.*t187.*2.026386450286474e-41;
t222 = -t220;
t227 = -t226;
t230 = t51+t223;
t231 = t49+t224;
t236 = -t234;
t237 = -t235;
t240 = t164.*t223.*5.123651579410472e-40;
t241 = t165.*t224.*5.123651579410472e-40;
t243 = t97+t98+t99+t100+t228;
t244 = t115+t116+t117+t118+t229;
t259 = t213.*t224.*2.561825789705236e-40;
t260 = t214.*t223.*2.561825789705236e-40;
t198 = t174+t176;
t203 = t7.*t197.*2.955836564032501e-1;
t242 = -t241;
t245 = t165.*t231.*1.024730315882094e-39;
t246 = t165.*t231.*5.123651579410472e-40;
t247 = t164.*t230.*1.024730315882094e-39;
t248 = t164.*t230.*5.123651579410472e-40;
t250 = t179.*t230.*5.123651579410472e-40;
t251 = t180.*t231.*5.123651579410472e-40;
t253 = t69+t70+t244;
t254 = t71+t72+t243;
t257 = t190.*t243.*2.026386450286474e-41;
t258 = t187.*t244.*2.026386450286474e-41;
t264 = t214.*t230.*2.561825789705236e-40;
t265 = t213.*t231.*2.561825789705236e-40;
t249 = -t248;
t252 = -t250;
t255 = t55+t253;
t256 = t56+t254;
t261 = t182.*t253.*1.013193225143237e-41;
t262 = t181.*t254.*1.013193225143237e-41;
t266 = -t264;
t267 = -t265;
t276 = t243.*t254.*5.065966125716185e-42;
t277 = t244.*t253.*5.065966125716185e-42;
t263 = -t261;
t268 = t181.*t256.*2.026386450286474e-41;
t269 = t182.*t255.*2.026386450286474e-41;
t270 = t181.*t256.*1.013193225143237e-41;
t271 = t182.*t255.*1.013193225143237e-41;
t273 = t190.*t256.*1.013193225143237e-41;
t274 = t187.*t255.*1.013193225143237e-41;
t278 = t244.*t255.*5.065966125716185e-42;
t279 = t243.*t256.*5.065966125716185e-42;
t272 = -t270;
t275 = -t273;
t280 = -t278;
t281 = -t279;
dL12_dq1 = [t149+t150+t160+t161+t168+t169+t185+t186+t203+t206+t21.*(t145+t225+t227+t251+t252+t274+t275+t21.*(t3.*8.68583195846977+t18.*7.822204368661843+t151+t162.*t178.*8.377702813447146e-40+t163.*t177.*8.377702813447146e-40+t179.*t183.*1.024730315882094e-39+t180.*t184.*1.024730315882094e-39+t187.*t192.*2.026386450286474e-41+t190.*t191.*2.026386450286474e-41+2.167104011981295e+1)+qd_t2.*t9.*4.342915979234885+qd_t1.*t20.*7.822204368661842+qd_t2.*t20.*3.911102184330921+qd_t3.*t20.*3.911102184330921-t177.*t200.*4.188851406723573e-40+t178.*t199.*4.188851406723573e-40+t183.*t223.*5.123651579410472e-40-t184.*t224.*5.123651579410472e-40+t191.*t254.*1.013193225143237e-41-t192.*t253.*1.013193225143237e-41);t149+t150+t157+t158+t159+t160+t161+t168+t169+t185+t186+t203+t206-t3.*t16.*2.171457989617442-t199.*t205.*2.094425703361786e-40-t200.*t204.*2.094425703361786e-40-t223.*t230.*2.561825789705236e-40-t224.*t231.*2.561825789705236e-40-t253.*t255.*5.065966125716185e-42-t254.*t256.*5.065966125716185e-42+t21.*(t153+t225+t227+t251+t252+t274+t275+qd_t1.*t9.*4.342915979234885+t162.*t199.*4.188851406723573e-40-t163.*t200.*4.188851406723573e-40+t179.*t223.*5.123651579410472e-40-t180.*t224.*5.123651579410472e-40-t187.*t253.*1.013193225143237e-41+t190.*t254.*1.013193225143237e-41+t21.*(t151+t162.^2.*8.377702813447146e-40+t163.^2.*8.377702813447146e-40+t179.^2.*1.024730315882094e-39+t180.^2.*1.024730315882094e-39+t187.^2.*2.026386450286474e-41+t190.^2.*2.026386450286474e-41+2.167104011981295e+1))+t199.^2.*2.094425703361786e-40+t200.^2.*2.094425703361786e-40+t223.^2.*2.561825789705236e-40+t224.^2.*2.561825789705236e-40+t253.^2.*5.065966125716185e-42+t254.^2.*5.065966125716185e-42+(t21.*(t145+t152+t162.*t199.*8.377702813447146e-40-t163.*t200.*8.377702813447146e-40-t162.*t205.*8.377702813447146e-40+t163.*t204.*8.377702813447146e-40+t179.*t223.*1.024730315882094e-39-t180.*t224.*1.024730315882094e-39-t179.*t230.*1.024730315882094e-39+t180.*t231.*1.024730315882094e-39-t187.*t253.*2.026386450286474e-41+t187.*t255.*2.026386450286474e-41+t190.*t254.*2.026386450286474e-41-t190.*t256.*2.026386450286474e-41))./2.0-qd_t1.*qd_t2.*t3.*2.171457989617442;t157+t158+t159+t160+t161+t168+t169+t185+t186+t203+t206+t232+t233+t236+t237+t259+t260+t266+t267+t276+t277+t280+t281+(t21.*(t152-t209+t210+t219-t221-t238+t239+t245-t247+t257-t258-t268+t269+qd_t1.*t10.*1.514284177462349e+1+qd_t2.*t10.*1.514284177462349e+1+qd_t3.*t10.*7.571420887311746))./2.0+t21.*(t153+t212+t215+t218+t222+t240+t242+t246+t249+t262+t263+t271+t272+t21.*(t4.*7.571420887311746+t193+t194+t207+t208+t216+t217+1.355133754936409e+1));t168+t169+t185+t186+t203+t206+t232+t233+t236+t237+t259+t260+t266+t267+t276+t277+t280+t281+t21.*(t212+t215+t218+t222+t240+t242+t246+t249+t262+t263+t271+t272+t21.*(t193+t194+t207+t208+t216+t217+2.732208963266266e-1))-(t21.*(t209-t210-t219+t221+t238-t239-t245+t247-t257+t258+t268-t269))./2.0;t185+t186+t203+t206+(t21.*(t31.*t231.*1.130329460430345e-20-t33.*t230.*1.130329460430345e-20+t179.*t196.*1.024730315882094e-39-t180.*t195.*1.024730315882094e-39-t187.*t229.*2.026386450286474e-41+t190.*t228.*2.026386450286474e-41-t166.*t256.*2.026386450286474e-41+t167.*t255.*2.026386450286474e-41))./2.0+t195.*t224.*2.561825789705236e-40+t196.*t223.*2.561825789705236e-40-t195.*t231.*2.561825789705236e-40-t196.*t230.*2.561825789705236e-40+t228.*t254.*5.065966125716185e-42+t229.*t253.*5.065966125716185e-42-t228.*t256.*5.065966125716185e-42-t229.*t255.*5.065966125716185e-42+t21.*(t31.*t224.*(-5.651647302151728e-21)+t33.*t223.*5.651647302151728e-21+t31.*t231.*5.651647302151728e-21-t33.*t230.*5.651647302151728e-21+t166.*t254.*1.013193225143237e-41-t167.*t253.*1.013193225143237e-41-t166.*t256.*1.013193225143237e-41+t167.*t255.*1.013193225143237e-41+t21.*(t31.*t180.*1.130329460430345e-20+t33.*t179.*1.130329460430345e-20+t167.*t187.*2.026386450286474e-41+t166.*t190.*2.026386450286474e-41+6.06510889254805e-2));t203+t206+t201.*t254.*5.065966125716185e-42+t202.*t253.*5.065966125716185e-42-t201.*t256.*5.065966125716185e-42-t202.*t255.*5.065966125716185e-42-(t21.*(t34.*t256.*5.352314471629427e-22-t35.*t255.*5.352314471629427e-22+t187.*t202.*2.026386450286474e-41-t190.*t201.*2.026386450286474e-41))./2.0+t21.*(t34.*t254.*2.676157235814713e-22-t35.*t253.*2.676157235814713e-22-t34.*t256.*2.676157235814713e-22+t35.*t255.*2.676157235814713e-22+t21.*(t35.*t187.*5.352314471629427e-22+t34.*t190.*5.352314471629427e-22))];
