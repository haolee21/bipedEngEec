function dV_dx = dV_dx(qd1,qd2,qd3,qd4,qd5,qd6,q2,q3,q4,q5,q6)
%DV_DX
%    DV_DX = DV_DX(QD1,QD2,QD3,QD4,QD5,QD6,Q2,Q3,Q4,Q5,Q6)

%    This function was generated by the Symbolic Math Toolbox version 8.4.
%    27-May-2020 01:00:05

t2 = cos(q2);
t3 = cos(q3);
t4 = cos(q5);
t5 = cos(q6);
t6 = sin(q2);
t7 = sin(q3);
t8 = sin(q5);
t9 = sin(q6);
t10 = q2+q3;
t11 = q3+q4;
t12 = q5+q6;
t13 = qd1.^2;
t14 = qd2.^2;
t15 = qd3.^2;
t16 = qd4.^2;
t17 = qd5.^2;
t18 = qd6.^2;
t19 = cos(t10);
t20 = cos(t11);
t21 = cos(t12);
t22 = q4+t10;
t23 = q5+t11;
t24 = sin(t10);
t25 = sin(t11);
t26 = sin(t12);
t30 = t11+t12;
t55 = qd1.*t9.*1.09078504542792e-1;
t56 = qd2.*t9.*1.09078504542792e-1;
t57 = qd3.*t9.*1.09078504542792e-1;
t58 = qd4.*t9.*1.09078504542792e-1;
t59 = qd5.*t9.*1.09078504542792e-1;
t60 = qd6.*t9.*1.09078504542792e-1;
t63 = qd1.*qd6.*t5.*1.09078504542792e-1;
t64 = qd2.*qd6.*t5.*1.09078504542792e-1;
t65 = qd3.*qd6.*t5.*1.09078504542792e-1;
t66 = qd4.*qd6.*t5.*1.09078504542792e-1;
t67 = qd5.*qd6.*t5.*1.09078504542792e-1;
t68 = qd1.*t6.*2.327631569418914e+1;
t69 = qd2.*t6.*2.327631569418914e+1;
t82 = t5.*t18.*5.4539252271396e-2;
t203 = qd1.*t7.*1.514284177462349e+1;
t204 = qd2.*t7.*1.514284177462349e+1;
t205 = qd3.*t7.*1.514284177462349e+1;
t206 = qd1.*qd3.*t3.*1.514284177462349e+1;
t207 = qd2.*qd3.*t3.*1.514284177462349e+1;
t211 = t3.*t15.*7.571420887311746;
t331 = qd1.*t8.*9.647505289768104e-1;
t332 = qd2.*t8.*9.647505289768104e-1;
t333 = qd3.*t8.*9.647505289768104e-1;
t334 = qd4.*t8.*9.647505289768104e-1;
t335 = qd5.*t8.*9.647505289768104e-1;
t345 = qd1.*qd5.*t4.*9.647505289768104e-1;
t346 = qd2.*qd5.*t4.*9.647505289768104e-1;
t347 = qd3.*qd5.*t4.*9.647505289768104e-1;
t348 = qd4.*qd5.*t4.*9.647505289768104e-1;
t349 = t4.*t17.*4.823752644884052e-1;
t27 = cos(t22);
t28 = cos(t23);
t29 = q5+t22;
t31 = sin(t22);
t32 = sin(t23);
t34 = cos(t30);
t35 = t12+t22;
t37 = sin(t30);
t40 = qd1.*t24.*1.564440873732369e+1;
t41 = qd2.*t24.*1.564440873732369e+1;
t42 = qd3.*t24.*1.564440873732369e+1;
t43 = qd1.*qd2.*t19.*1.564440873732369e+1;
t44 = qd1.*qd3.*t19.*1.564440873732369e+1;
t45 = qd2.*qd3.*t19.*1.564440873732369e+1;
t49 = t13.*t19.*7.822204368661843;
t50 = t14.*t19.*7.822204368661843;
t51 = t15.*t19.*7.822204368661843;
t70 = -t55;
t71 = -t56;
t72 = -t57;
t73 = -t58;
t74 = -t59;
t75 = -t60;
t76 = -t63;
t77 = -t64;
t78 = -t65;
t79 = -t66;
t80 = -t67;
t81 = -t69;
t87 = -t82;
t208 = -t203;
t209 = -t204;
t210 = -t205;
t212 = -t206;
t213 = -t207;
t214 = -t211;
t215 = qd1.*t25.*3.404826678658597;
t216 = qd2.*t25.*3.404826678658597;
t217 = qd3.*t25.*3.404826678658597;
t218 = qd4.*t25.*3.404826678658597;
t219 = qd1.*t26.*1.055813973565752e-1;
t220 = qd2.*t26.*1.055813973565752e-1;
t221 = qd3.*t26.*1.055813973565752e-1;
t222 = qd4.*t26.*1.055813973565752e-1;
t223 = qd5.*t26.*1.055813973565752e-1;
t224 = qd6.*t26.*1.055813973565752e-1;
t235 = qd1.*qd2.*t20.*3.404826678658597;
t236 = qd1.*qd3.*t20.*3.404826678658597;
t237 = qd1.*qd4.*t20.*3.404826678658597;
t238 = qd2.*qd3.*t20.*3.404826678658597;
t239 = qd2.*qd4.*t20.*3.404826678658597;
t240 = qd3.*qd4.*t20.*3.404826678658597;
t241 = qd1.*qd2.*t21.*1.055813973565752e-1;
t242 = qd1.*qd3.*t21.*1.055813973565752e-1;
t243 = qd1.*qd4.*t21.*1.055813973565752e-1;
t244 = qd2.*qd3.*t21.*1.055813973565752e-1;
t245 = qd1.*qd5.*t21.*1.055813973565752e-1;
t246 = qd2.*qd4.*t21.*1.055813973565752e-1;
t247 = qd1.*qd6.*t21.*1.055813973565752e-1;
t248 = qd2.*qd5.*t21.*1.055813973565752e-1;
t249 = qd3.*qd4.*t21.*1.055813973565752e-1;
t250 = qd2.*qd6.*t21.*1.055813973565752e-1;
t251 = qd3.*qd5.*t21.*1.055813973565752e-1;
t252 = qd3.*qd6.*t21.*1.055813973565752e-1;
t253 = qd4.*qd5.*t21.*1.055813973565752e-1;
t254 = qd4.*qd6.*t21.*1.055813973565752e-1;
t255 = qd5.*qd6.*t21.*1.055813973565752e-1;
t262 = t13.*t20.*1.702413339329298;
t263 = t14.*t20.*1.702413339329298;
t264 = t15.*t20.*1.702413339329298;
t265 = t16.*t20.*1.702413339329298;
t271 = t13.*t21.*5.279069867828761e-2;
t272 = t14.*t21.*5.279069867828761e-2;
t273 = t15.*t21.*5.279069867828761e-2;
t274 = t16.*t21.*5.279069867828761e-2;
t275 = t17.*t21.*5.279069867828761e-2;
t276 = t18.*t21.*5.279069867828761e-2;
t340 = -t331;
t341 = -t332;
t342 = -t333;
t343 = -t334;
t344 = -t335;
t350 = -t345;
t351 = -t346;
t352 = -t347;
t353 = -t348;
t359 = -t349;
t33 = cos(t29);
t36 = sin(t29);
t38 = sin(t35);
t39 = cos(t35);
t46 = -t40;
t47 = -t41;
t48 = -t42;
t52 = -t43;
t53 = -t44;
t54 = -t45;
t61 = -t50;
t62 = -t51;
t83 = qd1.*t31.*3.517602642454061;
t84 = qd2.*t31.*3.517602642454061;
t85 = qd3.*t31.*3.517602642454061;
t86 = qd4.*t31.*3.517602642454061;
t92 = qd1.*qd2.*t27.*3.517602642454061;
t93 = qd1.*qd3.*t27.*3.517602642454061;
t94 = qd1.*qd4.*t27.*3.517602642454061;
t95 = qd2.*qd3.*t27.*3.517602642454061;
t96 = qd2.*qd4.*t27.*3.517602642454061;
t97 = qd3.*qd4.*t27.*3.517602642454061;
t125 = t13.*t27.*1.75880132122703;
t126 = t14.*t27.*1.75880132122703;
t127 = t15.*t27.*1.75880132122703;
t128 = t16.*t27.*1.75880132122703;
t225 = -t215;
t226 = -t216;
t227 = -t217;
t228 = -t218;
t229 = -t219;
t230 = -t220;
t231 = -t221;
t232 = -t222;
t233 = -t223;
t234 = -t224;
t256 = qd1.*t37.*1.055813973565752e-1;
t257 = qd2.*t37.*1.055813973565752e-1;
t258 = qd3.*t37.*1.055813973565752e-1;
t259 = qd4.*t37.*1.055813973565752e-1;
t260 = qd5.*t37.*1.055813973565752e-1;
t261 = qd6.*t37.*1.055813973565752e-1;
t266 = -t236;
t267 = -t237;
t268 = -t238;
t269 = -t239;
t270 = -t240;
t277 = -t245;
t278 = -t247;
t279 = -t248;
t280 = -t250;
t281 = -t251;
t282 = -t252;
t283 = -t253;
t284 = -t254;
t285 = -t255;
t292 = -t264;
t293 = -t265;
t294 = qd1.*qd2.*t34.*1.055813973565752e-1;
t295 = qd1.*qd3.*t34.*1.055813973565752e-1;
t296 = qd1.*qd4.*t34.*1.055813973565752e-1;
t297 = qd2.*qd3.*t34.*1.055813973565752e-1;
t298 = qd1.*qd5.*t34.*1.055813973565752e-1;
t299 = qd2.*qd4.*t34.*1.055813973565752e-1;
t300 = qd1.*qd6.*t34.*1.055813973565752e-1;
t301 = qd2.*qd5.*t34.*1.055813973565752e-1;
t302 = qd3.*qd4.*t34.*1.055813973565752e-1;
t303 = qd2.*qd6.*t34.*1.055813973565752e-1;
t304 = qd3.*qd5.*t34.*1.055813973565752e-1;
t305 = qd3.*qd6.*t34.*1.055813973565752e-1;
t306 = qd4.*qd5.*t34.*1.055813973565752e-1;
t307 = qd4.*qd6.*t34.*1.055813973565752e-1;
t308 = qd5.*qd6.*t34.*1.055813973565752e-1;
t309 = -t275;
t310 = -t276;
t311 = t13.*t34.*5.279069867828761e-2;
t312 = t14.*t34.*5.279069867828761e-2;
t313 = t15.*t34.*5.279069867828761e-2;
t314 = t16.*t34.*5.279069867828761e-2;
t315 = t17.*t34.*5.279069867828761e-2;
t316 = t18.*t34.*5.279069867828761e-2;
t354 = qd1.*t32.*9.647505289768104e-1;
t355 = qd2.*t32.*9.647505289768104e-1;
t356 = qd3.*t32.*9.647505289768104e-1;
t357 = qd4.*t32.*9.647505289768104e-1;
t358 = qd5.*t32.*9.647505289768104e-1;
t360 = qd1.*qd2.*t28.*9.647505289768104e-1;
t361 = qd1.*qd3.*t28.*9.647505289768104e-1;
t362 = qd1.*qd4.*t28.*9.647505289768104e-1;
t363 = qd2.*qd3.*t28.*9.647505289768104e-1;
t364 = qd1.*qd5.*t28.*9.647505289768104e-1;
t365 = qd2.*qd4.*t28.*9.647505289768104e-1;
t366 = qd2.*qd5.*t28.*9.647505289768104e-1;
t367 = qd3.*qd4.*t28.*9.647505289768104e-1;
t368 = qd3.*qd5.*t28.*9.647505289768104e-1;
t369 = qd4.*qd5.*t28.*9.647505289768104e-1;
t384 = t13.*t28.*4.823752644884052e-1;
t385 = t14.*t28.*4.823752644884052e-1;
t386 = t15.*t28.*4.823752644884052e-1;
t387 = t16.*t28.*4.823752644884052e-1;
t388 = t17.*t28.*4.823752644884052e-1;
t394 = t55+t56+t57+t58+t59+t219+t220+t221+t222;
t396 = t75+t219+t220+t221+t222+t331+t332+t333+t334;
t88 = -t83;
t89 = -t84;
t90 = -t85;
t91 = -t86;
t98 = -t92;
t99 = -t93;
t100 = -t94;
t101 = -t95;
t102 = -t96;
t103 = -t97;
t104 = qd1.*t38.*1.09078504542792e-1;
t105 = qd2.*t38.*1.09078504542792e-1;
t106 = qd3.*t38.*1.09078504542792e-1;
t107 = qd4.*t38.*1.09078504542792e-1;
t108 = qd5.*t38.*1.09078504542792e-1;
t109 = qd6.*t38.*1.09078504542792e-1;
t110 = qd1.*qd2.*t39.*1.09078504542792e-1;
t111 = qd1.*qd3.*t39.*1.09078504542792e-1;
t112 = qd1.*qd4.*t39.*1.09078504542792e-1;
t113 = qd2.*qd3.*t39.*1.09078504542792e-1;
t114 = qd1.*qd5.*t39.*1.09078504542792e-1;
t115 = qd2.*qd4.*t39.*1.09078504542792e-1;
t116 = qd1.*qd6.*t39.*1.09078504542792e-1;
t117 = qd2.*qd5.*t39.*1.09078504542792e-1;
t118 = qd3.*qd4.*t39.*1.09078504542792e-1;
t119 = qd2.*qd6.*t39.*1.09078504542792e-1;
t120 = qd3.*qd5.*t39.*1.09078504542792e-1;
t121 = qd3.*qd6.*t39.*1.09078504542792e-1;
t122 = qd4.*qd5.*t39.*1.09078504542792e-1;
t123 = qd4.*qd6.*t39.*1.09078504542792e-1;
t124 = qd5.*qd6.*t39.*1.09078504542792e-1;
t150 = -t126;
t151 = -t127;
t152 = -t128;
t153 = t13.*t39.*5.4539252271396e-2;
t154 = t14.*t39.*5.4539252271396e-2;
t155 = t15.*t39.*5.4539252271396e-2;
t156 = t16.*t39.*5.4539252271396e-2;
t157 = t17.*t39.*5.4539252271396e-2;
t158 = t18.*t39.*5.4539252271396e-2;
t164 = qd1.*t36.*9.96705362804184e-1;
t165 = qd2.*t36.*9.96705362804184e-1;
t166 = qd3.*t36.*9.96705362804184e-1;
t167 = qd4.*t36.*9.96705362804184e-1;
t168 = qd5.*t36.*9.96705362804184e-1;
t169 = qd1.*qd2.*t33.*9.96705362804184e-1;
t170 = qd1.*qd3.*t33.*9.96705362804184e-1;
t171 = qd1.*qd4.*t33.*9.96705362804184e-1;
t172 = qd2.*qd3.*t33.*9.96705362804184e-1;
t173 = qd1.*qd5.*t33.*9.96705362804184e-1;
t174 = qd2.*qd4.*t33.*9.96705362804184e-1;
t175 = qd2.*qd5.*t33.*9.96705362804184e-1;
t176 = qd3.*qd4.*t33.*9.96705362804184e-1;
t177 = qd3.*qd5.*t33.*9.96705362804184e-1;
t178 = qd4.*qd5.*t33.*9.96705362804184e-1;
t194 = t13.*t33.*4.98352681402092e-1;
t195 = t14.*t33.*4.98352681402092e-1;
t196 = t15.*t33.*4.98352681402092e-1;
t197 = t16.*t33.*4.98352681402092e-1;
t198 = t17.*t33.*4.98352681402092e-1;
t286 = -t256;
t287 = -t257;
t288 = -t258;
t289 = -t259;
t290 = -t260;
t291 = -t261;
t317 = -t295;
t318 = -t296;
t319 = -t297;
t320 = -t298;
t321 = -t299;
t322 = -t300;
t323 = -t301;
t324 = -t302;
t325 = -t303;
t326 = -t304;
t327 = -t305;
t328 = -t306;
t329 = -t307;
t330 = -t308;
t336 = -t313;
t337 = -t314;
t338 = -t315;
t339 = -t316;
t370 = -t354;
t371 = -t355;
t372 = -t356;
t373 = -t357;
t374 = -t358;
t375 = -t361;
t376 = -t362;
t377 = -t363;
t378 = -t364;
t379 = -t365;
t380 = -t366;
t381 = -t367;
t382 = -t368;
t383 = -t369;
t389 = -t386;
t390 = -t387;
t391 = -t388;
t392 = t75+t233+t234+t344;
t397 = t70+t71+t72+t73+t74+t75+t229+t230+t231+t232+t233+t234;
t129 = -t104;
t130 = -t105;
t131 = -t106;
t132 = -t107;
t133 = -t108;
t134 = -t109;
t135 = -t110;
t136 = -t111;
t137 = -t112;
t138 = -t113;
t139 = -t114;
t140 = -t115;
t141 = -t116;
t142 = -t117;
t143 = -t118;
t144 = -t119;
t145 = -t120;
t146 = -t121;
t147 = -t122;
t148 = -t123;
t149 = -t124;
t159 = -t154;
t160 = -t155;
t161 = -t156;
t162 = -t157;
t163 = -t158;
t179 = -t164;
t180 = -t165;
t181 = -t166;
t182 = -t167;
t183 = -t168;
t184 = -t169;
t185 = -t170;
t186 = -t171;
t187 = -t172;
t188 = -t173;
t189 = -t174;
t190 = -t175;
t191 = -t176;
t192 = -t177;
t193 = -t178;
t199 = -t195;
t200 = -t196;
t201 = -t197;
t202 = -t198;
t393 = t153+t294+t311+t312;
t399 = t229+t230+t231+t232+t340+t341+t342+t343+t392;
t395 = t194+t360+t384+t385+t393;
t400 = t76+t77+t78+t79+t80+t87+t277+t278+t279+t280+t281+t282+t283+t284+t285+t309+t310+t393;
t398 = t125+t235+t262+t263+t395;
t401 = t277+t278+t279+t280+t281+t282+t283+t284+t285+t309+t310+t350+t351+t352+t353+t359+t395;
dV_dx = reshape([0.0,t52+t53+t54+t61+t62+t98+t99+t100+t101+t102+t103+t135+t136+t137+t138+t139+t140+t141+t142+t143+t144+t145+t146+t147+t148+t149+t150+t151+t152+t159+t160+t161+t162+t163+t184+t185+t186+t187+t188+t189+t190+t191+t192+t193+t199+t200+t201+t202-t2.*t14.*1.163815784709457e+1-qd1.*qd2.*t2.*2.327631569418914e+1,t52+t53+t54+t61+t62+t98+t99+t100+t101+t102+t103+t135+t136+t137+t138+t139+t140+t141+t142+t143+t144+t145+t146+t147+t148+t149+t150+t151+t152+t159+t160+t161+t162+t163+t184+t185+t186+t187+t188+t189+t190+t191+t192+t193+t199+t200+t201+t202+t212+t213+t214+t266+t267+t268+t269+t270+t292+t293+t317+t318+t319+t320+t321+t322+t323+t324+t325+t326+t327+t328+t329+t330+t336+t337+t338+t339+t375+t376+t377+t378+t379+t380+t381+t382+t383+t389+t390+t391,t98+t99+t100+t101+t102+t103+t135+t136+t137+t138+t139+t140+t141+t142+t143+t144+t145+t146+t147+t148+t149+t150+t151+t152+t159+t160+t161+t162+t163+t184+t185+t186+t187+t188+t189+t190+t191+t192+t193+t199+t200+t201+t202+t266+t267+t268+t269+t270+t292+t293+t317+t318+t319+t320+t321+t322+t323+t324+t325+t326+t327+t328+t329+t330+t336+t337+t338+t339+t375+t376+t377+t378+t379+t380+t381+t382+t383+t389+t390+t391,t135+t136+t137+t138+t139+t140+t141+t142+t143+t144+t145+t146+t147+t148+t149+t159+t160+t161+t162+t163+t184+t185+t186+t187+t188+t189+t190+t191+t192+t193+t199+t200+t201+t202+t277+t278+t279+t280+t281+t282+t283+t284+t285+t309+t310+t317+t318+t319+t320+t321+t322+t323+t324+t325+t326+t327+t328+t329+t330+t336+t337+t338+t339+t350+t351+t352+t353+t359+t375+t376+t377+t378+t379+t380+t381+t382+t383+t389+t390+t391,t76+t77+t78+t79+t80+t87+t135+t136+t137+t138+t139+t140+t141+t142+t143+t144+t145+t146+t147+t148+t149+t159+t160+t161+t162+t163+t277+t278+t279+t280+t281+t282+t283+t284+t285+t309+t310+t317+t318+t319+t320+t321+t322+t323+t324+t325+t326+t327+t328+t329+t330+t336+t337+t338+t339,t47+t48+t81+t89+t90+t91+t130+t131+t132+t133+t134+t180+t181+t182+t183+t210+t227+t228+t288+t289+t290+t291+t372+t373+t374+t392,t46+t47+t48-t68+t81+t88+t89+t90+t91+t129+t130+t131+t132+t133+t134+t179+t180+t181+t182+t183+t210+t227+t228+t288+t289+t290+t291+t372+t373+t374+t392,t46+t47+t48+t88+t89+t90+t91+t129+t130+t131+t132+t133+t134+t179+t180+t181+t182+t183+t208+t209+t210+t225+t226+t227+t228+t286+t287+t288+t289+t290+t291+t370+t371+t372+t373+t374+t392,t88+t89+t90+t91+t129+t130+t131+t132+t133+t134+t179+t180+t181+t182+t183+t225+t226+t227+t228+t286+t287+t288+t289+t290+t291+t370+t371+t372+t373+t374+t392,t129+t130+t131+t132+t133+t134+t179+t180+t181+t182+t183+t286+t287+t288+t289+t290+t291+t370+t371+t372+t373+t374+t399,t129+t130+t131+t132+t133+t134+t286+t287+t288+t289+t290+t291+t397,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t49+t125+t153+t194+t2.*t13.*1.163815784709457e+1,t49+t125+t153+t194+t212+t213+t214+t266+t267+t268+t269+t270+t292+t293+t317+t318+t319+t320+t321+t322+t323+t324+t325+t326+t327+t328+t329+t330+t336+t337+t338+t339+t375+t376+t377+t378+t379+t380+t381+t382+t383+t389+t390+t391,t125+t153+t194+t266+t267+t268+t269+t270+t292+t293+t317+t318+t319+t320+t321+t322+t323+t324+t325+t326+t327+t328+t329+t330+t336+t337+t338+t339+t375+t376+t377+t378+t379+t380+t381+t382+t383+t389+t390+t391,t153+t194+t277+t278+t279+t280+t281+t282+t283+t284+t285+t309+t310+t317+t318+t319+t320+t321+t322+t323+t324+t325+t326+t327+t328+t329+t330+t336+t337+t338+t339+t350+t351+t352+t353+t359+t375+t376+t377+t378+t379+t380+t381+t382+t383+t389+t390+t391,t76+t77+t78+t79+t80+t87+t153+t277+t278+t279+t280+t281+t282+t283+t284+t285+t309+t310+t317+t318+t319+t320+t321+t322+t323+t324+t325+t326+t327+t328+t329+t330+t336+t337+t338+t339,t40+t68+t83+t104+t164+t210+t227+t228+t288+t289+t290+t291+t372+t373+t374+t392,t210+t227+t228+t288+t289+t290+t291+t372+t373+t374+t392,t208+t209+t210+t225+t226+t227+t228+t286+t287+t288+t289+t290+t291+t370+t371+t372+t373+t374+t392,t225+t226+t227+t228+t286+t287+t288+t289+t290+t291+t370+t371+t372+t373+t374+t392,t286+t287+t288+t289+t290+t291+t370+t371+t372+t373+t374+t399,t286+t287+t288+t289+t290+t291+t397,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t49+t125+t153+t194,t49+t398+t3.*t13.*7.571420887311746+t3.*t14.*7.571420887311746+qd1.*qd2.*t3.*1.514284177462349e+1,t398,t401,t400,t40+t83+t104+t164+t203+t204+t215+t216+t256+t257+t354+t355+t392,t203+t204+t215+t216+t256+t257+t354+t355+t392,t392,t392,t399,t397,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t125+t153+t194,t398,t398,t401,t400,t83+t104+t164+t215+t216+t256+t257+t354+t355+t392,t215+t216+t256+t257+t354+t355+t392,t392,t392,t399,t397,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t153+t194,t395,t395,t241+t242+t243+t244+t246+t249+t271+t272+t273+t274+t395+t4.*t13.*4.823752644884052e-1+t4.*t14.*4.823752644884052e-1+t4.*t15.*4.823752644884052e-1+t4.*t16.*4.823752644884052e-1+qd1.*qd2.*t4.*9.647505289768104e-1+qd1.*qd3.*t4.*9.647505289768104e-1+qd1.*qd4.*t4.*9.647505289768104e-1+qd2.*qd3.*t4.*9.647505289768104e-1+qd2.*qd4.*t4.*9.647505289768104e-1+qd3.*qd4.*t4.*9.647505289768104e-1,t76+t77+t78+t79+t80+t87+t241+t242+t243+t244+t246+t249+t271+t272+t273+t274+t393,t104+t164+t256+t257+t354+t355+t396,t256+t257+t354+t355+t396,t396,t396,t75,t70+t71+t72+t73+t74+t75,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t153,t393,t393,t241+t242+t243+t244+t246+t249+t271+t272+t273+t274+t393,t241+t242+t243+t244+t246+t249+t271+t272+t273+t274+t393+t5.*t13.*5.4539252271396e-2+t5.*t14.*5.4539252271396e-2+t5.*t15.*5.4539252271396e-2+t5.*t16.*5.4539252271396e-2+t5.*t17.*5.4539252271396e-2+qd1.*qd2.*t5.*1.09078504542792e-1+qd1.*qd3.*t5.*1.09078504542792e-1+qd1.*qd4.*t5.*1.09078504542792e-1+qd2.*qd3.*t5.*1.09078504542792e-1+qd1.*qd5.*t5.*1.09078504542792e-1+qd2.*qd4.*t5.*1.09078504542792e-1+qd2.*qd5.*t5.*1.09078504542792e-1+qd3.*qd4.*t5.*1.09078504542792e-1+qd3.*qd5.*t5.*1.09078504542792e-1+qd4.*qd5.*t5.*1.09078504542792e-1,t104+t256+t257+t394,t256+t257+t394,t394,t394,t55+t56+t57+t58+t59,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[18,6]);
