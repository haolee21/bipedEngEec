function dMext_heel_dx5 = dMext_heel_dx5(in1)
%DMEXT_HEEL_DX5
%    DMEXT_HEEL_DX5 = DMEXT_HEEL_DX5(IN1)

%    This function was generated by the Symbolic Math Toolbox version 8.4.
%    27-May-2020 02:17:26

q1 = in1(:,1);
q2 = in1(:,2);
q3 = in1(:,3);
q4 = in1(:,4);
q5 = in1(:,5);
t2 = cos(q1);
t3 = sin(q1);
t4 = q1+q2;
t5 = t2.^2;
t6 = t3.^2;
t7 = cos(t4);
t8 = sin(t4);
t9 = q3+q4+t4;
t10 = t7.^2;
t11 = t8.^2;
t12 = cos(t9);
t13 = q5+t9;
t14 = sin(t9);
t53 = t2.*t7.*t8.*7.744424728713496e+57;
t54 = t3.*t7.*t8.*7.744424728713496e+57;
t213 = t2.*t3.*t7.*t8.*1.97286209645842e+77;
t15 = sin(t13);
t16 = cos(t13);
t17 = t12.^2;
t18 = t14.^2;
t45 = t3.*t10.*7.744424728713496e+57;
t46 = t2.*t11.*7.744424728713496e+57;
t49 = t6.*t12.*8.000938509076079e+57;
t50 = t5.*t14.*8.000938509076079e+57;
t55 = t2.*t3.*t12.*8.000938509076079e+57;
t56 = t2.*t3.*t14.*8.000938509076079e+57;
t75 = t11.*t12.*1.499226979701784e+58;
t76 = t10.*t14.*1.499226979701784e+58;
t79 = t3.*t7.*t12.*7.744424728713496e+57;
t80 = t3.*t7.*t12.*1.548884945742699e+58;
t81 = t2.*t8.*t12.*7.744424728713496e+57;
t82 = t2.*t7.*t14.*7.744424728713496e+57;
t83 = t2.*t7.*t14.*1.548884945742699e+58;
t84 = t3.*t8.*t12.*7.744424728713496e+57;
t85 = t3.*t8.*t12.*1.548884945742699e+58;
t86 = t3.*t7.*t14.*7.744424728713496e+57;
t87 = t2.*t8.*t14.*7.744424728713496e+57;
t88 = t2.*t8.*t14.*1.548884945742699e+58;
t129 = t2.*t12.*t14.*7.744424728713496e+57;
t130 = t3.*t12.*t14.*7.744424728713496e+57;
t135 = t7.*t8.*t12.*1.499226979701784e+58;
t136 = t7.*t8.*t14.*1.499226979701784e+58;
t170 = t7.*t12.*t14.*1.499226979701784e+58;
t172 = t8.*t12.*t14.*1.499226979701784e+58;
t211 = t6.*t10.*9.864310482292101e+76;
t212 = t5.*t11.*9.864310482292101e+76;
t216 = -t213;
t221 = t2.*t3.*t8.*t12.*1.97286209645842e+77;
t222 = t2.*t3.*t7.*t14.*1.97286209645842e+77;
t225 = t6.*t7.*t12.*1.97286209645842e+77;
t227 = t5.*t8.*t14.*1.97286209645842e+77;
t248 = t2.*t3.*t12.*t14.*5.918586289375261e+77;
t272 = t3.*t7.*t12.*t14.*3.819222454671376e+77;
t273 = t2.*t8.*t12.*t14.*3.819222454671376e+77;
t321 = t7.*t8.*t12.*t14.*7.393552841047993e+77;
t19 = t16.^2;
t20 = t15.^2;
t21 = t2.*t3.*t16.*3.244779102448822e+38;
t22 = t2.*t3.*t15.*3.244779102448822e+38;
t25 = t3.*t7.*t16.*3.140749987231091e+38;
t26 = t2.*t8.*t16.*3.140749987231091e+38;
t27 = t3.*t7.*t15.*3.140749987231091e+38;
t28 = t2.*t8.*t15.*3.140749987231091e+38;
t31 = t7.*t8.*t16.*6.080112186895894e+38;
t32 = t7.*t8.*t15.*6.080112186895894e+38;
t33 = t3.*t12.*t16.*3.140749987231091e+38;
t34 = t2.*t14.*t16.*3.140749987231091e+38;
t35 = t3.*t12.*t15.*3.140749987231091e+38;
t36 = t2.*t14.*t15.*3.140749987231091e+38;
t39 = t8.*t12.*t16.*6.080112186895894e+38;
t40 = t7.*t14.*t16.*6.080112186895894e+38;
t41 = t8.*t12.*t15.*6.080112186895894e+38;
t42 = t7.*t14.*t15.*6.080112186895894e+38;
t43 = t12.*t14.*t16.*1.216022437379179e+39;
t44 = t12.*t14.*t15.*1.216022437379179e+39;
t47 = -t45;
t48 = -t46;
t51 = t3.*t17.*7.744424728713496e+57;
t52 = t2.*t18.*7.744424728713496e+57;
t57 = t2.*t3.*t16.*9.663334543508915e+57;
t58 = t2.*t3.*t15.*9.663334543508915e+57;
t62 = t5.*t16.*9.663334543508915e+57;
t63 = t5.*t15.*9.663334543508915e+57;
t66 = t6.*t16.*9.663334543508915e+57;
t69 = -t55;
t70 = t6.*t15.*9.663334543508915e+57;
t72 = -t56;
t77 = t8.*t17.*1.499226979701784e+58;
t78 = t7.*t18.*1.499226979701784e+58;
t93 = t3.*t8.*t15.*9.353523579226171e+57;
t94 = t3.*t8.*t15.*1.870704715845234e+58;
t99 = t10.*t16.*1.810729059482018e+58;
t100 = t10.*t15.*1.810729059482018e+58;
t103 = t11.*t16.*1.810729059482018e+58;
t106 = -t80;
t107 = -t79;
t108 = -t81;
t109 = t11.*t15.*1.810729059482018e+58;
t112 = -t86;
t113 = -t88;
t114 = -t87;
t115 = t2.*t7.*t16.*9.353523579226171e+57;
t116 = t2.*t7.*t16.*1.870704715845234e+58;
t117 = t2.*t7.*t15.*9.353523579226171e+57;
t118 = t2.*t7.*t15.*1.870704715845234e+58;
t119 = t3.*t7.*t16.*9.353523579226171e+57;
t120 = t3.*t7.*t16.*1.870704715845234e+58;
t121 = t2.*t8.*t16.*9.353523579226171e+57;
t122 = t2.*t8.*t16.*1.870704715845234e+58;
t123 = t3.*t7.*t15.*9.353523579226171e+57;
t124 = t3.*t7.*t15.*1.870704715845234e+58;
t125 = t2.*t8.*t15.*9.353523579226171e+57;
t126 = t2.*t8.*t15.*1.870704715845234e+58;
t127 = t3.*t8.*t16.*9.353523579226171e+57;
t128 = t3.*t8.*t16.*1.870704715845234e+58;
t147 = -t135;
t148 = -t136;
t149 = t7.*t8.*t16.*1.810729059482018e+58;
t150 = t2.*t12.*t16.*9.353523579226171e+57;
t151 = t7.*t8.*t15.*1.810729059482018e+58;
t152 = t2.*t12.*t15.*9.353523579226171e+57;
t153 = t3.*t12.*t16.*9.353523579226171e+57;
t154 = t3.*t12.*t16.*1.870704715845234e+58;
t155 = t2.*t14.*t16.*9.353523579226171e+57;
t156 = t2.*t14.*t16.*1.870704715845234e+58;
t157 = t3.*t12.*t15.*9.353523579226171e+57;
t158 = t3.*t12.*t15.*1.870704715845234e+58;
t159 = t2.*t14.*t15.*9.353523579226171e+57;
t160 = t2.*t14.*t15.*1.870704715845234e+58;
t161 = t3.*t14.*t16.*9.353523579226171e+57;
t162 = t3.*t14.*t15.*9.353523579226171e+57;
t163 = t2.*t15.*t16.*1.129695315169956e+58;
t164 = t2.*t15.*t16.*2.259390630339913e+58;
t166 = t3.*t15.*t16.*1.129695315169956e+58;
t167 = t3.*t15.*t16.*2.259390630339913e+58;
t179 = t7.*t12.*t16.*1.810729059482018e+58;
t180 = t7.*t12.*t15.*1.810729059482018e+58;
t181 = t8.*t12.*t16.*1.810729059482018e+58;
t182 = t8.*t12.*t16.*3.621458118964036e+58;
t183 = t7.*t14.*t16.*1.810729059482018e+58;
t184 = t7.*t14.*t16.*3.621458118964036e+58;
t185 = t8.*t12.*t15.*1.810729059482018e+58;
t186 = t8.*t12.*t15.*3.621458118964036e+58;
t187 = t7.*t14.*t15.*1.810729059482018e+58;
t188 = t7.*t14.*t15.*3.621458118964036e+58;
t189 = t8.*t14.*t16.*1.810729059482018e+58;
t190 = t8.*t14.*t15.*1.810729059482018e+58;
t191 = t7.*t15.*t16.*2.18695352421207e+58;
t192 = t7.*t15.*t16.*3.280430286318105e+58;
t193 = t7.*t15.*t16.*4.37390704842414e+58;
t194 = t8.*t15.*t16.*2.18695352421207e+58;
t195 = t8.*t15.*t16.*3.280430286318105e+58;
t196 = t8.*t15.*t16.*4.37390704842414e+58;
t201 = t7.*t15.*t16.*6.560860572636209e+58;
t202 = t8.*t15.*t16.*6.560860572636209e+58;
t203 = t12.*t15.*t16.*1.093476762106035e+58;
t204 = t12.*t15.*t16.*2.18695352421207e+58;
t205 = t14.*t15.*t16.*1.093476762106035e+58;
t206 = t14.*t15.*t16.*2.18695352421207e+58;
t214 = t6.*t17.*2.95929314468763e+77;
t215 = t5.*t18.*2.95929314468763e+77;
t219 = t11.*t17.*3.696776420523996e+77;
t220 = t10.*t18.*3.696776420523996e+77;
t226 = t2.*t7.*t18.*3.819222454671376e+77;
t228 = t3.*t8.*t17.*3.819222454671376e+77;
t230 = t6.*t7.*t16.*2.382773773934124e+77;
t232 = t5.*t8.*t16.*2.382773773934124e+77;
t234 = t6.*t7.*t15.*2.382773773934124e+77;
t235 = t5.*t8.*t15.*2.382773773934124e+77;
t239 = -t221;
t240 = -t222;
t242 = t2.*t3.*t7.*t16.*2.382773773934124e+77;
t243 = t2.*t3.*t7.*t15.*2.382773773934124e+77;
t244 = t2.*t3.*t8.*t16.*2.382773773934124e+77;
t247 = t2.*t3.*t8.*t15.*2.382773773934124e+77;
t256 = t6.*t12.*t16.*7.148321321802371e+77;
t258 = t5.*t14.*t16.*7.148321321802371e+77;
t260 = t6.*t12.*t15.*7.148321321802371e+77;
t261 = t5.*t14.*t15.*7.148321321802371e+77;
t265 = -t248;
t266 = t2.*t3.*t12.*t16.*7.148321321802371e+77;
t267 = t2.*t3.*t12.*t15.*7.148321321802371e+77;
t268 = t2.*t3.*t14.*t16.*7.148321321802371e+77;
t271 = t2.*t3.*t14.*t15.*7.148321321802371e+77;
t282 = t11.*t12.*t15.*8.929749239681125e+77;
t283 = t10.*t14.*t15.*8.929749239681125e+77;
t287 = -t272;
t288 = -t273;
t289 = t2.*t3.*t15.*t16.*1.151141961303919e+78;
t291 = t3.*t7.*t12.*t16.*4.612761894583371e+77;
t292 = t2.*t8.*t12.*t16.*4.612761894583371e+77;
t293 = t2.*t7.*t14.*t16.*9.225523789166741e+77;
t294 = t3.*t7.*t12.*t15.*4.612761894583371e+77;
t295 = t2.*t8.*t12.*t15.*4.612761894583371e+77;
t296 = t2.*t7.*t14.*t15.*9.225523789166741e+77;
t297 = t3.*t8.*t12.*t16.*9.225523789166741e+77;
t298 = t3.*t7.*t14.*t16.*4.612761894583371e+77;
t299 = t2.*t8.*t14.*t16.*4.612761894583371e+77;
t300 = t3.*t8.*t12.*t15.*9.225523789166741e+77;
t301 = t3.*t7.*t14.*t15.*4.612761894583371e+77;
t302 = t2.*t8.*t14.*t15.*4.612761894583371e+77;
t303 = t5.*t15.*t16.*1.151141961303919e+78;
t304 = t6.*t15.*t16.*1.151141961303919e+78;
t306 = t11.*t12.*t16.*8.929749239681125e+77;
t308 = t10.*t14.*t16.*8.929749239681125e+77;
t322 = t3.*t7.*t15.*t16.*8.356768641518365e+77;
t323 = t2.*t8.*t15.*t16.*8.356768641518365e+77;
t325 = t7.*t8.*t12.*t16.*8.929749239681125e+77;
t326 = t7.*t8.*t12.*t15.*8.929749239681125e+77;
t327 = t7.*t8.*t14.*t16.*8.929749239681125e+77;
t328 = t7.*t8.*t14.*t15.*8.929749239681125e+77;
t329 = t10.*t15.*t16.*1.617769356584788e+78;
t330 = t11.*t15.*t16.*1.617769356584788e+78;
t331 = -t321;
t332 = t2.*t7.*t15.*t16.*1.671353728303673e+78;
t335 = t3.*t8.*t15.*t16.*1.671353728303673e+78;
t340 = t2.*t12.*t15.*t16.*5.571179094345577e+77;
t341 = t3.*t12.*t15.*t16.*2.785589547172788e+77;
t342 = t2.*t14.*t15.*t16.*2.785589547172788e+77;
t345 = t3.*t14.*t15.*t16.*5.571179094345577e+77;
t347 = t15.*t16.*t17.*1.078512904389859e+78;
t348 = t15.*t16.*t18.*1.078512904389859e+78;
t349 = t7.*t8.*t15.*t16.*1.617769356584788e+78;
t355 = t8.*t12.*t15.*t16.*5.392564521949295e+77;
t356 = t7.*t14.*t15.*t16.*5.392564521949295e+77;
t358 = t7.*t12.*t15.*t16.*1.078512904389859e+78;
t361 = t8.*t14.*t15.*t16.*1.078512904389859e+78;
t363 = t12.*t14.*t15.*t16.*1.078512904389859e+78;
t23 = t3.*t19.*3.793319722393432e+38;
t24 = t2.*t20.*3.793319722393432e+38;
t29 = t8.*t19.*7.34340828358964e+38;
t30 = t7.*t20.*7.34340828358964e+38;
t37 = t14.*t19.*1.468681656717928e+39;
t38 = t12.*t20.*1.468681656717928e+39;
t59 = -t51;
t60 = -t52;
t61 = t2.*t19.*1.129695315169956e+58;
t64 = t3.*t19.*1.129695315169956e+58;
t65 = t3.*t19.*4.518781260679825e+58;
t67 = t2.*t20.*1.129695315169956e+58;
t68 = t2.*t20.*4.518781260679825e+58;
t71 = t3.*t20.*1.129695315169956e+58;
t73 = -t57;
t74 = -t58;
t95 = -t77;
t96 = -t78;
t97 = t7.*t19.*2.18695352421207e+58;
t98 = t7.*t19.*3.280430286318105e+58;
t101 = t8.*t19.*2.18695352421207e+58;
t102 = t8.*t19.*3.280430286318105e+58;
t104 = t7.*t20.*2.18695352421207e+58;
t105 = t7.*t20.*3.280430286318105e+58;
t110 = t8.*t20.*2.18695352421207e+58;
t111 = t8.*t20.*3.280430286318105e+58;
t137 = -t120;
t138 = -t119;
t139 = -t121;
t140 = -t123;
t141 = -t126;
t142 = -t125;
t143 = t12.*t19.*1.093476762106035e+58;
t144 = t14.*t19.*1.093476762106035e+58;
t145 = t12.*t20.*1.093476762106035e+58;
t146 = t14.*t20.*1.093476762106035e+58;
t173 = -t149;
t174 = -t151;
t175 = -t154;
t176 = -t153;
t177 = -t160;
t178 = -t159;
t197 = -t182;
t198 = -t181;
t199 = -t188;
t200 = -t187;
t207 = -t204;
t208 = -t203;
t209 = -t206;
t210 = -t205;
t217 = t6.*t19.*5.755709806519594e+77;
t218 = t5.*t20.*5.755709806519594e+77;
t223 = t2.*t3.*t19.*1.151141961303919e+78;
t224 = t2.*t3.*t20.*1.151141961303919e+78;
t229 = t3.*t7.*t19.*8.356768641518365e+77;
t231 = t2.*t8.*t19.*8.356768641518365e+77;
t233 = t2.*t7.*t20.*8.356768641518365e+77;
t236 = t3.*t8.*t19.*8.356768641518365e+77;
t237 = t3.*t7.*t20.*8.356768641518365e+77;
t238 = t2.*t8.*t20.*8.356768641518365e+77;
t245 = t11.*t19.*8.088846782923942e+77;
t246 = t10.*t20.*8.088846782923942e+77;
t249 = -t232;
t252 = -t243;
t253 = -t244;
t254 = -t247;
t255 = t3.*t12.*t19.*2.785589547172788e+77;
t257 = t2.*t14.*t19.*2.785589547172788e+77;
t259 = t2.*t12.*t20.*2.785589547172788e+77;
t262 = t3.*t14.*t19.*2.785589547172788e+77;
t263 = t3.*t12.*t20.*2.785589547172788e+77;
t264 = t2.*t14.*t20.*2.785589547172788e+77;
t269 = t17.*t20.*5.392564521949295e+77;
t270 = t18.*t19.*5.392564521949295e+77;
t274 = t7.*t8.*t19.*1.617769356584788e+78;
t275 = t7.*t8.*t20.*1.617769356584788e+78;
t276 = -t258;
t279 = -t267;
t280 = -t268;
t281 = -t271;
t284 = t8.*t14.*t19.*5.392564521949295e+77;
t285 = t8.*t12.*t20.*5.392564521949295e+77;
t286 = t7.*t14.*t20.*5.392564521949295e+77;
t305 = t8.*t12.*t19.*5.392564521949295e+77;
t307 = t7.*t14.*t19.*5.392564521949295e+77;
t309 = t7.*t12.*t20.*5.392564521949295e+77;
t312 = -t289;
t313 = -t293;
t314 = -t294;
t315 = -t295;
t316 = -t298;
t317 = -t299;
t318 = -t301;
t319 = -t302;
t320 = -t303;
t324 = -t308;
t333 = -t322;
t334 = -t323;
t336 = -t326;
t337 = -t327;
t338 = -t328;
t339 = -t329;
t343 = t12.*t14.*t19.*1.078512904389859e+78;
t344 = -t332;
t346 = t12.*t14.*t20.*1.078512904389859e+78;
t350 = -t340;
t351 = -t341;
t352 = -t342;
t354 = -t347;
t357 = -t349;
t359 = -t355;
t360 = -t356;
t362 = -t358;
t364 = -t363;
t89 = -t65;
t90 = -t64;
t91 = -t68;
t92 = -t67;
t131 = -t102;
t132 = -t101;
t133 = -t105;
t134 = -t104;
t165 = -t143;
t168 = -t144;
t169 = -t145;
t171 = -t146;
t241 = -t224;
t250 = -t237;
t251 = -t238;
t277 = -t263;
t278 = -t264;
t290 = -t275;
t310 = -t285;
t311 = -t286;
t353 = -t346;
t370 = t211+t212+t214+t215+t216+t217+t218+t219+t220+t225+t226+t227+t228+t230+t233+t235+t236+t239+t240+t245+t246+t252+t253+t256+t259+t261+t262+t265+t269+t270+t279+t280+t283+t284+t287+t288+t296+t297+t306+t309+t312+t314+t315+t316+t317+t331+t333+t334+t336+t337+t351+t352+t357+t359+t360+t364;
t365 = t57+t70+t71+t90+t94+t109+t110+t119+t121+t132+t144+t149+t156+t162+t164+t171+t176+t184+t190+t193+t198+t207;
t366 = t58+t61+t62+t92+t97+t99+t116+t123+t125+t134+t145+t150+t151+t158+t165+t167+t178+t179+t186+t196+t200+t209;
t367 = t49+t60+t66+t72+t74+t75+t85+t92+t96+t103+t112+t114+t128+t130+t134+t140+t142+t145+t148+t157+t161+t166+t172+t174+t177+t185+t189+t194+t199+t210;
t368 = t50+t59+t63+t69+t73+t76+t83+t90+t95+t100+t107+t108+t118+t129+t132+t138+t139+t144+t147+t152+t155+t163+t170+t173+t175+t180+t183+t191+t197+t208;
t369 = t223+t229+t231+t234+t241+t242+t249+t250+t251+t254+t255+t257+t260+t266+t274+t276+t277+t278+t281+t282+t290+t291+t292+t300+t304+t305+t307+t310+t311+t313+t318+t319+t320+t324+t325+t330+t335+t338+t339+t343+t344+t345+t348+t350+t353+t354+t361+t362;
t371 = 1.0./t370;
t372 = t371.^2;
t373 = t366.*t371.*2.81474976710656e+19;
t374 = t365.*t371.*2.81474976710656e+19;
t375 = -t373;
t376 = -t374;
t377 = t368.*t369.*t372.*2.81474976710656e+19;
t378 = t367.*t369.*t372.*2.81474976710656e+19;
t379 = -t377;
t380 = t376+t378;
t381 = t375+t379;
dMext_heel_dx5 = reshape([t371.*(t91+t98+t115+t124+t133+t142+t143+t169+t179+t186+t200+t202+t206+t2.*t19.*4.518781260679825e+58+t2.*t12.*t16.*2.806057073767851e+58+t3.*t12.*t15.*5.612114147535703e+58-t2.*t14.*t15.*2.806057073767851e+58+t3.*t15.*t16.*9.03756252135965e+58).*2.81474976710656e+19+t369.*t372.*(t47+t53+t81+t82+t89+t95+t106+t117+t121+t131+t137+t168+t170+t180+t183+t192+t197+t203-t3.*t17.*2.323327418614049e+58+t2.*t12.*t14.*2.323327418614049e+58+t2.*t12.*t15.*2.806057073767851e+58-t3.*t12.*t16.*5.612114147535703e+58+t2.*t14.*t16.*2.806057073767851e+58+t2.*t15.*t16.*4.518781260679825e+58).*2.81474976710656e+19,t371.*(t89+t93+t111+t122+t131+t138+t146+t168+t184+t190+t198+t201+t204+t3.*t20.*4.518781260679825e+58-t3.*t12.*t16.*2.806057073767851e+58+t2.*t14.*t16.*5.612114147535703e+58+t3.*t14.*t15.*2.806057073767851e+58+t2.*t15.*t16.*9.03756252135965e+58).*2.81474976710656e+19-t369.*t372.*(t48+t54+t84+t86+t91+t96+t113+t123+t127+t133+t141+t169+t172+t185+t189+t195+t199+t205-t2.*t18.*2.323327418614049e+58+t3.*t12.*t14.*2.323327418614049e+58+t3.*t12.*t15.*2.806057073767851e+58-t2.*t14.*t15.*5.612114147535703e+58+t3.*t14.*t16.*2.806057073767851e+58+t3.*t15.*t16.*4.518781260679825e+58).*2.81474976710656e+19,t371.*(t58+t61+t62+t92-t98+t105+t115+t124+t142+t145+t150+t158+t165+t167+t178-t179-t186+t187-t202+t209).*-2.81474976710656e+19-t369.*t372.*(t47+t50+t53+t59+t63+t69+t73+t77+t81+t82+t90+t102+t106+t117+t121+t129+t137+t144+t152+t155+t163-t170+t175-t180+t182-t183-t192+t208+t5.*t8.*8.000938509076079e+57-t2.*t3.*t7.*8.000938509076079e+57).*2.81474976710656e+19,t371.*(t57+t70+t71+t90+t93+t102-t111+t122+t138+t144+t156+t162+t164+t171+t176+t181-t184-t190-t201+t207).*-2.81474976710656e+19+t369.*t372.*(t48+t49+t54+t60+t66+t72+t74+t78+t84+t86+t92+t105+t113+t123+t127+t130+t141+t145+t157+t161+t166-t172+t177-t185+t188-t189-t195+t210+t6.*t7.*8.000938509076079e+57-t2.*t3.*t8.*8.000938509076079e+57).*2.81474976710656e+19,t381,t380,t381,t380,t371.*(t22-t24+t27+t28-t30+t32+t35+t36-t38+t41+t42+t44+t2.*t19.*3.793319722393432e+38+t5.*t16.*3.244779102448822e+38+t7.*t19.*7.34340828358964e+38+t10.*t16.*6.080112186895894e+38+t12.*t19.*1.468681656717928e+39+t16.*t17.*1.216022437379179e+39+t2.*t7.*t16.*6.281499974462183e+38+t2.*t12.*t16.*6.281499974462183e+38+t3.*t15.*t16.*7.586639444786864e+38+t7.*t12.*t16.*1.216022437379179e+39+t8.*t15.*t16.*1.468681656717928e+39+t14.*t15.*t16.*2.937363313435856e+39).*-8.382656506659224e+38+t369.*t372.*(t21+t23+t25+t26+t29+t31+t33+t34+t37+t39+t40+t43-t5.*t15.*3.244779102448822e+38-t10.*t15.*6.080112186895894e+38-t15.*t17.*1.216022437379179e+39-t2.*t7.*t15.*6.281499974462183e+38-t2.*t12.*t15.*6.281499974462183e+38-t2.*t15.*t16.*3.793319722393432e+38-t7.*t12.*t15.*1.216022437379179e+39-t7.*t15.*t16.*7.34340828358964e+38-t12.*t15.*t16.*1.468681656717928e+39).*8.382656506659224e+38,t371.*(t21-t23+t25+t26-t29+t31+t33+t34-t37+t39+t40+t43+t6.*t15.*3.244779102448822e+38+t3.*t20.*3.793319722393432e+38+t11.*t15.*6.080112186895894e+38+t8.*t20.*7.34340828358964e+38+t15.*t18.*1.216022437379179e+39+t14.*t20.*1.468681656717928e+39+t3.*t8.*t15.*6.281499974462183e+38+t3.*t14.*t15.*6.281499974462183e+38+t2.*t15.*t16.*7.586639444786864e+38+t8.*t14.*t15.*1.216022437379179e+39+t7.*t15.*t16.*1.468681656717928e+39+t12.*t15.*t16.*2.937363313435856e+39).*-8.382656506659224e+38-t369.*t372.*(t22+t24+t27+t28+t30+t32+t35+t36+t38+t41+t42+t44-t6.*t16.*3.244779102448822e+38-t11.*t16.*6.080112186895894e+38-t16.*t18.*1.216022437379179e+39-t3.*t8.*t16.*6.281499974462183e+38-t3.*t14.*t16.*6.281499974462183e+38-t3.*t15.*t16.*3.793319722393432e+38-t8.*t14.*t16.*1.216022437379179e+39-t8.*t15.*t16.*7.34340828358964e+38-t14.*t15.*t16.*1.468681656717928e+39).*8.382656506659224e+38,0.0,0.0],[2,6]);
