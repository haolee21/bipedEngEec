function out1 = six_M(q2,q3,q4,q5,q6)
%SIX_M
%    OUT1 = SIX_M(Q2,Q3,Q4,Q5,Q6)

%    This function was generated by the Symbolic Math Toolbox version 8.4.
%    27-May-2020 01:00:00

t2 = cos(q2);
t3 = cos(q3);
t4 = cos(q5);
t5 = cos(q6);
t6 = q2+q3;
t7 = q3+q4;
t8 = q5+q6;
t9 = cos(t6);
t10 = cos(t7);
t11 = cos(t8);
t12 = q4+t6;
t13 = q5+t7;
t17 = t7+t8;
t23 = t5.*1.09078504542792e-1;
t24 = t5.*5.4539252271396e-2;
t25 = t2.*1.163815784709457e+1;
t29 = t3.*7.571420887311746;
t30 = t3.*1.514284177462349e+1;
t37 = t4.*4.823752644884052e-1;
t38 = t4.*9.647505289768104e-1;
t14 = cos(t12);
t15 = cos(t13);
t16 = q5+t12;
t19 = cos(t17);
t20 = t8+t12;
t22 = t9.*7.822204368661843;
t31 = t10.*1.702413339329298;
t32 = t10.*3.404826678658597;
t33 = t11.*5.279069867828761e-2;
t34 = t11.*1.055813973565752e-1;
t41 = t24+1.413712088292136e-2;
t18 = cos(t16);
t21 = cos(t20);
t26 = t14.*1.75880132122703;
t35 = t19.*5.279069867828761e-2;
t36 = t19.*1.055813973565752e-1;
t39 = t15.*4.823752644884052e-1;
t40 = t15.*9.647505289768104e-1;
t42 = t33+t41;
t44 = t23+t34+t38+1.551409756943671;
t45 = t23+t33+t37+3.492238326654392e-1;
t27 = t21.*5.4539252271396e-2;
t28 = t18.*4.98352681402092e-1;
t43 = t35+t42;
t47 = t35+t39+t45;
t48 = t31+t35+t39+t44;
t50 = t23+t29+t31+t34+t35+t38+t39+9.456211101891361;
t46 = t27+t43;
t49 = t27+t28+t47;
t51 = t26+t27+t28+t48;
t52 = t22+t26+t27+t28+t50;
t53 = t22+t23+t25+t26+t27+t28+t30+t32+t34+t36+t38+t40+2.022101682332396e+1;
out1 = reshape([t2.*2.327631569418914e+1+t9.*1.564440873732369e+1+t14.*3.517602642454061+t18.*9.96705362804184e-1+t21.*1.09078504542792e-1+t23+t30+t32+t34+t36+t38+t40+3.340009651583829e+1,t53,t52,t51,t49,t46,t53,t23+t30+t32+t34+t36+t38+t40+2.022101682332396e+1,t50,t48,t47,t43,t52,t50,t23+t34+t38+9.456211101891361,t44,t45,t42,t51,t48,t44,t44,t45,t42,t49,t47,t45,t45,t23+3.492238326654392e-1,t41,t46,t43,t42,t42,t41,1.413712088292136e-2],[6,6]);
