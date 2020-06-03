function Crow = coriolis_row_1(rob,in2,in3)
%% CORIOLIS_ROW_1 - Computation of the robot specific Coriolis matrix row for joint 1 of 5. 
% ========================================================================= 
%    
%    Crow = coriolis_row_1(rob,q,qd) 
%    Crow = rob.coriolis_row_1(q,qd) 
%    
%  Description:: 
%    Given a full set of joint variables and their first order temporal derivatives this function computes the 
%    Coriolis matrix row number 1 of 5 for noname. 
%    
%  Input:: 
%    rob: robot object of noname specific class 
%    qd:  5-element vector of generalized 
%    q:  5-element vector of generalized 
%    
%  Output:: 
%    Crow:  [1x5] row of the robot Coriolis matrix 
%    
%  Example:: 
%    --- 
%    
%  Known Bugs:: 
%    --- 
%    
%  TODO:: 
%    --- 
%    
%  References:: 
%    1) Robot Modeling and Control - Spong, Hutchinson, Vidyasagar 
%    2) Modelling and Control of Robot Manipulators - Sciavicco, Siciliano 
%    3) Introduction to Robotics, Mechanics and Control - Craig 
%    4) Modeling, Identification & Control of Robots - Khalil & Dombre 
%    
%  Authors:: 
%    This is an autogenerated function. 
%    Code generator written by: 
%    Joern Malzahn 
%    2012 RST, Technische Universitaet Dortmund, Germany 
%    http://www.rst.e-technik.tu-dortmund.de 
%    
%  See also coriolis.
%    
    
% Copyright (C) 1993-2020, by Peter I. Corke 
% Copyright (C) 2012-2020, by Joern Malzahn 
% 
% This file has been automatically generated with The Robotics Toolbox for Matlab (RTB). 
% 
% RTB and code generated with RTB is free software: you can redistribute it and/or modify 
% it under the terms of the GNU Lesser General Public License as published by 
% the Free Software Foundation, either version 3 of the License, or 
% (at your option) any later version. 
%  
% RTB is distributed in the hope that it will be useful, 
% but WITHOUT ANY WARRANTY; without even the implied warranty of 
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
% GNU Lesser General Public License for more details. 
%  
% You should have received a copy of the GNU Leser General Public License 
% along with RTB.  If not, see <http://www.gnu.org/licenses/>. 
% 
% http://www.petercorke.com 
% 
% The code generation module emerged during the work on a project funded by 
% the German Research Foundation (DFG, BE1569/7-1). The authors gratefully  
% acknowledge the financial support. 

%% Bugfix
%  In some versions the symbolic toolbox writes the constant $pi$ in
%  capital letters. This way autogenerated functions might not work properly.
%  To fix this issue a local variable is introduced:
PI = pi;
   




%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    22-Apr-2020 01:41:29

q2 = in2(:,2);
q3 = in2(:,3);
q4 = in2(:,4);
q5 = in2(:,5);
qd1 = in3(:,1);
qd2 = in3(:,2);
qd3 = in3(:,3);
qd4 = in3(:,4);
qd5 = in3(:,5);
t2 = sin(q2);
t3 = sin(q3);
t4 = sin(q5);
t5 = q2+q3;
t6 = q3+q4;
t7 = q4+t5;
t8 = q5+t6;
t9 = sin(t5);
t10 = sin(t6);
t15 = qd2.*t2.*1.49191875e+1;
t17 = qd5.*t4.*4.380075e-1;
t19 = qd3.*t3.*1.181061e+1;
t11 = q5+t7;
t12 = sin(t7);
t13 = sin(t8);
t16 = -t15;
t18 = -t17;
t20 = -t19;
t21 = qd1.*t10.*1.6202025;
t22 = qd2.*t10.*1.6202025;
t23 = qd3.*t10.*1.6202025;
t24 = qd4.*t10.*1.6202025;
t25 = qd1.*t9.*1.181061e+1;
t26 = qd2.*t9.*1.181061e+1;
t27 = qd3.*t9.*1.181061e+1;
t14 = sin(t11);
t28 = qd1.*t13.*4.380075e-1;
t29 = qd2.*t13.*4.380075e-1;
t30 = qd3.*t13.*4.380075e-1;
t31 = qd4.*t13.*4.380075e-1;
t32 = qd5.*t13.*4.380075e-1;
t38 = -t21;
t39 = -t22;
t40 = -t23;
t41 = -t24;
t42 = -t25;
t43 = -t26;
t44 = -t27;
t50 = qd1.*t12.*1.6202025;
t51 = qd2.*t12.*1.6202025;
t52 = qd3.*t12.*1.6202025;
t53 = qd4.*t12.*1.6202025;
t33 = qd1.*t14.*4.380075e-1;
t34 = qd2.*t14.*4.380075e-1;
t35 = qd3.*t14.*4.380075e-1;
t36 = qd4.*t14.*4.380075e-1;
t37 = qd5.*t14.*4.380075e-1;
t45 = -t28;
t46 = -t29;
t47 = -t30;
t48 = -t31;
t49 = -t32;
t59 = -t50;
t60 = -t51;
t61 = -t52;
t62 = -t53;
t54 = -t33;
t55 = -t34;
t56 = -t35;
t57 = -t36;
t58 = -t37;
Crow = [t16+t18+t20+t40+t41+t43+t44+t47+t48+t49+t55+t56+t57+t58+t60+t61+t62,t16+t18+t20+t40+t41+t42+t43+t44+t47+t48+t49+t54+t55+t56+t57+t58+t59+t60+t61+t62-qd1.*t2.*1.49191875e+1,t18+t20+t38+t39+t40+t41+t42+t43+t44+t45+t46+t47+t48+t49+t54+t55+t56+t57+t58+t59+t60+t61+t62-qd1.*t3.*1.181061e+1-qd2.*t3.*1.181061e+1,t18+t38+t39+t40+t41+t45+t46+t47+t48+t49+t54+t55+t56+t57+t58+t59+t60+t61+t62,(t4+t13+t14).*(qd1+qd2+qd3+qd4+qd5).*(-4.380075e-1)];
