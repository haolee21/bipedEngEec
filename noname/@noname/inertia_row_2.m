function Irow = inertia_row_2(rob,in2)
%% INERTIA_ROW_2 - Computation of the robot specific inertia matrix row for corresponding to joint 2 of 6. 
% ========================================================================= 
%    
%    Irow = inertia_row_2(rob,q) 
%    Irow = rob.inertia_row_2(q) 
%    
%  Description:: 
%    Given a full set of joint variables this function computes the 
%    inertia matrix row number 2 of 6 for noname. 
%    
%  Input:: 
%    rob: robot object of noname specific class 
%    q:  6-element vector of generalized 
%         coordinates 
%    Angles have to be given in radians! 
%    
%  Output:: 
%    Irow:  [1x6] row of the robot inertia matrix 
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
%    This is an autogenerated function! 
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
   




%    This function was generated by the Symbolic Math Toolbox version 8.4.
%    03-Apr-2020 15:04:31

q2 = in2(:,2);
q3 = in2(:,3);
q4 = in2(:,4);
q5 = in2(:,5);
t2 = cos(q3);
t3 = cos(q4);
t4 = cos(q5);
t5 = sin(q3);
t6 = sin(q4);
t7 = sin(q5);
t8 = t2.*t3.*(9.0./2.0e+1);
t9 = t2.*t6.*(9.0./2.0e+1);
t10 = t3.*t5.*(9.0./2.0e+1);
t11 = t5.*t6.*(9.0./2.0e+1);
t13 = t2.*t3.*3.3075;
t14 = t2.*t6.*3.3075;
t15 = t3.*t5.*3.3075;
t16 = t5.*t6.*3.3075;
t18 = t2.*2.62458e+1;
t19 = t5.*2.62458e+1;
t20 = t2.*t3.*7.441875e-1;
t21 = t5.*t6.*7.441875e-1;
t23 = t2.*1.181061e+1;
t12 = -t11;
t17 = -t16;
t22 = -t21;
t24 = t9+t10;
t25 = t8+t12+9.0./2.0e+1;
t26 = t4.*t24.*4.326;
t27 = t7.*t24.*4.326;
t29 = t7.*t24.*9.7335e-1;
t28 = -t27;
t30 = -t29;
t31 = t4.*t25.*4.326;
t32 = t7.*t25.*4.326;
t33 = t4.*t25.*9.7335e-1;
t34 = t26+t32;
t38 = t28+t31+9.7335e-1;
t35 = t4.*t34;
t36 = t7.*t34;
t39 = t4.*t38;
t40 = t7.*t38;
t37 = t36.*(9.0./2.0e+1);
t41 = -t40;
t42 = t39.*(9.0./2.0e+1);
t46 = t13+t17+t36+t39+1.65375;
t43 = t14+t15+t35+t41;
t47 = t3.*t46;
t48 = t6.*t46;
t44 = t3.*t43;
t45 = t6.*t43;
t49 = -t48;
t50 = t19+t44+t49;
t52 = t18+t45+t47+2.62458e+1;
t51 = t5.*t50.*(9.0./2.0e+1);
t53 = t2.*t52.*(9.0./2.0e+1);
Irow = [t20+t22+t23+t30+t33+t37+t42+t51+t53+cos(q2).*(t2.*t52+t5.*t50+1.65375).*(9.0./2.0e+1)+sin(q2).*(t2.*t50-t5.*t52).*(9.0./2.0e+1)+1.7031735e+1,t20+t22+t23+t30+t33+t37+t42+t51+t53+1.7031735e+1,t20+t22+t23+t30+t33+t37+t42+1.653561e+1,t20+t22+t30+t33+t37+t42+7.8813e-1,t30+t33+2.92005e-1,0.0];
