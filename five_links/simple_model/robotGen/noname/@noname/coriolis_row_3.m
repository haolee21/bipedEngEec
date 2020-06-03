function Crow = coriolis_row_3(rob,in2,in3)
%% CORIOLIS_ROW_3 - Computation of the robot specific Coriolis matrix row for joint 3 of 5. 
% ========================================================================= 
%    
%    Crow = coriolis_row_3(rob,q,qd) 
%    Crow = rob.coriolis_row_3(q,qd) 
%    
%  Description:: 
%    Given a full set of joint variables and their first order temporal derivatives this function computes the 
%    Coriolis matrix row number 3 of 5 for noname. 
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
%    22-Apr-2020 01:41:30

q2 = in2(:,2);
q3 = in2(:,3);
q4 = in2(:,4);
q5 = in2(:,5);
qd1 = in3(:,1);
qd2 = in3(:,2);
qd3 = in3(:,3);
qd4 = in3(:,4);
qd5 = in3(:,5);
t2 = sin(q3);
t3 = sin(q5);
t4 = q3+q4;
t5 = q5+t4;
t6 = sin(t4);
t8 = qd5.*t3.*4.380075e-1;
t10 = qd1.*t2.*1.181061e+1;
t11 = qd2.*t2.*1.181061e+1;
t7 = sin(t5);
t9 = -t8;
t12 = qd1.*t6.*1.6202025;
t13 = qd2.*t6.*1.6202025;
t14 = qd1.*t7.*4.380075e-1;
t15 = qd2.*t7.*4.380075e-1;
Crow = [t9+t10+t11+t12+t13+t14+t15+qd1.*sin(q2+q3).*1.181061e+1+qd1.*sin(q2+t4).*1.6202025+qd1.*sin(q2+t5).*4.380075e-1,t9+t10+t11+t12+t13+t14+t15,t9,t9,t3.*(qd1+qd2+qd3+qd4+qd5).*(-4.380075e-1)];
