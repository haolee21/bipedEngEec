function Crow = coriolis_row_4(rob,in2,in3)
%% CORIOLIS_ROW_4 - Computation of the robot specific Coriolis matrix row for joint 4 of 6. 
% ========================================================================= 
%    
%    Crow = coriolis_row_4(rob,q,qd) 
%    Crow = rob.coriolis_row_4(q,qd) 
%    
%  Description:: 
%    Given a full set of joint variables and their first order temporal derivatives this function computes the 
%    Coriolis matrix row number 4 of 6 for noname. 
%    
%  Input:: 
%    rob: robot object of noname specific class 
%    qd:  6-element vector of generalized 
%    q:  6-element vector of generalized 
%    
%  Output:: 
%    Crow:  [1x6] row of the robot Coriolis matrix 
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
   




%    This function was generated by the Symbolic Math Toolbox version 8.4.
%    03-Apr-2020 15:11:08

q2 = in2(:,2);
q3 = in2(:,3);
q4 = in2(:,4);
q5 = in2(:,5);
qd1 = in3(:,1);
qd2 = in3(:,2);
qd3 = in3(:,3);
qd4 = in3(:,4);
qd5 = in3(:,5);
t2 = sin(q5);
t3 = q3+q4;
t4 = q5+t3;
t5 = sin(t3);
t7 = qd5.*t2.*4.380075e-1;
t6 = sin(t4);
t8 = -t7;
t9 = qd1.*t5.*1.6202025;
t10 = qd2.*t5.*1.6202025;
t11 = qd1.*t6.*4.380075e-1;
t12 = qd2.*t6.*4.380075e-1;
Crow = [t8+t9+t10+t11+t12+qd1.*sin(q2+t3).*1.6202025+qd1.*sin(q2+t4).*4.380075e-1,t8+t9+t10+t11+t12,t8,t8,t2.*(qd1+qd2+qd3+qd4+qd5).*(-4.380075e-1),0.0];
