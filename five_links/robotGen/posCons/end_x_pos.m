function end_x_pos = end_x_pos(q1,q2,q3,q4,q5)
%END_X_POS
%    END_X_POS = END_X_POS(Q1,Q2,Q3,Q4,Q5)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    21-Apr-2020 22:27:04

t2 = q1+q2;
t3 = q3+q4+t2;
end_x_pos = cos(q5+t3).*(9.0./2.0e+1)+cos(q1).*(9.0./2.0e+1)+cos(t2).*(9.0./2.0e+1)+cos(t3).*(9.0./2.0e+1);
