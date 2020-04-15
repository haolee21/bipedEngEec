%% this script will calculate the dynamic equation of the model

syms q1 q2 q3 q4 q5 dq1 dq2 dq3 dq4 dq5 real
syms fext_x fext_y fext_z text_1 text_2 text_3 real%external force/torque
syms t1 t2 t3 t4 t5 real %joint torque


M = five_M(q2,q3,q4,q5);
G = five_G(q1,q2,q3,q4,q5);
V = five_V(q2,q3,q4,q5,dq1,dq2,dq3,dq4,dq5);
J = five_J(q1,q2,q3,q4,q5);
tau = [t1,t2,t3,t4,t5].';
ddq = M\(tau-V-G.'-J.'*[fext_x, fext_y, fext_z,text_1, text_2, text_3].');