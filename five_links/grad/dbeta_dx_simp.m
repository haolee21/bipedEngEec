%% Generate function for d(beta)/dx


syms q1 q2 q3 q4 q5 dq1 dq2 dq3 dq4 dq5

M = five_M_simp(q2,q3,q4,q5);
V = five_V(q2,q3,q4,q5,dq1,dq2,dq3,dq4,dq5)*[dq1;dq2;dq3;dq4;dq5];
G = five_G(q1,q2,q3,q4,q5);


dV_dx=simplify([diff(V,q1),diff(V,q2),diff(V,q3),diff(V,q4),diff(V,q5),diff(V,dq1),diff(V,dq2),diff(V,dq3),diff(V,dq4),diff(V,dq5),zeros(5,5)].'); %transpose because of den layout

dG_dx = [diff(G,q1);diff(G,q2);diff(G,q3);diff(G,q4);diff(G,q5);zeros(10,5)].';

dM_dx1 = zeros(5);
dM_dx2 = diff(M,q2);
dM_dx3 = diff(M,q3);
dM_dx4 = diff(M,q4);
dM_dx5 = diff(M,q5);
dM_dx6 = zeros(5);
dM_dx7 = zeros(5);
dM_dx8 = zeros(5);
dM_dx9 = zeros(5);
dM_dx10 = zeros(5);
dM_dx11 = zeros(5);
dM_dx12 = zeros(5);
dM_dx13 = zeros(5);
dM_dx14 = zeros(5);
dM_dx15 = zeros(5);

% matlabFunction(dV_dx,'file','dV_dx');
% matlabFunction(dG_dx,'file','dG_dx');
matlabFunction(dM_dx2,'file','dM_dx2');
matlabFunction(dM_dx3,'file','dM_dx3');
matlabFunction(dM_dx4,'file','dM_dx4');
matlabFunction(dM_dx5,'file','dM_dx5');