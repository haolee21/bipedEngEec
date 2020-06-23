function [tau_ext,dTau_ext,Fn,Fs,dFn,dFs] = toe_grf(x,p)
x = x(1,:); %we leave for future delay function

J = six_J(x(1),x(2),x(3),x(4),x(5),x(6));
dJdx1 = dJ_q1(x);
dJdx2 = dJ_q2(x);
dJdx3 = dJ_q3(x);
dJdx4 = dJ_q4(x);
dJdx5 = dJ_q5(x);
dJdx6 = dJ_q6(x);




Fn = Fn_toe(x,p.toe_th,p.k,p.cmax,p.dmax,p.us,p.ud);
Fs = -Fs_toe(x,p.toe_th,p.k,p.cmax,p.dmax,p.us,p.ud);
dFs = -dFs_toe(x,p.toe_th,p.k,p.cmax,p.dmax,p.us,p.ud);
dFn = dFn_toe(x,p.toe_th,p.k,p.cmax,p.dmax,p.us,p.ud);
tau_ext = [Fs,Fn,zeros(1,4)]*J;


Fext_diag = repmat({[Fs,Fn,zeros(1,4)]},1,3*p.numJ);
dTau_ext = real(blkdiag(Fext_diag{:})*[dJdx1;dJdx2;dJdx3;dJdx4;dJdx5;dJdx6;zeros(p.numJ*2*p.numJ,p.numJ)]...
        +[dFs,dFn,zeros(2*p.numJ,4);zeros(p.numJ)]*J);



end