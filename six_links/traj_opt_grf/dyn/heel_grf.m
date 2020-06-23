function [tau_ext,dTau_ext,Fn,Fs,dFn,dFs] = heel_grf(x,p)
x = x(1,:); %we leave for future delay function

J2 = six_J2(x(1),x(2),x(3),x(4),x(5),x(6));
dJdx1 = dJ2_q1(x);
dJdx2 = dJ2_q2(x);
dJdx3 = dJ2_q3(x);
dJdx4 = dJ2_q4(x);
dJdx5 = dJ2_q5(x);
dJdx6 = dJ2_q6(x);



Fn = Fn_heel(x,p.toe_th,p.k,p.cmax,p.dmax,p.us,p.ud);    
dFn = [dFn_heel1(x,p.toe_th,p.k,p.cmax,p.dmax,p.us,p.ud);
       dFn_heel2(x,p.toe_th,p.k,p.cmax,p.dmax,p.us,p.ud);
       dFn_heel3(x,p.toe_th,p.k,p.cmax,p.dmax,p.us,p.ud);
       dFn_heel4(x,p.toe_th,p.k,p.cmax,p.dmax,p.us,p.ud);
       dFn_heel5(x,p.toe_th,p.k,p.cmax,p.dmax,p.us,p.ud);
       dFn_heel6(x,p.toe_th,p.k,p.cmax,p.dmax,p.us,p.ud);
       dFn_heel7(x,p.toe_th,p.k,p.cmax,p.dmax,p.us,p.ud);
       dFn_heel8(x,p.toe_th,p.k,p.cmax,p.dmax,p.us,p.ud);
       dFn_heel9(x,p.toe_th,p.k,p.cmax,p.dmax,p.us,p.ud);
       dFn_heel10(x,p.toe_th,p.k,p.cmax,p.dmax,p.us,p.ud);
       dFn_heel11(x,p.toe_th,p.k,p.cmax,p.dmax,p.us,p.ud);
       dFn_heel12(x,p.toe_th,p.k,p.cmax,p.dmax,p.us,p.ud)];
Fs = Fs_heel(x,p.toe_th,p.k,p.cmax,p.dmax,p.us,p.ud);
dFs = [dFs_heel1(x,p.toe_th,p.k,p.cmax,p.dmax,p.us,p.ud);
       dFs_heel2(x,p.toe_th,p.k,p.cmax,p.dmax,p.us,p.ud);
       dFs_heel3(x,p.toe_th,p.k,p.cmax,p.dmax,p.us,p.ud);
       dFs_heel4(x,p.toe_th,p.k,p.cmax,p.dmax,p.us,p.ud);
       dFs_heel5(x,p.toe_th,p.k,p.cmax,p.dmax,p.us,p.ud);
       dFs_heel6(x,p.toe_th,p.k,p.cmax,p.dmax,p.us,p.ud);
       dFs_heel7(x,p.toe_th,p.k,p.cmax,p.dmax,p.us,p.ud);
       dFs_heel8(x,p.toe_th,p.k,p.cmax,p.dmax,p.us,p.ud);
       dFs_heel9(x,p.toe_th,p.k,p.cmax,p.dmax,p.us,p.ud);
       dFs_heel10(x,p.toe_th,p.k,p.cmax,p.dmax,p.us,p.ud);
       dFs_heel11(x,p.toe_th,p.k,p.cmax,p.dmax,p.us,p.ud);
       dFs_heel12(x,p.toe_th,p.k,p.cmax,p.dmax,p.us,p.ud)];
tau_ext = [-Fs,Fn,zeros(1,4)]*J2;
Fext_diag = repmat({[-Fs,Fn,zeros(1,4)]},1,3*p.numJ);
dTau_ext = real(blkdiag(Fext_diag{:})*[dJdx1;dJdx2;dJdx3;dJdx4;dJdx5;dJdx6;zeros(p.numJ*2*p.numJ,p.numJ)]...
          +[-dFs,dFn,zeros(2*p.numJ,4);zeros(p.numJ)]*J2);



end
