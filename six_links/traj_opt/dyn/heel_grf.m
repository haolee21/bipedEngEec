function [tau_ext,dTau_ext,Fn,Fs,dFn,dFs] = heel_grf(x,p)
x = x(1,:); %we leave for future delay function
ypos = heel_y_pos(x);
J2 = six_J2(x(1),x(2),x(3),x(4),x(5),x(6));
dJdx1 = dJ2_q1(x);
dJdx2 = dJ2_q2(x);
dJdx3 = dJ2_q3(x);
dJdx4 = dJ2_q4(x);
dJdx5 = dJ2_q5(x);
dJdx6 = dJ2_q6(x);
end_vel = x(p.numJ+1:2*p.numJ)*J2.';

fs_dir =1;

    
if(p.toe_th-ypos<p.dmax)
    Fn = Fn_heel1(x,p.toe_th);
    dFn = dFn_heel1(x,p.toe_th);
    if(abs(end_vel(1))<0.05)
        Fs = Fs_heel1_s(x,p.toe_th);
        dFs = dFs_heel1_s(x,p.toe_th);
    else
        Fs = Fs_heel1_d(x,p.toe_th);
        dFs = dFs_heel1_d(x,p.toe_th);
    end
else
    Fn = Fn_heel2(x,p.toe_th);
    dFn = dFn_heel2(x,p.toe_th);
    if(abs(end_vel(1))<0.05)
        Fs = Fs_heel2_s(x,p.toe_th);
        dFs = dFs_heel2_s(x,p.toe_th);
    else
        Fs = Fs_heel2_d(x,p.toe_th);
        dFs = dFs_heel2_d(x,p.toe_th);
    end
end
if(Fn>0)
    if(end_vel(1)>0)
        fs_dir=-1;
    end
    tau_ext = [fs_dir*Fs,Fn,zeros(1,4)]*J2;
    
    Fext_diag = repmat({[fs_dir*Fs,Fn,zeros(1,4)]},1,3*p.numJ);
    dTau_ext = real(blkdiag(Fext_diag{:})*[dJdx1;dJdx2;dJdx3;dJdx4;dJdx5;dJdx6;zeros(p.numJ*2*p.numJ,p.numJ)]...
        +[fs_dir*dFs,dFn,zeros(2*p.numJ,4);zeros(p.numJ)]*J2);
else
    tau_ext = zeros(1,p.numJ);
    dTau_ext = zeros(3*p.numJ,p.numJ);
    Fn =0;
    Fs=0;
    dFn = zeros(p.numJ*2,1);
    dFs = zeros(p.numJ*2,1);
end

Fs = fs_dir*Fs;

end
