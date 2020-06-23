function [c,ceq,gradc,gradceq]=grf_cons(x,p)

if(~isnumeric(x))
    
    ceq_Fy_toe = sym(zeros(size(x,2),1));
    gradceq_Fy_toe = sym(zeros(size(x,1),size(x,2),size(x,2)));
    ceq_Fy_heel = sym(zeros(size(x,2),1));
    gradceq_Fy_heel = sym(zeros(size(x,1),size(x,2),size(x,2)));
    c_Fs_toe = sym(zeros(size(x,2),1));
    c_Fs_heel = sym(zeros(size(x,2),1));
    gradc_Fs_toe = sym(zeros(size(x,1),size(x,2),size(x,2)));
    gradc_Fs_heel = sym(zeros(size(x,1),size(x,2),size(x,2)));
else
    
    ceq_Fy_toe = zeros(size(x,2),1);
    gradceq_Fy_toe = zeros(size(x,1),size(x,2),size(x,2));
    ceq_Fy_heel = zeros(size(x,2),1);
    gradceq_Fy_heel = zeros(size(x,1),size(x,2),size(x,2));
    c_Fs_toe = zeros(size(x,2),1);
    c_Fs_heel = zeros(size(x,2),1);
    gradc_Fs_toe = zeros(size(x,1),size(x,2),size(x,2));
    gradc_Fs_heel = zeros(size(x,1),size(x,2),size(x,2));
    
end

for i=1:size(x,2)
    cur_s = x(1:p.numJ*2,i).';
    cur_F_toe = x(p.numJ*3+1:p.numJ*3+2,i).';
    cur_F_heel = x(p.numJ*3+3:p.numJ*3+4,i).';
    ceq_Fy_toe(i,1) = Fn_toe_ceq1(cur_s,cur_F_toe,p.toe_th,p.k,p.cmax,p.dmax);
    ceq_Fy_heel(i,1) = Fn_heel_ceq1(cur_s,cur_F_heel,p.toe_th,p.k,p.cmax,p.dmax);
    gradceq_Fy_toe(:,i,i) = dFn_toe_ceq1(cur_s,cur_F_toe,p.toe_th,p.k,p.cmax,p.dmax);
    gradceq_Fy_heel(:,i,i) = dFn_heel_ceq1(cur_s,cur_F_heel,p.toe_th,p.k,p.cmax,p.dmax);
    
    c_Fs_toe(i,1) = Fs_toe_c1(cur_s,cur_F_toe,p.us,p.ud);
    c_Fs_heel(i,1) = Fs_heel_c1(cur_s,cur_F_heel,p.us,p.ud);
    
    gradc_Fs_toe(:,i,i) = [dFs_toe_c1_1(cur_s,cur_F_toe,p.us,p.ud);
                           dFs_toe_c1_2(cur_s,cur_F_toe,p.us,p.ud);
                           dFs_toe_c1_3(cur_s,cur_F_toe,p.us,p.ud);
                           dFs_toe_c1_4(cur_s,cur_F_toe,p.us,p.ud);
                           dFs_toe_c1_5(cur_s,cur_F_toe,p.us,p.ud);
                           dFs_toe_c1_6(cur_s,cur_F_toe,p.us,p.ud);
                           dFs_toe_c1_7(cur_s,cur_F_toe,p.us,p.ud);
                           dFs_toe_c1_8(cur_s,cur_F_toe,p.us,p.ud);
                           dFs_toe_c1_9(cur_s,cur_F_toe,p.us,p.ud);
                           dFs_toe_c1_10(cur_s,cur_F_toe,p.us,p.ud);
                           dFs_toe_c1_11(cur_s,cur_F_toe,p.us,p.ud);
                           dFs_toe_c1_12(cur_s,cur_F_toe,p.us,p.ud);
                           zeros(p.numJ,1);
                           dFs_toe_c1_19(cur_s,cur_F_toe,p.us,p.ud);
                           dFs_toe_c1_20(cur_s,cur_F_toe,p.us,p.ud);
                           zeros(2,1)];
                       
    gradc_Fs_heel(:,i,i) = [dFs_heel_c1_1(cur_s,cur_F_heel,p.us,p.ud);
                           dFs_heel_c1_2(cur_s,cur_F_heel,p.us,p.ud);
                           dFs_heel_c1_3(cur_s,cur_F_heel,p.us,p.ud);
                           dFs_heel_c1_4(cur_s,cur_F_heel,p.us,p.ud);
                           dFs_heel_c1_5(cur_s,cur_F_heel,p.us,p.ud);
                           dFs_heel_c1_6(cur_s,cur_F_heel,p.us,p.ud);
                           dFs_heel_c1_7(cur_s,cur_F_heel,p.us,p.ud);
                           dFs_heel_c1_8(cur_s,cur_F_heel,p.us,p.ud);
                           dFs_heel_c1_9(cur_s,cur_F_heel,p.us,p.ud);
                           dFs_heel_c1_10(cur_s,cur_F_heel,p.us,p.ud);
                           dFs_heel_c1_11(cur_s,cur_F_heel,p.us,p.ud);
                           dFs_heel_c1_12(cur_s,cur_F_heel,p.us,p.ud);
                           zeros(p.numJ,1);
                           zeros(2,1);
                           dFs_heel_c1_21(cur_s,cur_F_heel,p.us,p.ud);
                           dFs_heel_c1_22(cur_s,cur_F_heel,p.us,p.ud)];
                         
           
                       
    
    
end
c = [c_Fs_toe;c_Fs_heel];
gradc = [reshape(gradc_Fs_toe,[size(x,1)*size(x,2),size(x,2)]),reshape(gradc_Fs_heel,[size(x,1)*size(x,2),size(x,2)])];

ceq = [ceq_Fy_toe;ceq_Fy_heel];
gradceq = [reshape(gradceq_Fy_toe,[size(x,1)*size(x,2),size(x,2)]),reshape(gradceq_Fy_heel,[size(x,1)*size(x,2),size(x,2)])];


end