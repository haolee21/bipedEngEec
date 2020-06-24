function Hout = hessianfcn(x,lambda,p)
Hout = zeros(size(x,1));
ineq_idx =1;
eq_idx=1;

%hipvel_neq constraints
for i=1:floor(p.gaitT/p.sampT)+1 %hipVelCon, total 101 constraints, will be affexted by 
    curHout = zeros(size(x,1));
    curX = x((i-1)*(p.numJ*3+4)+1:i*(p.numJ*3+4));
    cur_h=h_hipVelCon(curX.');
    cur_h=cur_h+cur_h.'-diag(diag(cur_h));
    
    curHout((i-1)*(p.numJ*3+4)+1:(i-1)*(p.numJ*3+4)+p.numJ*3,(i-1)*(p.numJ*3+4)+1:(i-1)*(p.numJ*3+4)+p.numJ*3)=cur_h;
    Hout = Hout+lambda.ineqnonlin(ineq_idx)*curHout;
    ineq_idx = ineq_idx+1;
end

%grf_neq constraints
for i=1:floor(p.gaitT/p.sampT)+1
    %toe constraints
    curHout = zeros(size(x,1));
    curX = x((i-1)*(p.numJ*3+4)+1:i*(p.numJ*3+4));
    cur_h = hess_grf_c_toe(curX.',p.toe_th,p.dmax,p.cmax,p.k,p.us,p.ud);
    curHout((i-1)*(p.numJ*3+4)+1:i*(p.numJ*3+4),(i-1)*(p.numJ*3+4)+1:i*(p.numJ*3+4))=cur_h;
    Hout = Hout+lambda.ineqnonlin(ineq_idx)*curHout;
    ineq_idx = ineq_idx+1; 
    %heel constraints
    curHout = zeros(size(x,1));
    cur_h = hess_grf_c_heel(curX.',p.toe_th,p.dmax,p.cmax,p.k,p.us,p.ud);
    curHout((i-1)*(p.numJ*3+4)+1:i*(p.numJ*3+4),(i-1)*(p.numJ*3+4)+1:i*(p.numJ*3+4))=cur_h;
    Hout = Hout+lambda.ineqnonlin(ineq_idx)*curHout;
    ineq_idx = ineq_idx+1; 
end

%eq constraints
% dynamic constraints
for i=1:floor(p.gaitT/p.sampT)%iterate for timeSteps -1
    x1 = x((i-1)*(p.numJ*3+4)+1:i*(p.numJ*3+4));
    x2 = x(i*(p.numJ*3+4)+1:(i+1)*(p.numJ*3+4));
    [f1,df1s,df1u,df1f] = f_x2(x1,p);
    [f2,df2s,df2u,df2f] = f_x2(x2,p);
%     s_half = 0.5*(x1(1:p.numJ*2)+x2(1:p.numJ*2))+p.sampT/8*(f1-f2);
%     u_half = 0.5*(x1(2*p.numJ+1:3*p.numJ)+x2(2*p.numJ+1:3*p.numJ));
%     F_half = 0.5*(x1(3*p.numJ+1:3*p.numJ+4)+x2(3*p.numJ+1:3*p.numJ+4));
%     x_half = [s_half;u_half;F_half];
%     [~,df_half_s,~,~] = f_x2(x_half,p);
    
    
    % I will probably cause confusion here, but here f's dimension should
    % be p.numJ, since we only have 1/M for ddq, and dq's hessian is zero
    dM_dxx1 = dM_dxx(x1.');
    dM_dxx2 = dM_dxx(x2.');
    dMdxx_f1 = zeros(p.numJ,p.numJ,p.numJ);
    dMdxx_f2 = zeros(p.numJ,p.numJ,p.numJ);
    for it1=1:p.numJ
        for it2=1:p.numJ
            dMdxx_f1(it1,it2,:)=reshape(dM_dxx1(it1,it2,:,:),[p.numJ,p.numJ])*f1(p.numJ+1:p.numJ*2,:); %I only solve dM_dxx for theta, but we actually have qdot, yet they will just be zeros
            dMdxx_f2(it1,it2,:)=reshape(dM_dxx2(it1,it2,:,:),[p.numJ,p.numJ])*f2(p.numJ+1:p.numJ*2,:);
        end
    end
    
    
    % solve the hessian of M*ddq (since we cannot get 1/M symbolically)
     
    Mddq_d1= Mddq_dxx(x1.',p.ank_stiff,p.knee_stiff);
    Mddq_d2= Mddq_dxx(x2.',p.ank_stiff,p.knee_stiff);
    
%     Mddq_dxx1 = Mddq_d1(1:2*p.numJ,1:2*p.numJ,:);
%     Mddq_dxx2 = Mddq_d2(1:2*p.numJ,1:2*p.numJ,:);
%     
%     Mddq_dux1 = Mddq_d1(2*p.numJ+1:3*p.numJ,1:2*p.numJ,:);
%     Mddq_dux2 = Mddq_d2(2*p.numJ+1:3*p.numJ,1:2*p.numJ,:);
%     
%     Mddq_dFx1 = Mddq_d1(3*p.numJ+1:3*p.numJ+4,1:2*p.numJ,:);
%     Mddq_dFx2 = Mddq_d2(3*p.numJ+1:3*p.numJ+4,1:2*p.numJ,:);
    
    dMdx2_1 = dM_dx2(x1(2),x1(3),x1(4),x1(5),x1(6));
    dMdx3_1 = dM_dx3(x1(2),x1(3),x1(4),x1(5),x1(6));
    dMdx4_1 = dM_dx4(x1(2),x1(3),x1(4),x1(5),x1(6));
    dMdx5_1 = dM_dx5(x1(2),x1(3),x1(4),x1(5),x1(6));
    dMdx6_1 = dM_dx6(x1(2),x1(3),x1(4),x1(5),x1(6));
    
    dMdx2_2 = dM_dx2(x2(2),x2(3),x2(4),x2(5),x2(6));
    dMdx3_2 = dM_dx3(x2(2),x2(3),x2(4),x2(5),x2(6));
    dMdx4_2 = dM_dx4(x2(2),x2(3),x2(4),x2(5),x2(6));
    dMdx5_2 = dM_dx5(x2(2),x2(3),x2(4),x2(5),x2(6));
    dMdx6_2 = dM_dx6(x2(2),x2(3),x2(4),x2(5),x2(6));
    
    %% solve the cross terms (dM/dx * df/dx)
    % former one
    cross_term1 = zeros(p.numJ,p.numJ,p.numJ);
    cross_term2 = zeros(p.numJ,p.numJ,p.numJ);
    cross_term1(2,:,:) = dMdx2_1*df1s(p.numJ+1:2*p.numJ,p.numJ+1:2*p.numJ).';
    cross_term1(3,:,:) = dMdx3_1*df1s(p.numJ+1:2*p.numJ,p.numJ+1:2*p.numJ).';
    cross_term1(4,:,:) = dMdx4_1*df1s(p.numJ+1:2*p.numJ,p.numJ+1:2*p.numJ).';
    cross_term1(5,:,:) = dMdx5_1*df1s(p.numJ+1:2*p.numJ,p.numJ+1:2*p.numJ).';
    cross_term1(6,:,:) = dMdx6_1*df1s(p.numJ+1:2*p.numJ,p.numJ+1:2*p.numJ).';
    
    cross_term2(2,:,:) = dMdx2_2*df2s(p.numJ+1:2*p.numJ,p.numJ+1:2*p.numJ).';
    cross_term2(3,:,:) = dMdx3_2*df2s(p.numJ+1:2*p.numJ,p.numJ+1:2*p.numJ).';
    cross_term2(4,:,:) = dMdx4_2*df2s(p.numJ+1:2*p.numJ,p.numJ+1:2*p.numJ).';
    cross_term2(5,:,:) = dMdx5_2*df2s(p.numJ+1:2*p.numJ,p.numJ+1:2*p.numJ).';
    cross_term2(6,:,:) = dMdx6_2*df2s(p.numJ+1:2*p.numJ,p.numJ+1:2*p.numJ).';
    
    cross_term1 = permute(cross_term1,[3,1,2]);
    cross_term2 = permute(cross_term2,[3,1,2]);
    df_dss1 = zeros(3*p.numJ+4,3*p.numJ+4,p.numJ);
    df_dss2 = zeros(3*p.numJ+4,3*p.numJ+4,p.numJ);
    df_dss1(1:p.numJ,1:p.numJ,:) = dMdxx_f1+cross_term1;  %in the equations these are positive, but later I will subtract them
    df_dss2(1:p.numJ,1:p.numJ,:) = dMdxx_f2+cross_term2;
    
    % later one
%     cross_term1 = zeros(p.numJ,p.numJ,p.numJ);
%     cross_term2 = zeros(p.numJ,p.numJ,p.numJ);
%     cross_term1(:,:,2) = df1s(p.numJ+1:2*p.numJ,p.numJ+1:2*p.numJ)*dMdx2_1;
%     cross_term1(:,:,3) = df1s(p.numJ+1:2*p.numJ,p.numJ+1:2*p.numJ)*dMdx3_1;
%     cross_term1(:,:,4) = df1s(p.numJ+1:2*p.numJ,p.numJ+1:2*p.numJ)*dMdx4_1;
%     cross_term1(:,:,5) = df1s(p.numJ+1:2*p.numJ,p.numJ+1:2*p.numJ)*dMdx5_1;
%     cross_term1(:,:,6) = df1s(p.numJ+1:2*p.numJ,p.numJ+1:2*p.numJ)*dMdx6_1;
%     
%     cross_term2(2,:,:) = df2s(p.numJ+1:2*p.numJ,p.numJ+1:2*p.numJ)*dMdx2_2;
%     cross_term2(3,:,:) = df2s(p.numJ+1:2*p.numJ,p.numJ+1:2*p.numJ)*dMdx3_2;
%     cross_term2(4,:,:) = df2s(p.numJ+1:2*p.numJ,p.numJ+1:2*p.numJ)*dMdx4_2;
%     cross_term2(5,:,:) = df2s(p.numJ+1:2*p.numJ,p.numJ+1:2*p.numJ)*dMdx5_2;
%     cross_term2(6,:,:) = df2s(p.numJ+1:2*p.numJ,p.numJ+1:2*p.numJ)*dMdx6_2;
    cross_term1 = permute(cross_term1,[2,1,3]);
    cross_term2 = permute(cross_term2,[2,1,3]);
    df_dss1(1:p.numJ,1:p.numJ,:) = df_dss1(1:p.numJ,1:p.numJ,:)+cross_term1;  %in the equations these are positive, but later I will subtract them
    df_dss2(1:p.numJ,1:p.numJ,:) = df_dss1(1:p.numJ,1:p.numJ,:)+cross_term2;
    
    %% solve the cross term for df_dux  dM/dx(df/du)^T
    cross_term1= zeros(p.numJ,p.numJ,p.numJ);
    cross_term2= zeros(p.numJ,p.numJ,p.numJ);
    
    cross_term1(:,:,2) = dMdx2_1*df1u(:,p.numJ+1:2*p.numJ).';
    cross_term1(:,:,3) = dMdx3_1*df1u(:,p.numJ+1:2*p.numJ).';
    cross_term1(:,:,4) = dMdx4_1*df1u(:,p.numJ+1:2*p.numJ).';
    cross_term1(:,:,5) = dMdx5_1*df1u(:,p.numJ+1:2*p.numJ).';
    cross_term1(:,:,6) = dMdx6_1*df1u(:,p.numJ+1:2*p.numJ).';
    
    cross_term2(:,:,2) = dMdx2_2*df2u(:,p.numJ+1:2*p.numJ).';
    cross_term2(:,:,3) = dMdx3_2*df2u(:,p.numJ+1:2*p.numJ).';
    cross_term2(:,:,4) = dMdx4_2*df2u(:,p.numJ+1:2*p.numJ).';
    cross_term2(:,:,5) = dMdx5_2*df2u(:,p.numJ+1:2*p.numJ).';
    cross_term2(:,:,6) = dMdx6_2*df2u(:,p.numJ+1:2*p.numJ).';
    
    cross_term1 = permute(cross_term1,[3,1,2]);
    cross_term2 = permute(cross_term2,[3,1,2]);
    
    df_dss1(2*p.numJ+1:3*p.numJ,1:p.numJ,:)=cross_term1;
    df_dss2(2*p.numJ+1:3*p.numJ,1:p.numJ,:)=cross_term2;
    % remember, if there is df/dux, there is also df/dxu
    df_dss1(1:p.numJ,2*p.numJ+1:3*p.numJ,:)=permute(cross_term1,[2,1,3]);
    df_dss2(1:p.numJ,2*p.numJ+1:3*p.numJ,:)=permute(cross_term2,[2,1,3]);
    %% solve the cross term for df_dFx  dM/dx(df/dF)^T
    cross_term1= zeros(p.numJ,p.numJ,4);
    cross_term2= zeros(p.numJ,p.numJ,4);
    
    cross_term1(2,:,:) = dMdx2_1*df1f(:,p.numJ+1:2*p.numJ).';
    cross_term1(3,:,:) = dMdx3_1*df1f(:,p.numJ+1:2*p.numJ).';
    cross_term1(4,:,:) = dMdx4_1*df1f(:,p.numJ+1:2*p.numJ).';
    cross_term1(5,:,:) = dMdx5_1*df1f(:,p.numJ+1:2*p.numJ).';
    cross_term1(6,:,:) = dMdx6_1*df1f(:,p.numJ+1:2*p.numJ).';
    
    cross_term2(2,:,:) = dMdx2_2*df2f(:,p.numJ+1:2*p.numJ).';
    cross_term2(3,:,:) = dMdx3_2*df2f(:,p.numJ+1:2*p.numJ).';
    cross_term2(4,:,:) = dMdx4_2*df2f(:,p.numJ+1:2*p.numJ).';
    cross_term2(5,:,:) = dMdx5_2*df2f(:,p.numJ+1:2*p.numJ).';
    cross_term2(6,:,:) = dMdx6_2*df2f(:,p.numJ+1:2*p.numJ).';
    
    cross_term1 = permute(cross_term1,[3,1,2]);
    cross_term2 = permute(cross_term2,[3,1,2]);
    df_dss1(3*p.numJ+1:3*p.numJ+4,1:p.numJ,:) = cross_term1;
    df_dss2(3*p.numJ+1:3*p.numJ+4,1:p.numJ,:) = cross_term2;
    
    df_dss1(1:p.numJ,3*p.numJ+1:3*p.numJ+4,:)=permute(cross_term1,[2,1,3]);
    df_dss2(1:p.numJ,3*p.numJ+1:3*p.numJ+4,:)=permute(cross_term2,[2,1,3]);
    
    
    df1 = Mddq_d1 - df_dss1;
    df2 = Mddq_d2 - df_dss2;
    
    
    
    
    %the last step is divide by M
    M1 = six_M(x1(2),x1(3),x1(4),x1(5),x1(6));
    M2 = six_M(x2(2),x2(3),x2(4),x2(5),x2(6));
    for it1=1:2*p.numJ
        for it2=1:2*p.numJ
            df1(it1,it2,:) = reshape(df1(it1,it2,:),[1,p.numJ])/M1;
            df2(it1,it2,:) = reshape(df2(it1,it2,:),[1,p.numJ])/M2;
        end
        
    end
    
    %in the original f function, f is (2*p.numJ), we omit it
    %before since the first p.numJ ones are zeros, now we add them back
    df1_new = zeros(3*p.numJ+4,3*p.numJ+4,2*p.numJ);
    df2_new = zeros(3*p.numJ+4,3*p.numJ+4,2*p.numJ);
    df1_new(:,:,p.numJ+1:2*p.numJ) = df1;
    df2_new(:,:,p.numJ+1:2*p.numJ)=df2;
    
    % now we have df_dxx terms (22x22), yet, it is not the hessian of the
    % constraints
    
% %     c1 = -p.sampT/6*(df1+
%     cross_term1 = zeros(3*p.numJ+4,3*p.numJ+4,2*p.numJ);
%     cross_term2 = zeros(3*p.numJ+4,3*p.numJ+4,2*p.numJ);
%     
%     for it1=1:2*p.numJ
%         for it2=1:2*p.numJ
%             cross_term1(it1,it2,:) = reshape(df1_new(it1,it2,:),[1,2*p.numJ])*df_half_s;
%             cross_term2(it1,it2,:) = reshape(df2_new(it1,it2,:),[1,2*p.numJ])*df_half_s;
%         end
%     end
%     
%     cross_term1_1 = zeros(3*p.numJ+4,3*p.numJ+4,2*p.numJ);
%     cross_term2_1 = zeros(3*p.numJ+4,3*p.numJ+4,2*p.numJ);
%     for it1=1:2*p.numJ
%         cross_term1_1(it1,:,:) = reshape(df1_new(it1,:,:),[3*p.numJ+4,2*p.numJ])*df_half_s;
%         cross_term2_1(it1,:,:) = reshape(df2_new(it1,:,:),[3*p.numJ+4,2*p.numJ])*df_half_s;
%     end
%     
%     
%     dcdxx1 = -p.sampT/6*(df1_new+p.sampT/8*cross_term1);
%     dcdxx2 = -p.sampT/6*(df2_new-p.sampT/8*cross_term2);
    

    %here I use normal traptodiz method
    
    dcdxx1 = -p.sampT/2*df1_new;
    dcdxx2 = -p.sampT/2*df2_new;


    
    
    for i_con = 1:p.numJ*2 %each time step has numJ constraints
        
        curHout = zeros(size(x,1));
        curHout((i-1)*(3*p.numJ+4)+1:i*(3*p.numJ+4),(i-1)*(3*p.numJ+4)+1:i*(3*p.numJ+4))=dcdxx1(:,:,i_con);
        curHout(i*(3*p.numJ+4)+1:(i+1)*(3*p.numJ+4),i*(3*p.numJ+4)+1:(i+1)*(3*p.numJ+4))=dcdxx2(:,:,i_con);
        Hout = Hout+lambda.eqnonlin(eq_idx)*curHout;
        eq_idx = eq_idx+1;
        
    end
    
    
    
    
    
end

% grf eqn constraints

for i=1:floor(p.gaitT/p.sampT)+1
    %toe constraints
    curHout = zeros(size(x,1));
    curX = x((i-1)*(p.numJ*3+4)+1:i*(p.numJ*3+4));
    cur_h = hess_grf_ceq_toe(curX.',p.toe_th,p.dmax,p.cmax,p.k,p.us,p.ud);
    curHout((i-1)*(p.numJ*3+4)+1:i*(p.numJ*3+4),(i-1)*(p.numJ*3+4)+1:i*(p.numJ*3+4))=cur_h;
    Hout = Hout+lambda.eqnonlin(eq_idx)*curHout;
    eq_idx = eq_idx+1; 
    %heel constraints
    curHout = zeros(size(x,1));
    cur_h = hess_grf_ceq_heel(curX.',p.toe_th,p.dmax,p.cmax,p.k,p.us,p.ud);
    curHout((i-1)*(p.numJ*3+4)+1:i*(p.numJ*3+4),(i-1)*(p.numJ*3+4)+1:i*(p.numJ*3+4))=cur_h;
    Hout = Hout+lambda.eqnonlin(eq_idx)*curHout;
    eq_idx = eq_idx+1; 
end


%hip_len constraints
curHout = zeros(size(x,1));
x1 = x(1:p.numJ*3+4);
x2 = x(end-p.numJ*3-3:end);
cur_h = h_hipLen(x1.',x2.',p.hipLen);

curHout(1:3*p.numJ+4,1:3*p.numJ+4) = cur_h(1:3*p.numJ+4,1:3*p.numJ+4);
curHout(end-3*p.numJ-3:end,end-3*p.numJ-3:end) = cur_h(3*p.numJ+5:end,3*p.numJ+5:end);
Hout = Hout+lambda.eqnonlin(eq_idx)*curHout;
    
   
    







end