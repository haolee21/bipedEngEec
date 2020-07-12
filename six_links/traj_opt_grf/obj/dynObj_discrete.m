function [dObj,dObjGrad]=dynObj_discrete(x,p)

dObj=0;

dObjGrad = zeros(size(x,1),size(x,2));
for i=1:size(x,2)-2
    q1 = x(1:p.numJ,i);
    q2 = x(1:p.numJ,i+1);
    q3 = x(1:p.numJ,i+2);
    
    u1 = x(p.numJ+1:p.numJ*2,i);
    u2 = x(p.numJ+1:p.numJ*2,i+1);
    u3 = x(p.numJ+1:p.numJ*2,i+2);
    
    f_toe1 = x(p.numJ*2+1:p.numJ*2+2,i);
    f_toe2 = x(p.numJ*2+1:p.numJ*2+2,i+1);
    f_toe3 = x(p.numJ*2+1:p.numJ*2+2,i+2);
    
    f_heel1 = x(p.numJ*2+3:p.numJ*2+4,i);
    f_heel2 = x(p.numJ*2+3:p.numJ*2+4,i+1);
    f_heel3 = x(p.numJ*2+3:p.numJ*2+4,i+2);
    
    s_toe1 = x(p.numJ*2+5,i);
    s_toe2 = x(p.numJ*2+5,i+1);
    s_toe3 = x(p.numJ*2+5,i+2);
    
    s_heel1= x(p.numJ*2+6,i);
    s_heel2= x(p.numJ*2+6,i+1);
    s_heel3= x(p.numJ*2+6,i+2);
    
    tend_ank1_1 = [pi/2-q1(1,1),0,0,0,0,0]*p.ank_stiff;
    tend_ank1_2 = [pi/2-q2(1,1),0,0,0,0,0]*p.ank_stiff;
    tend_ank1_3 = [pi/2-q3(1,1),0,0,0,0,0]*p.ank_stiff;
    
    tend_ank2_1 = [0,0,0,0,0,-pi/2-q1(6,1)]*p.ank_stiff;
    tend_ank2_2 = [0,0,0,0,0,-pi/2-q2(6,1)]*p.ank_stiff;
    tend_ank2_3 = [0,0,0,0,0,-pi/2-q3(6,1)]*p.ank_stiff;
    
    tend_kne1_1 = [0,-q1(2,1),0,0,0,0]*p.knee_stiff;
    tend_kne1_2 = [0,-q2(2,1),0,0,0,0]*p.knee_stiff;
    tend_kne1_3 = [0,-q3(2,1),0,0,0,0]*p.knee_stiff;
    
    tend_kne2_1 = [0,0,0,0,-q1(5,1),0]*p.knee_stiff;
    tend_kne2_2 = [0,0,0,0,-q2(5,1),0]*p.knee_stiff;
    tend_kne2_3 = [0,0,0,0,-q3(5,1),0]*p.knee_stiff;
    
    dL1 = [dL1_1((q2.'+q3.')/2,(q3.'-q2.')/p.sampT,p.sampT);...
        dL1_2((q2.'+q3.')/2,(q3.'-q2.')/p.sampT,p.sampT);...
        dL1_3((q2.'+q3.')/2,(q3.'-q2.')/p.sampT,p.sampT);...
        dL1_4((q2.'+q3.')/2,(q3.'-q2.')/p.sampT,p.sampT);...
        dL1_5((q2.'+q3.')/2,(q3.'-q2.')/p.sampT,p.sampT);...
        dL1_6((q2.'+q3.')/2,(q3.'-q2.')/p.sampT,p.sampT)]*p.sampT;
    
    dL2 = [dL2_1((q1.'+q2.')/2,(q2.'-q1.')/p.sampT,p.sampT);...
        dL2_2((q1.'+q2.')/2,(q2.'-q1.')/p.sampT,p.sampT);...
        dL2_3((q1.'+q2.')/2,(q2.'-q1.')/p.sampT,p.sampT);...
        dL2_4((q1.'+q2.')/2,(q2.'-q1.')/p.sampT,p.sampT);...
        dL2_5((q1.'+q2.')/2,(q2.'-q1.')/p.sampT,p.sampT);...
        dL2_6((q1.'+q2.')/2,(q2.'-q1.')/p.sampT,p.sampT)]*p.sampT;
    
    u_sum = (u1+2*u2+u3)*p.sampT/4;
    
    tau_toe_1 = Tau_toe(q1.',f_toe1.',s_toe1);
    tau_toe_2 = Tau_toe(q2.',f_toe2.',s_toe2);
    tau_toe_3 = Tau_toe(q3.',f_toe3.',s_toe3);
    
    tau_heel_1 = Tau_heel(q1.',f_heel1.',s_heel1);
    tau_heel_2 = Tau_heel(q2.',f_heel2.',s_heel2);
    tau_heel_3 = Tau_heel(q3.',f_heel3.',s_heel3);
    
    tau_toe = (tau_toe_1+2*tau_toe_2+tau_toe_3)*p.sampT/4;
    tau_heel = (tau_heel_1+2*tau_heel_2+tau_heel_3)*p.sampT/4;
    

    tend_ank1 = (tend_ank1_1+2*tend_ank1_2+tend_ank1_3).'/4*p.sampT;
    tend_ank2 = (tend_ank2_1+2*tend_ank2_2+tend_ank2_3).'/4*p.sampT;
    
    tend_kne1 = (tend_kne1_1+2*tend_kne1_2+tend_kne1_3).'/4*p.sampT;
    tend_kne2 = (tend_kne2_1+2*tend_kne2_2+tend_kne2_3).'/4*p.sampT;
    
    dq =(q3-q1)/2;
    fri_tau = p.joint_fri*eye(p.numJ)*dq; %we simplify many terms here, dq1 = (q2-q1)/sampT,  dq2 = (q3-q2)/sampT, thus, the average is (q3-q1)/2/sampT
    
    v = dL1+dL2+u_sum+tau_toe+tau_heel+tend_ank1+tend_ank2+tend_kne1+tend_kne2-fri_tau;
    
    dObj = dObj + v.'*v*0.5;
    
    %now we solve the gradient
    
    if nargout>1
        % the gradient is actually dvdx*v, dvdx is the gradient if you
        % treat this dynamic as constraint
        grad_temp = zeros(size(x,1),size(x,2),p.numJ);
        
        % gradient related to Lagrangian
        %q1 is related to dL2, first input
        grad_temp(1:p.numJ,i,:)=[dL21_dq1((q1.'+q2.')/2,(q2.'-q1.')/p.sampT,p.sampT),...
                                 dL22_dq1((q1.'+q2.')/2,(q2.'-q1.')/p.sampT,p.sampT),...
                                 dL23_dq1((q1.'+q2.')/2,(q2.'-q1.')/p.sampT,p.sampT),...
                                 dL24_dq1((q1.'+q2.')/2,(q2.'-q1.')/p.sampT,p.sampT),...
                                 dL25_dq1((q1.'+q2.')/2,(q2.'-q1.')/p.sampT,p.sampT),...
                                 dL26_dq1((q1.'+q2.')/2,(q2.'-q1.')/p.sampT,p.sampT)]*p.sampT;
        %q2 is related to dL2, first input, dL1, 2nd input
        grad_temp(1:p.numJ,i+1,:)= [dL21_dq2((q1.'+q2.')/2,(q2.'-q1.')/p.sampT,p.sampT)+dL11_dq1((q2.'+q3.')/2,(q3.'-q2.')/p.sampT,p.sampT),...
                                  dL22_dq2((q1.'+q2.')/2,(q2.'-q1.')/p.sampT,p.sampT)+dL12_dq1((q2.'+q3.')/2,(q3.'-q2.')/p.sampT,p.sampT),...
                                  dL23_dq2((q1.'+q2.')/2,(q2.'-q1.')/p.sampT,p.sampT)+dL13_dq1((q2.'+q3.')/2,(q3.'-q2.')/p.sampT,p.sampT),...
                                  dL24_dq2((q1.'+q2.')/2,(q2.'-q1.')/p.sampT,p.sampT)+dL14_dq1((q2.'+q3.')/2,(q3.'-q2.')/p.sampT,p.sampT),...
                                  dL25_dq2((q1.'+q2.')/2,(q2.'-q1.')/p.sampT,p.sampT)+dL15_dq1((q2.'+q3.')/2,(q3.'-q2.')/p.sampT,p.sampT),...
                                  dL26_dq2((q1.'+q2.')/2,(q2.'-q1.')/p.sampT,p.sampT)+dL16_dq1((q2.'+q3.')/2,(q3.'-q2.')/p.sampT,p.sampT)]*p.sampT;
        %q3 is related to dL2, 2nd input
        grad_temp(1:p.numJ,i+2,:)=[dL11_dq2((q2.'+q3.')/2,(q3.'-q2.')/p.sampT,p.sampT),...
                                   dL12_dq2((q2.'+q3.')/2,(q3.'-q2.')/p.sampT,p.sampT),...
                                   dL13_dq2((q2.'+q3.')/2,(q3.'-q2.')/p.sampT,p.sampT),...
                                   dL14_dq2((q2.'+q3.')/2,(q3.'-q2.')/p.sampT,p.sampT),...
                                   dL15_dq2((q2.'+q3.')/2,(q3.'-q2.')/p.sampT,p.sampT),...
                                   dL16_dq2((q2.'+q3.')/2,(q3.'-q2.')/p.sampT,p.sampT)]*p.sampT;
        
        % gradient related to tendons
        % font ankle
        grad_temp(1,i,1) = grad_temp(1,i,(i-1)*p.numJ+1)-p.ank_stiff/4*p.sampT;
        grad_temp(1,i+1,1) = grad_temp(1,i+1,(i-1)*p.numJ+1)-2*p.ank_stiff/4*p.sampT;
        grad_temp(1,i+2,1) = grad_temp(1,i+2,(i-1)*p.numJ+1)-p.ank_stiff/4*p.sampT;
        % back ankle
        grad_temp(6,i,6) = grad_temp(6,i,i*p.numJ)-p.ank_stiff/4*p.sampT;
        grad_temp(6,i+1,6) = grad_temp(6,i+1,i*p.numJ)-2*p.ank_stiff/4*p.sampT;
        grad_temp(6,i+2,6) = grad_temp(6,i+2,i*p.numJ)-p.ank_stiff/4*p.sampT;
        % front knee
        grad_temp(2,i,2) = grad_temp(2,i,(i-1)*p.numJ+2)-p.knee_stiff/4*p.sampT;
        grad_temp(2,i+1,2) = grad_temp(2,i+1,(i-1)*p.numJ+2)-2*p.knee_stiff/4*p.sampT;
        grad_temp(2,i+2,2) = grad_temp(2,i+2,(i-1)*p.numJ+2)-p.knee_stiff/4*p.sampT;
        % back knee
        grad_temp(5,i,5) =   grad_temp(5,i,(i-1)*p.numJ+5)-p.knee_stiff/4*p.sampT;
        grad_temp(5,i+1,5) = grad_temp(5,i+1,(i-1)*p.numJ+5)-2*p.knee_stiff/4*p.sampT;
        grad_temp(5,i+2,5) = grad_temp(5,i+2,(i-1)*p.numJ+5)-p.knee_stiff/4*p.sampT;
        
        
        % gradient related to u
        grad_temp(p.numJ+1:2*p.numJ,i,:) =eye(p.numJ)/4*p.sampT;
        grad_temp(p.numJ+1:2*p.numJ,i+1,:) = 2*eye(p.numJ)/4*p.sampT;
        grad_temp(p.numJ+1:2*p.numJ,i+2,:) = eye(p.numJ)/4*p.sampT;
        
        % gradient related to f_toe
        
        grad_temp(:,i,:) = reshape(grad_temp(:,i,:),[p.numJ*2+6,p.numJ])+dTau_toe(q1.',f_toe1.',s_toe1)/4*p.sampT;
        grad_temp(:,i+1,:) = reshape(grad_temp(:,i+1,:),[p.numJ*2+6,p.numJ])+2*dTau_toe(q2.',f_toe2.',s_toe2)/4*p.sampT;
        grad_temp(:,i+2,:) = reshape(grad_temp(:,i+2,:),[p.numJ*2+6,p.numJ])+dTau_toe(q3.',f_toe3.',s_toe3)/4*p.sampT;
        
        % gradient related to f_heel
        grad_temp(:,i,:) = reshape(grad_temp(:,i,:),[p.numJ*2+6,p.numJ])+dTau_heel(q1.',f_heel1.',s_heel1)/4*p.sampT;
        grad_temp(:,i+1,:) = reshape(grad_temp(:,i+1,:),[p.numJ*2+6,p.numJ])+2*dTau_heel(q2.',f_heel2.',s_heel2)/4*p.sampT;
        grad_temp(:,i+2,:) = reshape(grad_temp(:,i+2,:),[p.numJ*2+6,p.numJ])+dTau_heel(q3.',f_heel3.',s_heel3)/4*p.sampT;
        
        %gradient related to joint friction
        grad_temp(1:p.numJ,i,:)=reshape(grad_temp(1:p.numJ,i,:),[p.numJ,p.numJ])+eye(p.numJ)*p.joint_fri/2;
        grad_temp(1:p.numJ,i+2,:)=reshape(grad_temp(1:p.numJ,i+2,:),[p.numJ,p.numJ])-eye(p.numJ)*p.joint_fri/2;
        
        % multiply the vector (gradient of 0.5*v.'*v is dv/dx*v
        for idx=1:p.numJ
            dObjGrad(:,i) = dObjGrad(:,i)+reshape(grad_temp(:,i,:),[size(x,1),p.numJ])*v;
            dObjGrad(:,i+1)= dObjGrad(:,i+1)+reshape(grad_temp(:,i+1,:),[size(x,1),p.numJ])*v;
            dObjGrad(:,i+2) = dObjGrad(:,i+2)+reshape(grad_temp(:,i+2,:),[size(x,1),p.numJ])*v;
        end
        
    end
    
end




end