function [u,J_toe,J_heel,c_toe,c_heel] = u_no_ext(q1,q2,q3,p)

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





tend_ank1 = (tend_ank1_1+2*tend_ank1_2+tend_ank1_3).'/4*p.sampT;
tend_ank2 = (tend_ank2_1+2*tend_ank2_2+tend_ank2_3).'/4*p.sampT;

tend_kne1 = (tend_kne1_1+2*tend_kne1_2+tend_kne1_3).'/4*p.sampT;
tend_kne2 = (tend_kne2_1+2*tend_kne2_2+tend_kne2_3).'/4*p.sampT;

dq =(q3-q1)/2;
fri_tau = p.joint_fri*eye(p.numJ)*dq; %we simplify many terms here, dq1 = (q2-q1)/sampT,  dq2 = (q3-q2)/sampT, thus, the average is (q3-q1)/2/sampT

u=-dL1-dL2-tend_ank1-tend_ank2-tend_kne1-tend_kne2+fri_tau;


J_toe1=six_J(q1(1),q1(2),q1(3),q1(4),q1(5),q1(6));
J_toe2=six_J(q2(1),q2(2),q2(3),q2(4),q2(5),q2(6));
J_toe3=six_J(q3(1),q3(2),q3(3),q3(4),q3(5),q3(6));
J_toe=(J_toe1+2*J_toe2+J_toe3)/4;
J_toe = J_toe(1:2,:);

J_heel1=six_J2(q1(1),q1(2),q1(3),q1(4),q1(5),q1(6));
J_heel2=six_J2(q2(1),q2(2),q2(3),q2(4),q2(5),q2(6));
J_heel3=six_J2(q3(1),q3(2),q3(3),q3(4),q3(5),q3(6));
J_heel=(J_heel1+2*J_heel2+J_heel3)/4;
J_heel = J_heel(1:2,:);

toe_vel = J_toe*dq;
heel_vel = J_heel*dq;

toe_h1=end_y_pos(q1.');
toe_h2=end_y_pos(q2.');
toe_h3=end_y_pos(q3.');
heel_h1=heel_y_pos(q1.');
heel_h2=heel_y_pos(q2.');
heel_h3=heel_y_pos(q3.');

toe_h = (toe_h1+2*toe_h2+toe_h3)/4;
heel_h=(heel_h1+2*heel_h2+heel_h3)/4;

c_toe = p.k0/p.k1;
c_heel=p.k0/p.k1;

if(toe_h<p.toe_th)
    cur_c=p.k*(p.toe_th-toe_h)^2-p.cmax*(p.toe_th-toe_h)/p.dmax*toe_vel(2);
    if(cur_c>0)
        c_toe = p.k0/(cur_c^2*0.1+p.k1);
    end
    
end
if(heel_h<p.toe_th)
    cur_c = p.k*(p.toe_th-heel_h)^2-p.cmax*(p.toe_th-heel_h)/p.dmax*heel_vel(2);
    if(cur_c>0)
        c_heel=p.k0/(cur_c^2*0.1+p.k1);
    end
end
end