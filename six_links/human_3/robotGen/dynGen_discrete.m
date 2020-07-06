function dyn = dynGen_discrete(robot,end_eff,sampT)
%%
numJ = length(robot.links);
q1 = sym('q1',[1,numJ]);
q2 = sym('q2',[1,numJ]);
% syms q1 [1 numJ] 
% syms q2 [1 numJ]
assume(q1,'real');
assume(q2,'real');
q_t = 0.5*(q1+q2); %turn q into symbolic expression
qd_t = (q2-q1)/sampT;

% symbols for output

% syms G_sym  % this is for easier to find G matrix, the robot should already have the g value

% calculate total L, K

curRot = eye(4);
% g = robot.gravity*G_sym;
g = robot.gravity;
V_cur = sym(zeros(1,numJ));
K_mat = sym(zeros(1,numJ));

% calculate joint base velocity/rotational matrix  (calculate first in
% order to use parallel for later)
j_vel = sym(zeros(3,numJ));
rotM = cell(1,numJ);
w_all = cell(1,numJ);
cur_w = [0;0;0];
for i=2:numJ % base has no translational velocity
    cur_q = sym(zeros(1,i-1));
    for k=1:i-1
        cur_q(1,k) = q_t(k);
    end
    curRot =curRot*turnRTtoMatrix(robot.A(i-1,cur_q));
    rotM{1,i-1} = curRot;
    
    
    if robot.links(i).jointtype=='R'
        cur_w = curRot(1:3,1:3)*[0;0;qd_t(i-1)]+cur_w;
        j_vel(:,i) = j_vel(:,i-1)+ simplify(curRot(1:3,1:3)*cross(cur_w,[robot.links(i).a;0;0]));
        w_all{1,i-1}=cur_w;
    end
end
% in the above loop we have not calculate the last rotational matrix, which
% will be done here
cur_q = sym(zeros(1,numJ));
for k=1:numJ
    cur_q(1,k) = q_t(k);
end
curRot = curRot*turnRTtoMatrix(robot.A(numJ,cur_q));
rotM{1,numJ} = curRot;
w_all{1,numJ} = cur_w + curRot(1:3,1:3)*[0;0;qd_t(numJ)];



%% if we only need new jacobian, we can ignore this cell
% 
dL1 = sym(zeros(numJ,numJ));
dL2 = sym(zeros(numJ,numJ));


% for some reason in parfor I cannot called diff for symbolic, here I
% self-defined one
% syms x
% diff_sym = matlabFunction(diff(f1(x),x));

for i=1:numJ
%     potential energy
    lc = [robot.links(i).r,1].';
    
    lc_base = rotM{1,i}*lc;
    V_cur = robot.links(i).m*g.'*lc_base(1:3,1);
    
    
%     kinetic energy
    if robot.links(i).jointtype == 'R' % if this is revolute joint
%         rotational energy
        cur_w = w_all{1,i};
        T_cur = simplify(0.5*cur_w.'*(robot.links(i).I)*cur_w);
%         translational energy
        cur_v = j_vel(:,i)+rotM{1,i}(1:3,1:3)*cross(cur_w,lc(1:3,1));
        T_cur = T_cur+0.5*robot.links(i).m*(cur_v.'*cur_v);
        
        L = T_cur-V_cur;
        
    end

%     % for discrete Lagrangian, we calculate the derivate with respect to q1, q2
%     % for L = L(q1,q2) 
    for c=1:i %there is no point to calculate joint more than i since it will be zero (e.g., link 2's energy is not related to joint 5)
% %         
        dL1(c,i) = diff(L,q1(1,c));
        dL2(c,i) = diff(L,q2(1,c));

    end
        
     

    
    
    
end
% add up L1, L2 terms
dL1 = sum(dL1,2);
dL2 = sum(dL2,2);
% % build M,V,G matrix

% totL = L1-L2;
% dyn.M = sym(zeros(numJ));
% dyn.G = sym(zeros(1,numJ));
% dyn.V = sym(zeros(1,numJ));
% for c=1:numJ
%     for i =1:numJ
%         dyn.M(c,i) = simplify(diff(totL(c,1),dqd(1,i)));
%         
%     end
%     dyn.G(1,c) = diff(totL(c,1),G_sym);
%     dyn.V(1,c) = simplify(totL(c,1)-dyn.M(c,:)*dqd.'-dyn.G(1,c)*G_sym);
%     
%     
% end

dyn.dL1 =dL1;
dyn.dL2 = dL2;


%%
% to make the calculation simplier, we did not add external force, thus
% ,not jacobian can be found
% for here, we directly use forward kinematic to find the jacobian

% yet, jacobian only focus on position at single time frame, for here we
% use q1


% Jacobian at the hip, this is for making sure it won't go backward
Hip = turnRTtoMatrix(robot.A(1:1:3,q1))*[0;0;0;1];
hipPos = Hip(1:3,1);

dyn.J_hip = sym(zeros(3,numJ));
for i=1:numJ
    for k=1:3
        dyn.J_hip(k,i) = simplify(diff(hipPos(k,1),q1(i)));
    end
    
end
% 
% 
% % Jacobian at the toe
endT  = turnRTtoMatrix(robot.A(1:1:numJ,q1))*[end_eff(1);end_eff(2);0;1]; % end_eff is [l_foot,l_heel,l_calf], for the toe, it is [l_foot,l_heel,0] position of ankle joint
endPos = endT(1:3,1);

dyn.J = sym(zeros(3,numJ));
for i=1:numJ
    for k=1:3
        dyn.J(k,i) = simplify(diff(endPos(k,1),q1(i)));
        
    end
    
end
% %% this is for a special jacobian on heel, since external force also act on the heel
endT2 = turnRTtoMatrix(robot.A(1:1:numJ,q1))*[0;end_eff(3);0;1];  % heel is [0,l_heel,0] position of ankle joint
endPos2 = endT2(1:3,1);
dyn.J2 = sym(zeros(3,numJ));
% end_w2 = subs(w_all{1,end},[q_t,qd_t],[q,qd]);
for i=1:numJ
    for k=1:3
        dyn.J2(k,i) = simplify(diff(endPos2(k,1),q1(i)));
%         dyn.J2(k+3,i) = simplify(diff(end_w2(k,1),qd(i))); 
%         
    end
%     
end
end