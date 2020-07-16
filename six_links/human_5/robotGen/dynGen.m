function dyn = dynGen(robot,end_eff)
%%
numJ = length(robot.links);
syms qt(t) [1 numJ];
syms qdt(t) [1 numJ];
q_t = qt(t); %turn q into symbolic expression
qd_t = qdt(t);
assume(q_t,'real');
assume(qd_t,'real');
% symbols for output
syms q [1 numJ] 
syms qd [1 numJ] 
syms dqd [1 numJ]
syms G_sym  % this is for easier to find G matrix, the robot should already have the g value

% calculate total L, K

curRot = eye(4);
g = robot.gravity*G_sym;
T_mat = sym(zeros(1,numJ));
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
L1 = sym(zeros(numJ,numJ));
L2 = sym(zeros(numJ,numJ));


for i=1:numJ
%     potential energy
    lc = [robot.links(i).r,1].';
    
    lc_base = rotM{1,i}*lc;
    T_mat(1,i) = simplify(robot.links(i).m*g.'*lc_base(1:3,1));
    
    
%     kinetic energy
    if robot.links(i).jointtype == 'R' % if this is revolute joint
%         rotational energy
        cur_w = w_all{1,i};
        K_mat(1,i) = K_mat(1,i)+simplify(0.5*cur_w.'*(robot.links(i).I)*cur_w);
%         translational energy
        cur_v = simplify(j_vel(:,i)+rotM{1,i}(1:3,1:3)*cross(cur_w,lc(1:3,1)));
        K_mat(1,i) = simplify(K_mat(1,i)+0.5*robot.links(i).m*(cur_v.'*cur_v));
        
        L = K_mat(1,i)-T_mat(1,i);
%         calculate d/dt(d/d_qdot(L))
        for c=1:i 
%             we will have to substitute qd_t with qd first to calculate
%             d/qd_dot since qd_t is a symfun, not sym
            L1_row = subs(L,qd_t(1,c),qd(1,c));
            L1_row = diff(L1_row,qd(1,c));
            L1_row = subs(L1_row,qd(1,c),qd_t(1,c));
            L1_row = diff(L1_row,t);
            
            L2_row = subs(L,q_t(1,c),q(1,c));
            L2_row = diff(L2_row,q(1,c));
            L1(c,i) = L1_row;
            L2(c,i) = L2_row;
        end

        for c=1:i
            for k=1:size(L1,1)
                L1(k,i) = subs(L1(k,i),diff(q_t(1,c),t),qd(1,c));
                L1(k,i) = subs(L1(k,i),diff(qd_t(1,c),t),dqd(1,c));    
                L1(k,i) = subs(L1(k,i),q_t(1,c),q(1,c));
                L1(k,i) = subs(L1(k,i),qd_t(1,c),qd(1,c));
                L1(k,i) = simplify(L1(k,i));
            
                L2(k,i) = subs(L2(k,i),diff(q_t(1,c),t),qd(1,c));
                L2(k,i) = subs(L2(k,i),diff(qd_t(1,c),t),dqd(1,c));    
                L2(k,i) = subs(L2(k,i),q_t(1,c),q(1,c));
                L2(k,i) = subs(L2(k,i),qd_t(1,c),qd(1,c));
                L2(k,i) = simplify(L2(k,i));
            end
        end

    end
    
    
end
% add up L1, L2 terms
% 
% 
% % build M,V,G matrix
L1 = sum(L1,2);
L2 = sum(L2,2);
totL = L1-L2;
dyn.M = sym(zeros(numJ));
dyn.G = sym(zeros(1,numJ));
dyn.V = sym(zeros(1,numJ));
for c=1:numJ
    for i =1:numJ
        dyn.M(c,i) = simplify(diff(totL(c,1),dqd(1,i)));
        
    end
    dyn.G(1,c) = diff(totL(c,1),G_sym);
    dyn.V(1,c) = simplify(totL(c,1)-dyn.M(c,:)*dqd.'-dyn.G(1,c)*G_sym);
    
    
end




%%
% to make the calculation simplier, we did not add external force, thus
% ,not jacobian can be found
% for here, we directly use forward kinematic to find the jacobian

% Jacobian at the hip, this is for making sure it won't go backward
Hip = turnRTtoMatrix(robot.A(1:1:3,q))*[0;0;0;1];
hipPos = Hip(1:3,1);
hip_w = subs(w_all{1,3},[q_t,qd_t],[q,qd]);
dyn.J_hip = sym(zeros(6,numJ));
for i=1:numJ
    for k=1:3
        dyn.J_hip(k,i) = simplify(diff(hipPos(k,1),q(i)));
        dyn.J_hip(k+3,i) = simplify(diff(hip_w(k,1),qd(i))); 
    end
    
end


% Jacobian at the toe
endT  = turnRTtoMatrix(robot.A(1:1:numJ,q))*[end_eff(1);end_eff(2);0;1]; % end_eff is [l_foot,l_heel,l_calf], for the toe, it is [l_foot,l_heel,0] position of ankle joint
endPos = endT(1:3,1);
end_w = subs(w_all{1,end},[q_t,qd_t],[q,qd]);
dyn.J = sym(zeros(6,numJ));
for i=1:numJ
    for k=1:3
        dyn.J(k,i) = simplify(diff(endPos(k,1),q(i)));
        dyn.J(k+3,i) = simplify(diff(end_w(k,1),qd(i))); 
    end
    
end
%% this is for a special jacobian on heel, since external force also act on the heel
endT2 = turnRTtoMatrix(robot.A(1:1:numJ,q))*[0;end_eff(3);0;1];  % heel is [0,l_heel,0] position of ankle joint
endPos2 = endT2(1:3,1);
dyn.J2 = sym(zeros(6,numJ));
end_w2 = subs(w_all{1,end},[q_t,qd_t],[q,qd]);
for i=1:numJ
    for k=1:3
        dyn.J2(k,i) = simplify(diff(endPos2(k,1),q(i)));
        dyn.J2(k+3,i) = simplify(diff(end_w2(k,1),qd(i))); 
        
    end
    
end
end