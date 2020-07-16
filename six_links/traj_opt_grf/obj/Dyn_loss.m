function loss = Dyn_loss(x,p)
loss = 0;
% calculate the forward dynamic
q = x(1:p.numJ,:);
dq = gradient(q)/p.sampT;
ddq = gradient(dq)/p.sampT;
c_toe = x(p.numJ+1,:);
c_heel = x(p.numJ+2,:);
tau = zeros(p.numJ,size(x,2));
for i =1:size(x,2)
    J_toe = six_J(q(1,i),q(2,i),q(3,i),q(4,i),q(5,i),q(6,i));
    J_heel = six_J2(q(1,i),q(2,i),q(3,i),q(4,i),q(5,i),q(6,i));
    J.toe{1,i} = J_toe(1:2,:);
    J.heel{1,i} = J_heel(1:2,:);
    tau(:,i) = u_no_ext([q(:,i);dq(:,i)],ddq(:,i),p);
end


tau_vec = reshape(tau,[size(tau,1)*size(tau,2),1]);

JT_toe_cell = cellfun(@transpose,J.toe,'UniformOutput',false);
JT_heel_cell = cellfun(@transpose,J.heel,'UniformOutput',false);
JT_toe_big = blkdiag(JT_toe_cell{:});
JT_heel_big = blkdiag(JT_heel_cell{:});
R = diag(repmat(p.jointW,1,size(x,2)));

W_toe_cell = arrayfun(p.c_matFun,c_toe,'UniformOutput',false);
W_toe = blkdiag(W_toe_cell{:});
W_heel_cell = arrayfun(p.c_matFun,c_heel,'UniformOutput',false);
W_heel = blkdiag(W_heel_cell{:});

cvx_solver sedumi
cvx_begin quiet
    variables f_toe(2*size(x,2),1) f_heel(2*size(x,2),1) u(p.numJ*size(x,2),1)
    minimize( norm(JT_toe_big*f_toe+JT_heel_big*f_heel+u-tau_vec)+u.'*R*u+f_toe.'*W_toe*f_toe+f_heel.'*W_heel*f_heel)
    subject to
        p.fri_mat*f_toe<=0;
        p.fri_mat*f_heel<=0;

%             M1*f_toe ==0;
%             M1*f_heel ==0;
%             M2*f_toe ==0;
%             M2*f_heel==0;
% %
%             u(1,1) == -u(param.numJ*size(x,2),1);
%             u(2,1) == -u(param.numJ*size(x,2)-1,1);
%             u(3,1) == -u(param.numJ*size(x,2)-2,1);
%             u(4,1) == -u(param.numJ*size(x,2)-3,1);
%             u(5,1) == -u(param.numJ*size(x,2)-4,1);
%             u(6,1) == -u(param.numJ*size(x,2)-5,1);
cvx_end
% reshape to fit the format later
f_toe = reshape(f_toe,[2,size(x,2)]);
f_heel = zeros(2,size(x,2));
%     f_heel= reshape(f_heel,[2,size(x,2)]);
u = reshape(u,[p.numJ,size(x,2)]);

for i=1:size(x,2)
    err = J.toe{1,i}.'*f_toe(:,i)+J.heel{1,i}.'*f_heel(:,i)+u(:,i)-tau(:,i);
    loss = loss +err.'*err*10^-3;
    
end


end