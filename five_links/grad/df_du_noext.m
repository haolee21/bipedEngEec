function grad = df_du_noext(q1,q2,q3,q4,q5,dq1,dq2,dq3,dq4,dq5,u)
% u is the joint torque with dim=5
% external force is excluded since we will calculate it later

M = five_M(q2,q3,q4,sq);
grad = M\eye(size(M,1));
end