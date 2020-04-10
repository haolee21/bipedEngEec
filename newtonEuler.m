function dynEqn = newtonEuler(mod)
%function tou = newtonEuler(rMat,p,pc,m,Ic,var,var_dot,var_dot2,total_var)
for i=1:length(mod.total_var) 
    syms(char(mod.total_var(i))); 
    assume(mod.total_var(i),'real');
end
w=sym(zeros(3,length(mod.var_dot)+1));
w_dot=sym(zeros(3,length(mod.var_dot)+1));
w_dot2=sym(zeros(3,length(mod.var_dot)+1));
v_dot=sym(zeros(3,length(mod.var_dot)+1));
vc_dot=sym(zeros(3,length(mod.var)));
F=sym(zeros(3,length(mod.var_dot)));
f=sym(zeros(3,length(mod.var_dot)+1));
N=sym(zeros(3,length(mod.var_dot)+1));
n=sym(zeros(3,length(mod.var_dot)+1));
%initialize v_dot(1) = [0,g,0]';
g = mod.g; % I am just lazy here, did not want to change all the g terms 
v_dot(:,1)=vpa([0 -g 0].');

f(:,1)=[fx;fy;fz];
for i=1:length(mod.rMat)-1
    rotM_temp=mod.rMat{1,i};
    w(:,i+1)=rotM_temp.'*w(:,i)+[0 0 mod.var_dot(i)].';
    w_dot(:,i+1)=rotM_temp.'*w_dot(:,i)+rotM_temp.'*cross(w(:,i),[0 0 mod.var_dot(i)].')+[0 0 mod.var_dot2(i)].';
    v_dot(:,i+1)=rotM_temp.'*(cross(w_dot(:,i),mod.p(:,i))+cross(w(:,i),cross(w(:,i),mod.p(:,i)))+v_dot(:,i));
    vc_dot(:,i+1)=cross(w_dot(:,i+1),mod.pc(:,i))+cross(w(:,i+1),cross(w(:,i+1),mod.pc(:,i)))+v_dot(:,i+1);
    F(:,i+1)=mod.m(i)*vc_dot(:,i+1);
    N(:,i+1)=mod.Ic{1,i}*w_dot(:,i+1)+cross(w(:,i+1),mod.Ic{1,i}*w(:,i+1));
end
%now solve inward iteration
mod.rMat=fliplr(mod.rMat);
N=fliplr(N);
F=fliplr(F);
mod.pc=fliplr(mod.pc);
for i=1:length(mod.rMat)-1
    rotM_temp=mod.rMat{1,i};
    f(:,i+1)=rotM_temp*f(:,i)+F(:,i);
    n(:,i+1)=N(:,i)+rotM_temp*n(:,i)+cross(mod.pc(:,i),F(:,i))+cross(mod.p(:,i),rotM_temp*f(:,i));
end
n=fliplr(n);
for i=1:length(mod.rMat)-1
    tau(1,i)=[0 0 1]*n(:,i);
end

tau = vpa(tau,5);
dynEqn.tau = tau;


varLen = length(mod.var_dot);
M = sym(zeros(varLen));
G = sym(zeros(varLen,1));
F = sym(zeros(varLen,3));

for i=1:length(varLen)
    for j=1:varLen
        M(i,j) = simplify(diff(tau(1,i),mod.var_dot2(1,j)));
    end
    for j =1:length(mod.fext)
        F(i,j) = simplify(diff(tau(1,i),mod.fext(j)));
    end
end
G = simplify(diff(tau,g).');

dynEqn.M=M;

dynEqn.F = F;

V = tau.'-dynEqn.M*mod.var_dot2.'-G*g-F*[fx;fy;fz];
% for i =1:length(mod.var)
%     V(i,1) = simplify(V(i,1));
%     tau(1,i) = simplify(tau(1,i));
%     
% end
dynEqn.V = V;
dynEqn.G=G;
% dynEqn.tau = tau;
end