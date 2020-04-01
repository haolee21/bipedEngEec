function tou = newtonEuler(rMat,p,pc,m,Ic,var,var_dot,var_dot2,total_var)
for i=1:length(total_var) 
    syms(char(total_var(i))); 
end
w=sym(zeros(3,length(var_dot)+1));
w_dot=sym(zeros(3,length(var_dot)+1));
w_dot2=sym(zeros(3,length(var_dot)+1));
v_dot=sym(zeros(3,length(var_dot)+1));
vc_dot=sym(zeros(3,length(var)));
F=sym(zeros(3,length(var_dot)));
f=sym(zeros(3,length(var_dot)+1));
N=sym(zeros(3,length(var_dot)+1));
n=sym(zeros(3,length(var_dot)+1));
%initialize v_dot(1) = [0,g,0]';
v_dot(:,1)=vpa([-g 0 0].');

f(:,1)=[fx;fy;fz];
for i=1:length(rMat)-1
    rotM_temp=rMat{1,i};
    w(:,i+1)=rotM_temp.'*w(:,i)+[0 0 var_dot(i)].';
    w_dot(:,i+1)=rotM_temp.'*w_dot(:,i)+rotM_temp.'*cross(w(:,i),[0 0 var_dot(i)].')+[0 0 var_dot2(i)].';
    v_dot(:,i+1)=rotM_temp.'*(cross(w_dot(:,i),p(:,i))+cross(w(:,i),cross(w(:,i),p(:,i)))+v_dot(:,i));
    vc_dot(:,i+1)=cross(w_dot(:,i+1),pc(:,i))+cross(w(:,i+1),cross(w(:,i+1),pc(:,i)))+v_dot(:,i+1);
    F(:,i+1)=m(i)*vc_dot(:,i+1);
    N(:,i+1)=Ic{1,i}*w_dot(:,i+1)+cross(w(:,i+1),Ic{1,i}*w(:,i+1));
end
%now solve inward iteration
rMat=fliplr(rMat);
N=fliplr(N);
F=fliplr(F);
pc=fliplr(pc);
for i=1:length(rMat)-1
    rotM_temp=rMat{1,i};
    f(:,i+1)=rotM_temp*f(:,i)+F(:,i);
    n(:,i+1)=N(:,i)+rotM_temp*n(:,i)+cross(pc(:,i),F(:,i))+cross(p(:,i),rotM_temp*f(:,i));
end
n=fliplr(n);
for i=1:length(rMat)-1
    tou(1,i)=[0 0 1]*n(:,i);
end
tou=vpa(tou,3);