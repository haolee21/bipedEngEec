function [dx,dxGrad] = dynamics(x,tauext,jointNum)

q = x(1:jointNum,:);
dq = x(jointNum+1:2*jointNum,:);
u = x



end