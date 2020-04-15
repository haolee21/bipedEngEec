function [dx,curTauExt] = dyn_dq(x,tauext,jointNum,toe_th)
% this function will calculate the dq terms 
% return a row with 2*jointNum elements

% x: variables that need to be solved by fmincon, [q,dq,u]
% tauext: torque causes by external 

% in order to use the function "rowfun" to make it efficient, the inputs
% are all row vectors

q = x(1:jointNum);
dq = x(jointNum+1:2*jointNum);
u = x(2*jointNum+1:3*jointNum);


M = five_M(q(2),q(3),q(4),q(5));
V = five_V(q(2),q(3),q(4),q(5),dq(1),dq(2),dq(3),dq(4),dq(5))*dq.';
G = five_G(q(1),q(2),q(3),q(4),q(5));

toe_pos = endPos(q(1),q(2),q(3),q(4),q(5));

% here we calculate the ddq and external torque
% to simplify the calculation, we calculate them at the same time

term1 = u.'-V-G;

if(toe_pos(2)<toe_th)
    curTauExt = (M\term1).';
else
    curTauExt = zeros(1,jointNum); 
end
ddq = M\(u.'-V-G.'+tauext.');

dx = [dq,ddq.'];




end