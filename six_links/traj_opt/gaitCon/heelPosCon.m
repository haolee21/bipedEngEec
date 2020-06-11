function [c,gradc] = heelPosCon(x,p)

start_i = floor(size(x,2)*0.4);
end_i = floor(size(x,2)*0.6);
c = -ones(1,end_i-start_i+1);
% c = gpuArray(c);
gradc = zeros(size(x,1),size(x,2),size(c,2));
% gradc = gpuArray(gradc);

%% we made some modification, originally we only need ypos>0, in order to create ground clearance, we need in the middle > param.toe_th
% 30% gait -> 70% gait, higher than 2*p.toe_th since it is a 0.5* tanh+0.5 function,
% tanh(0) = 0.5


for i=1:end_i-start_i+1
    curX = x(:,start_i+i-1);
    c(i)=p.gndclear-heel_y_pos(curX.');
    gradc(:,i+start_i-1,i) = [-1*heel_y_grad(curX.'),zeros(1,p.numJ*2)].';
end
c = c.';
gradc = reshape(gradc,[size(x,1)*size(x,2),size(c,1)]); %c was tranpose in the previous line



end