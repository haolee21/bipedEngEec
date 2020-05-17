function [c,gradc] = yposCon(x,p)

c = -ones(1,size(x,2));
% c = gpuArray(c);
gradc = zeros(size(x,1),size(x,2),size(x,2));
% gradc = gpuArray(gradc);

%% we made some modification, originally we only need ypos>0, in order to create ground clearance, we need in the middle > param.toe_th
% 30% gait -> 70% gait, higher than 2*p.toe_th since it is a 0.5* tanh+0.5 function,
% tanh(0) = 0.5
start_i = floor(size(x,2)*0.3);
end_i = floor(size(x,2)*0.7);

for i=start_i:end_i
    curX = x(:,i);
    c(i)=p.gndclear-end_y_pos(curX.');
    gradc(:,i,i) = [-1*end_y_grad(curX.'),zeros(1,p.numJ*2)].';
end
c = c.';
gradc = reshape(gradc,[size(x,1)*size(x,2),size(x,2)]);



end