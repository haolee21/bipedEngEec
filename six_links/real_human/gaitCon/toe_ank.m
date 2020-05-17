function [c,gradc] = toe_ank(x,p)

c = -ones(1,size(x,2));
gradc = zeros(size(x,1),size(x,2),size(x,2));

start_i = floor(size(x,2)*0.2);
end_i = floor(size(x,2)*0.8);

%% toe should be lower than heel at the beginning
for i =1:start_i
    curX = x(:,i);
    c(i) =end_y_pos(curX.')- heel_y_pos(curX.');
    gradc(:,i,i)=[end_y_grad(curX.')-heel_y_grad(curX.'),zeros(1,p.numJ*2)].';
end

%% ankle should be lower than heel at the end
for i=end_i:size(x,2)
    curX=x(:,i);
    c(i) = heel_y_pos(curX.')-end_y_pos(curX.');
    gradc(:,i,i)=[heel_y_grad(curX.')-end_y_grad(curX.'),zeros(1,p.numJ*2)].';
end
c = c.';
gradc = reshape(gradc,[size(x,1)*size(x,2),size(x,2)]);

end