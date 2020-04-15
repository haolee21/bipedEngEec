function [c,gradc] = yposCon(x)

c = gpuArray(zeros(1,size(x,2)));
gradc = gpuArray(zeros(size(x,1),size(x,2),size(x,2)));
for i=1:size(x,2)
    curX = x(:,i);
    c(i)=-1*end_y_pos(curX.');
    gradc(:,i,i) = [end_y_poss_grad(curX),zeros(1,10)].';
end
c = c.';
gradc = reshape(gradc,[size(x,1)*size(x,2),size(x,2)]);


end