function [c,gradc]=headPosCon(x,p)

c = zeros(1,size(x,2));
% c = gpuArray(c);
gradc = zeros(size(x,1),size(x,2),size(x,2));
% gradc = gpuArray(gradc);
for i=1:size(x,2)
    curX = x(:,i);
    c(i)=p.head_h-head_y_pos(curX.');
    gradc(:,i,i) = [-1*head_y_grad(curX.'),zeros(1,p.numJ*3-3)].';
end
c = c.';
gradc = reshape(gradc,[size(x,1)*size(x,2),size(x,2)]);

end