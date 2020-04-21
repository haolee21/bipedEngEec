function [c,gradc]=headPosCon(x,p)

c = zeros(1,size(x,2));
% c = gpuArray(c);
gradc = zeros(size(x,1),size(x,2),size(x,2));
% gradc = gpuArray(gradc);
for i=1:size(x,2)
    curX = x(:,i);
    c(i)=p.head_h-five_link_head_h(curX.');
    gradc(:,i,i) = [-1*five_link_head_grad(curX),zeros(1,12)].';
end
c = c.';
gradc = reshape(gradc,[size(x,1)*size(x,2),size(x,2)]);

end