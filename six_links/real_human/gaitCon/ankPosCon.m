function [c,gradc]=ankPosCon(x,p)
c = -ones(1,size(x,2));
gradc = zeros(size(x,1),size(x,2),size(x,2));

start_i = floor(size(x,2)*0.2);
end_i = floor(size(x,2)*0.8);
for i=start_i:end_i
    curX=x(:,i);
    c(i)=p.gndclear-ank_y_pos(curX.');
    gradc(:,i,i)=[-1*ank_y_grad(curX.'),zeros(1,p.numJ*2)].';
end
c=c.';
gradc=reshape(gradc,[size(x,1)*size(x,2),size(x,2)]);



end