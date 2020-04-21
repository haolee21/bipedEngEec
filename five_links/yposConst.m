function [c,c_grad] = yposConst(x)
% the input x is a gpuArray

idxList = 1:1:size(x,1);

cfun =@(i)end_y_pos(x,i);% this i is for arrayfun, end_y_pos does not take i as input
gradfun = @(i)end_y_poss_grad(x,i);

c = arrayfun(cfun,idxList);
c_grad = reshape(cell2mat(arrayfun(gradfun,idxList,'UniformOutput',false)),[5,size(x,1)]);

end