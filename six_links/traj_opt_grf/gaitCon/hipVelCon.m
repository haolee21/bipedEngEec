function [c,gradc]=hipVelCon(x,p)
% hip velocity has to be negative since it is going on the negative
% direction

% c_cell = cellfun(@(x)hip_x_vel(x,p),num2cell(x,1));
% c = cell2mat(c_cell);
% gradc_cell = cellfun(@(x)dHip_x_vel(x,p),num2cell(x,1));
% gradc = cell2mat(gradc_cell);
c = zeros(size(x,2),1);
gradc = zeros(size(x,1),size(x,2),size(x,2));

if(~isnumeric(x))
    c = sym(zeros(size(x,2),1));
    gradc=sym(zeros(size(x,1),size(x,2),size(x,2)));
end

for i =1:size(x,2)
    c(i,1) = hip_x_vel(x(:,i).');
    gradc(1:2*p.numJ,i,i)=dHip_x_vel(x(:,i).');
    
end
gradc = reshape(gradc,[size(x,1)*size(x,2),size(x,2)]);
end