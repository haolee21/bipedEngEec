function [ceq,gradceq]=fri_cons(x,p)


ceq = zeros(size(x,2),1);
gradceq = zeros(size(x,1),size(x,2),size(x,2));


for i=1:size(x,2)
    ceq(i,1) = fri_toe(x(1:p.numJ*2,i).',x(p.numJ*3+5,i))+fri_heel(x(1:p.numJ*2,i).',x(p.numJ*3+6,i));
    
    gradceq(:,i,i) = dfri_toe(x(1:p.numJ*2,i).',x(p.numJ*3+5,i))+dfri_heel(x(1:p.numJ*2,i).',x(p.numJ*3+6,i));
    
end

gradceq = reshape(gradceq,[size(x,1)*size(x,2),size(x,2)]);
end