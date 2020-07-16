function [c,gradc]=slack_cons(x,p)

c_toe = zeros(size(x,2),1);
gradc_toe = zeros(size(x,1),size(x,2),size(x,2));
c_heel =  zeros(size(x,2),1);
gradc_heel = zeros(size(x,1),size(x,2),size(x,2));

for i=1:size(x,2)
    c_toe(i,1)=s_toe_con(x(1:p.numJ,i).',x(3*p.numJ+5,i),p.toe_th);
    gradc_toe(:,i,i)=ds_toe_con(x(1:p.numJ,i).',x(3*p.numJ+5,i),p.toe_th);
    c_heel(i,1)=s_heel_con(x(1:p.numJ,i).',x(3*p.numJ+5,i),p.toe_th);
    gradc_heel(:,i,i)=ds_heel_con(x(1:p.numJ,i).',x(3*p.numJ+5,i),p.toe_th);
    
end
c=[c_toe;c_heel];
gradc = [reshape(gradc_toe,[size(x,1)*size(x,2),size(x,2)]),reshape(gradc_heel,[size(x,1)*size(x,2),size(x,2)])];



end