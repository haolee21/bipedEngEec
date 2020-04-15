function [c,ceq,gradc,gradceq] = five_link_nonlcon(x,p)

xgpu = gpuArray(x);
[c_gpu,gradc_gpu]=yposCon(xgpu);
[ceq_gpu1,gradceq_gpu1]=dynConst(xgpu,p);

[ceq2,gradceq2]=gaitLenCon(x,p);


[c,ceq1,gradc,gradceq1] = gather(c_gpu,ceq_gpu1,gradc_gpu,gradceq_gpu1);

ceq = [ceq1;ceq2];
gradceq = [gradceq1,gradceq2];

end