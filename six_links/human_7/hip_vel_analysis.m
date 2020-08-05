addpath('robotGen/dyn');
addpath('robotGen/grad');
addpath('robotGen/grf/discrete');

data1 = load('07280017').result;
hip_vx = zeros(1,size(data1.x,2)-1);
x = data1.x;
p = data1.param;
for i=1:size(data1.x,2)-1
    q1 = x(1:p.numJ,i);
    q2 = x(1:p.numJ,i+1);
    q =(q1+q2)/2;
    dq = (q2-q1)/p.sampT;
    hip_vx(i)=hip_vel_x(q.',dq.',p.sampT);
    
end