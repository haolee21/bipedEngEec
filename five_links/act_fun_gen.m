% create the sigma function that we use as activation function

syms q1 q2 q3 q4 q5 th

ypos = end_y_pos([q1,q2,q3,q4,q5]);

sigma = 0.8*(0.5*tanh(400*(th-ypos))+0.5);

matlabFunction(sigma,'file','sigma_out');

dsigma_dq = [diff(sigma,q1);diff(sigma,q2);diff(sigma,q3);diff(sigma,q4);diff(sigma,q5)];

matlabFunction(dsigma_dq,'file','dsigma_dq');