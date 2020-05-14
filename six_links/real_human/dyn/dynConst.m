function [ceq,ceq_grad] = dynConst(x,param)
% calculate the constraints of following system dynamic

numJ = param.numJ; %here we have to define number of joints since this is related to the dynamic equations



%% equality constraints 
% e.g., q(1) - q(0) = dt*0.5*(f1+f2)

q = x(1:numJ,:); % the first constraints are only calculated by q and qdot
dq = x(numJ+1:numJ*2,:);

ceq_p11 = q(:,2:end)-q(:,1:end-1);
ceq_p11 = reshape(ceq_p11,size(ceq_p11,1)*size(ceq_p11,2),1);
ceq_p12 = dq(:,2:end)-dq(:,1:end-1);
ceq_p12 = reshape(ceq_p12,size(ceq_p12,1)*size(ceq_p12,2),1);

ceq_p1 = [ceq_p11;ceq_p12];


ceq_p21 = reshape(dq(:,1:end-1)+dq(:,2:end),[size(ceq_p11,1),1]);





x_extend = [x(:,2:end);x(:,1:end-1)];%because parallel toolbox will be faster if we only use one column


df_dxt = zeros(size(x,2),size(x,1),numJ);
% df_dxt = gpuArray(df_dxt);

df_dxt_k = zeros(size(x,2)-1,size(x,1),numJ);
% df_dxt_k = gpuArray(df_dxt_k);
% the following loop only solve from f(x2,x1), f(x1) is ignored, we need to
% do it first
ceq_p22 = zeros(size(x,2),numJ);% this dimension is not correct, but we need to make it simplier in order to use parfor, it is [f1,f2,f3....], not [f1+f2,f2+f3] for now
% ceq_p22 = gpuArray(ceq_p22);


[out,grad_t,~] = f_x(x(:,1).',param,1);
df_dxt(1,:,:)=grad_t;
ceq_p22(1,:) = out;


for i=1:size(x,2)-1
    
    [out,grad_t,grad_tk] = f_x(x_extend(:,i).',param,i);
    df_dxt(i+1,:,:) = grad_t;
    df_dxt_k(i,:,:) = grad_tk;
    ceq_p22(i+1,:)=out;  %although it should be 2n(m-1)x1 vector, to use parallel toolbox, we assign to different column first
    
end
% reshape ceq_p22
ceq_p22 = ceq_p22(1:end-1,:)+ceq_p22(2:end,:);
ceq_p22 = reshape(ceq_p22.',[size(ceq_p21,1),size(ceq_p21,2)]);

ceq = ceq_p1-0.5*param.sampT*[ceq_p21;ceq_p22];

% form the gradient
ceq_grad = zeros(size(x,1),size(x,2),(size(x,2)-1)*2*numJ); %total gradient, 3n x m x 2n(m-1)
% ceq_grad = gpuArray(ceq_grad); 

%calculate the first gradient since it is a special case
p_temp = zeros(size(x,1),2*numJ*(size(x,2)-1));
% p_temp = gpuArray(p_temp);

%p1 terms
p_temp(1:numJ,1:numJ) =p_temp(1:numJ,1:numJ)-eye(numJ);
p_temp(numJ+1:2*numJ,(size(x,2)-1)*numJ+1:size(x,2)*numJ) = p_temp(numJ+1:2*numJ,(size(x,2)-1)*numJ+1:size(x,2)*numJ)-eye(numJ);
%p2 terms
p_temp(numJ+1:2*numJ,1:numJ) = p_temp(numJ+1:2*numJ,1:numJ)-0.5*param.sampT*eye(numJ);
p_temp(:,(size(x,2)-1)*numJ+1:(size(x,2)+1)*numJ) = p_temp(:,(size(x,2)-1)*numJ+1:(size(x,2)+1)*numJ)-0.5*param.sampT*[reshape(df_dxt(1,:,:),[3,1]*numJ)+reshape(df_dxt_k(1,:,:),[3,1]*numJ),reshape(df_dxt_k(1,:,:),[3,1]*numJ)];
ceq_grad(:,1,:) = p_temp;

sampT = param.sampT;
clear p_temp;
for i=2:size(x,2)-2 % the first and the last 2 are special cases
    p_temp = zeros(size(x,1),2*numJ*(size(x,2)-1));
%     p_temp = gpuArray(p_temp);
    
    %p1 terms
    p_temp(1:numJ,(i-2)*numJ+1:i*numJ) = [eye(numJ),-eye(numJ)];
    p_temp(numJ+1:2*numJ,(i+size(x,2)-3)*numJ+1:(i+size(x,2)-1)*numJ) = [eye(numJ),-eye(numJ)];
    
    %p2 terms
    p_temp(:,(i-2)*numJ+1:i*numJ) = p_temp(:,(i-2)*numJ+1:i*numJ)-0.5*sampT* [zeros(numJ),zeros(numJ);eye(numJ),eye(numJ);zeros(numJ),zeros(numJ)];
    p_temp(:,(i+size(x,2)-3)*numJ+1:(i+size(x,2))*numJ) =p_temp(:,(i+size(x,2)-3)*numJ+1:(i+size(x,2))*numJ)-0.5*sampT*[reshape(df_dxt(i,:,:),[3,1]*numJ),reshape(df_dxt_k(i,:,:),[3,1]*numJ)+reshape(df_dxt(i,:,:),[3,1]*numJ),reshape(df_dxt_k(i,:,:),[3,1]*numJ)];

    ceq_grad(:,i,:) = p_temp;
    
end

% add the gradient of the last 2 time step

%   m-1
% the p1 terms are similar, we take it out because we need to use parfor 
p_temp = zeros(size(x,1),2*numJ*(size(x,2)-1));
% p_temp = gpuArray(p_temp);

p_temp(1:numJ,(size(x,2)-3)*numJ+1:(size(x,2)-1)*numJ) =p_temp(1:numJ,(size(x,2)-3)*numJ+1:(size(x,2)-1)*numJ)+ [eye(numJ),-eye(numJ)];
p_temp(numJ+1:2*numJ,(2*size(x,2)-4)*numJ+1:(2*size(x,2)-2)*numJ) = p_temp(numJ+1:2*numJ,(2*size(x,2)-4)*numJ+1:(2*size(x,2)-2)*numJ)+[eye(numJ),-eye(numJ)];
%   p2 terms
p_temp(:,(size(x,2)-3)*numJ+1:(size(x,2)-1)*numJ) = p_temp(:,(size(x,2)-3)*numJ+1:(size(x,2)-1)*numJ)-0.5*sampT*[zeros(numJ),zeros(numJ);eye(numJ),eye(numJ);zeros(numJ),zeros(numJ)];
p_temp(:,(2*size(x,2)-4)*numJ+1:(2*size(x,2)-2)*numJ) =p_temp(:,(2*size(x,2)-4)*numJ+1:(2*size(x,2)-2)*numJ)-0.5*sampT*[reshape(df_dxt(size(x,2)-1,:,:),[3,1]*numJ),reshape(df_dxt_k(size(x,2)-1,:,:),[3,1]*numJ)+reshape(df_dxt(size(x,2)-1,:,:),[3,1]*numJ)];

ceq_grad(:,size(x,2)-1,:) = p_temp;
clear p_temp;
% m
% p1 terms
p_temp = zeros(size(x,1),2*numJ*(size(x,2)-1));
% p_temp = gpuArray(p_temp);

p_temp(1:numJ,(size(x,2)-2)*numJ+1:(size(x,2)-1)*numJ) = p_temp(1:numJ,(size(x,2)-2)*numJ+1:(size(x,2)-1)*numJ)+eye(numJ);
p_temp(numJ+1:2*numJ,(2*size(x,2)-3)*numJ+1:(2*size(x,2)-2)*numJ) = p_temp(numJ+1:2*numJ,(2*size(x,2)-3)*numJ+1:(2*size(x,2)-2)*numJ)+eye(numJ);
% p2 terms
p_temp(numJ+1:2*numJ,(size(x,2)-2)*numJ+1:(size(x,2)-1)*numJ) = p_temp(numJ+1:2*numJ,(size(x,2)-2)*numJ+1:(size(x,2)-1)*numJ)-0.5*param.sampT* eye(numJ);
p_temp(:,(2*size(x,2)-3)*numJ+1:(2*size(x,2)-2)*numJ) =p_temp(:,(2*size(x,2)-3)*numJ+1:(2*size(x,2)-2)*numJ)-0.5*param.sampT*reshape(df_dxt(size(x,2),:,:),[3,1]*numJ);
ceq_grad(:,size(x,2),:) = p_temp;


%reshape for calculation
ceq_grad = reshape(ceq_grad,[size(x,1)*size(x,2),size(ceq_grad,3)]);
end