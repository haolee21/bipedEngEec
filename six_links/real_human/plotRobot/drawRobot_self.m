function drawRobot_self(sol,p,fig)
% drawRobot(q,p)
%
% This function draws the robot with configuration q and parameters p
%
% INPUTS:
%   q = [5, 1] = column vector of a single robot configuration
%   p = parameter struct
%


% Compute the points that will be used for plotting

%% define mass for each segments, should be input param oneday
m_foot = 1.05;
m_calf = 3.52;
m_thigh = 7.77;
m_torso = 34.66+6.07;
m_tot = m_foot+m_calf*2+m_thigh*2+m_torso;
startHip = -100;
%%
for frame=1:size(sol,2)
    P = getRobotPos(sol(1,frame),sol(2,frame),sol(3,frame),sol(4,frame),sol(5,frame),sol(6,frame));
   
    G = six_G(sol(1,frame),sol(2,frame),sol(3,frame),sol(4,frame),sol(5,frame),sol(6,frame));
    
   
%     x = [0; P(1:2:end,frame)];
%     y = [0; P(2:2:end,frame)];
    P1 = P(1:2,1);
    P2 = P(1:2,2);
    P3 = P(1:2,3);
    P4 = P(1:2,4);
    P5 = P(1:2,5);
    P6 = P(1:2,6);
    P7 = P(1:2,7);
    x = [0,P1(1),P2(1),P3(1),P4(1),P5(1),P7(1),P6(1)]; % P7 is heel, P6 is toe
    y = [0,P1(2),P2(2),P3(2),P4(2),P5(2),P7(2),P6(2)];
    
    G1 = P(4:5,1);
    G2 = P(4:5,2);
    G3 = P(4:5,3);
    G4 = P(4:5,4);
    G5 = P(4:5,5);
    G6 = P(4:5,6);
    
    com_x =( G1(1)*m_calf+G2(1)*m_thigh+G3(1)*m_torso+G4(1)*m_thigh+G5(1)*m_calf+G6(1)*m_foot)/m_tot;
    
    %plot grf
     if(frame<floor(p.gaitT/p.sampT/2))
        M_ext = Mext_toe(sol(:,frame).');
        sigma = sigma_toe(sol(1:p.numJ,frame).',p.toe_th);
        fext_loc = P6;
        
        
     else
        M_ext = Mext_heel(sol(:,frame).');
        sigma = sigma_heel(sol(1:p.numJ,frame).',p.toe_th);
        fext_loc = P7;
  
    end
    fext = sigma*M_ext*(G.'-sol(p.numJ*2+1:p.numJ*3,frame))/500; %divide 500 is just for ploting 
    
    
    if(fext(2)<0)
        fext = zeros(2,1);
    end
    
    
    
    % Heuristics:
    L = 2;  % Maximum extended leg length
    xBnd = L*[-1.2,1.2];
    yBnd = [-0.5*L,1.5*L ];
    
    % Colors:
    colorGround = [118,62,12]/255;
    colorStance = [200,60,60]/255;
    colorSwing = [60,60,200]/255;
    colorTorso = [160, 80, 160]/255;
    
  
    % Set up the figure
    hold off;
    clf(fig);
    % Plot the ground:
    plot(xBnd,[p.toe_th,p.toe_th],'LineWidth',6,'Color',colorGround);
    
    hold on;
    
    % Plot the links:
    plot(x(1:3),y(1:3),'LineWidth',4,'Color',colorStance);
    plot(x(3:4),y(3:4),'LineWidth',4,'Color',colorTorso);
    plot(x([3,5,6,7,8]),y([3,5,6,7,8]),'LineWidth',4,'Color',colorSwing);
    
    % Plot the joints:
    plot(0, 0,'k.','MarkerSize',30);
    plot(P1(1), P1(2),'k.','MarkerSize',30);
    plot(P2(1), P2(2),'k.','MarkerSize',30);
    plot(P3(1), P3(2),'k.','MarkerSize',30);
    plot(P4(1), P4(2),'k.','MarkerSize',30);
    plot(P5(1), P5(2),'k.','MarkerSize',30);
    plot(P6(1), P6(2),'k.','MarkerSize',30);
    
    % Plot the CoM:
    plot(G1(1), G1(2),'ko','MarkerSize',8,'LineWidth',2);
    plot(G2(1), G2(2),'ko','MarkerSize',8,'LineWidth',2);
    plot(G3(1), G3(2),'ko','MarkerSize',8,'LineWidth',2);
    plot(G4(1), G4(2),'ko','MarkerSize',8,'LineWidth',2);
    plot(G5(1), G5(2),'ko','MarkerSize',8,'LineWidth',2);
    plot(G6(1), G6(2),'ko','MarkerSize',8,'LineWidth',2);
    
    plot(com_x,0,'go','MarkerSize',8,'LineWidth',2);
    
    quiver(fext_loc(1),fext_loc(2),fext(1),fext(2),'color',[0.2,1,0],'LineStyle','--','LineWidth',2);
    yline(p.toe_th,'-','Toe Height');
    yline(p.head_h,'-','Head Height');
    if(startHip<-50)
        startHip = P2(1);
    end
    xline(startHip,'-','Start');
    xline(startHip-p.hipLen,'-','End');
    % Format the axis:
    axis([xBnd,yBnd]); axis equal; axis off;
    pause(0.002);
end
end