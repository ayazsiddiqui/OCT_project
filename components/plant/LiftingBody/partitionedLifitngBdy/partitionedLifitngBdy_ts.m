% clear;
clc

load('partDsgn1_lookupTables.mat')

vFlow   = [1 0 0];
angVelBdy  = [0 0 0];
euler = [0;0;0];

vCM = [0 0.0 0.0];

fluidDensity = 1;

ctrlSurfDefl = 0;

% sim('partitionedLifitngBdy_th')

% FBdy.Data

% figure('Position',[1          41        1920         963],'Units','Pixels')
% % axes
% % set(gca,'NextPlot','add')
% % grid on
%
% plot3([0 velCMBdy(1)],...
%     [0 velCMBdy(2)],...
%     [0 velCMBdy(3)],...
%     'LineStyle','-','Color','k','LineWidth',2,'DisplayName','Velocity')
% grid on
% hold on
% xlim([-1 1])
% ylim([-1 1])
% zlim([-1 1])
% axis equal
%
% plot3([0 angVelBdy(1)],...
%     [0 angVelBdy(2)],...
%     [0 angVelBdy(3)],...
%     'LineStyle','-','Color','k','LineWidth',2,'DisplayName','Angular Velocity')
%
% plot3([-velWindBdy(1) 0],...
%     [-velWindBdy(2) 0],...
%     [-velWindBdy(3) 0],...
%     'LineStyle','-','Color','b','LineWidth',2,'DisplayName','Wind')
%
% plot3([aeroCentPosVec(1) aeroCentPosVec(1)+FLift.Data(1)],...
%     [aeroCentPosVec(2) aeroCentPosVec(2)+FLift.Data(2)],...
%     [aeroCentPosVec(3) aeroCentPosVec(3)+FLift.Data(3)],...
%     'LineStyle','-','Color','g','LineWidth',2,'DisplayName','Lift')
%
% plot3([aeroCentPosVec(1) aeroCentPosVec(1)+FDrag.Data(1)],...
%     [aeroCentPosVec(2) aeroCentPosVec(2)+FDrag.Data(2)],...
%     [aeroCentPosVec(3) aeroCentPosVec(3)+FDrag.Data(3)],...
%     'LineStyle','-','Color','r','LineWidth',2,'DisplayName','Drag')
%
% plot3([aeroCentPosVec(1) aeroCentPosVec(1)+FBdy.Data(1)],...
%     [aeroCentPosVec(2) aeroCentPosVec(2)+FBdy.Data(2)],...
%     [aeroCentPosVec(3) aeroCentPosVec(3)+FBdy.Data(3)],...
%     'LineStyle','--','Color','b','LineWidth',2,'Marker','x','DisplayName','F Net')
%
% plot3([0 MBdy.Data(1)],...
%     [0 MBdy.Data(2)],...
%     [0 MBdy.Data(3)],...
%     'LineStyle','--','Color','g','LineWidth',2,'Marker','o','DisplayName','Moment')
%
% legend
%
%
%
%
