clear;clc

load('partDsgn1_lookupTables.mat')

velCMBdy   = [0 0 0];
angVelBdy  = [0 0 0];
velWindBdy = [0 0.1 0.0];

fluidDensity = 1;

ctrlSurfDefl = 0;

aeroSurf(1).refArea        = 1;
aeroSurf(1).aeroCentPosVec = [0.1 -1 0];
aeroSurf(1).spanUnitVec    = [0 1 0];
aeroSurf(1).chordUnitVec   = [1 0 0];
aeroSurf(1).fluidDensity   = 1;
aeroSurf(1).CL = aeroStruct(1).CL;
aeroSurf(1).CD = aeroStruct(1).CD;
aeroSurf(1).alpha = aeroStruct(1).alpha;
aeroSurf(1).GainCL = aeroStruct(1).GainCL;
aeroSurf(1).GainCD = aeroStruct(1).GainCD;

aeroSurf(2).refArea        = 1;
aeroSurf(2).aeroCentPosVec = [0.1 1 0];
aeroSurf(2).spanUnitVec    = [0 1 0];
aeroSurf(2).chordUnitVec   = [1 0 0];
aeroSurf(2).fluidDensity   = 1;
aeroSurf(2).CL = aeroStruct(2).CL;
aeroSurf(2).CD = aeroStruct(2).CD;
aeroSurf(2).alpha = aeroStruct(2).alpha;
aeroSurf(2).GainCL = aeroStruct(2).GainCL;
aeroSurf(2).GainCD = aeroStruct(2).GainCD;

aeroSurf(3).refArea        = 1;
aeroSurf(3).aeroCentPosVec = [4 0 0];
aeroSurf(3).spanUnitVec    = [0 1 0];
aeroSurf(3).chordUnitVec   = [1 0 0];
aeroSurf(3).fluidDensity   = 1;
aeroSurf(3).CL = aeroStruct(3).CL;
aeroSurf(3).CD = aeroStruct(3).CD;
aeroSurf(3).alpha = aeroStruct(3).alpha;
aeroSurf(3).GainCL = aeroStruct(3).GainCL;
aeroSurf(3).GainCD = aeroStruct(3).GainCD;

aeroSurf(4).refArea        = 1;
aeroSurf(4).aeroCentPosVec = [4 0 0.1];
aeroSurf(4).spanUnitVec    = [0 0 1];
aeroSurf(4).chordUnitVec   = [1 0 0];
aeroSurf(4).fluidDensity   = 1;
aeroSurf(4).CL = aeroStruct(4).CL;
aeroSurf(4).CD = aeroStruct(4).CD;
aeroSurf(4).alpha = aeroStruct(4).alpha;
aeroSurf(4).GainCL = aeroStruct(4).GainCL;
aeroSurf(4).GainCD = aeroStruct(4).GainCD;

sim('partitionedLifitngBdy_th')

FBdy.Data

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
