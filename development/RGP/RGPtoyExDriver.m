clear 
% close all
clc

size_grid = 201;
basisDesignSpace = linspace(0,20,size_grid)*(pi/180);
desPtSeq = [0 0 1 2 3 4 15 5 7 6 8 12 8 8 8 8 8]'*(pi/180);
PerfmIndxLib = [0.3 0.28 0.22 0.35 0.18 0.4 0.33 0.42 0.24 0.26 0.19 0.25 0.24 0.29 0.3 0.31 0.19]';
% thetaLib = [2 4 3 1 0 15 5 7 6 8 12 8 8 8]'*(pi/180);
% PerfmIndxLib = [0.3 0.22 0.35 0.18 0.4 0.33 0.42 0.24 0.26 0.19 0.25 0.25 0.22 0.22]';
valPtCounter = length(desPtSeq);
% sigmaSq = 0.025;        % Approximate variance level due to noise
sigmaSq = 0.12;        % Approximate variance level due to noise
hypPara = 0.01;
predVarVec = ones(size_grid,2);

for ii = 1:length(basisDesignSpace)
    for jj = 1:length(basisDesignSpace)
        basisCov(ii,jj) = covFuncEval(basisDesignSpace(ii),basisDesignSpace(jj),hypPara); 
    end 
end 

% Initialize covaraince and mean function
invbasisCov = ones(size_grid,size_grid)/(basisCov);
covUpdateF_prev = basisCov;
meanFuncF_prev = zeros(length(basisDesignSpace),1);
meanFuncBasis = zeros(length(basisDesignSpace),1);

%-------------------------------------------------------------------------
valPtCounterRel = valPtCounter;
for jj = 1:valPtCounterRel
    desParaCur = desPtSeq(jj);
    responseCurr = PerfmIndxLib(jj);
    
    for ii = 1:size_grid
        covCurPara(1,ii) = covFuncEval(desParaCur,basisDesignSpace(ii),hypPara);
    end 
    covCurrPoint = covFuncEval(desParaCur,desParaCur,hypPara);
    
%     % Distance calculation to linearly interpolate between nearest 2 points
%     % to the current point
%     distCur2Grid = sort(abs(basisDesignSpace-desParaCur));   
%     firstDistIndex = find(abs(basisDesignSpace-desParaCur) == distCur2Grid(1));
%     secDistIndex = find(abs(basisDesignSpace-desParaCur) == distCur2Grid(2));
%     if length(secDistIndex) == 2
%         % If the current point corresponds to a grid point, use the mean at
%         % that point, not an interpolation between 2 points
%         trueIndex = round(mean(secDistIndex));
%         meanCurPara = meanFuncF_prev(trueIndex);
%     else 
%         CurParaInterp = [basisDesignSpace(firstDistIndex) basisDesignSpace(secDistIndex)];
%         meanCurParaInterp = [meanFuncF_prev(firstDistIndex) meanFuncF_prev(secDistIndex)];
%         % Acutal function call to linearly interpolate between closest points
%         meanCurPara = interp1(CurParaInterp,meanCurParaInterp,desParaCur);
%     end 
%         
    % Inference part of the algorithm
    covUpdateVec = covCurPara/(basisCov+0.00001*eye(size_grid,size_grid));   % This is J_k in the original paper
       
%     meanFuncPupdate = meanCurPara + covUpdateVec*(meanFuncF_prev - meanFuncBasis);
    meanFuncPupdate = covUpdateVec*meanFuncF_prev;
    covUpdatedP = covCurrPoint + covUpdateVec*(covUpdateF_prev - basisCov)*covUpdateVec';
    
    % Update portion of the algorithm
    kalGainMatrix = covUpdateF_prev*covUpdateVec'*(covUpdatedP + sigmaSq)^-1;
    meanFuncFupdate = meanFuncF_prev + kalGainMatrix*(responseCurr - meanFuncPupdate);
    covUpdatedF = covUpdateF_prev - kalGainMatrix*covUpdateVec*covUpdateF_prev;
    
    % Keep current vector/matrix values for future iterations
    meanFuncF_prev = meanFuncFupdate;
    covUpdateF_prev = covUpdatedF;
    
    % Calculation of prediction variance based on updated covariance matrix
    predVar(:,jj) = diag(covUpdatedF);
    
end 

%--------------------------------------------------------------------------
% Calculation of true predictive mean and variance from traditional GP
% modelling

for ii = 1:valPtCounterRel
    for jj = 1:valPtCounterRel
        trueCovMatrix(ii,jj) = covFuncEval(desPtSeq(ii),desPtSeq(jj),hypPara); 
    end 
end 

for ii = 1:size_grid
    for jj = 1:valPtCounterRel
        covRowVectorTrue(jj,ii) = covFuncEval(basisDesignSpace(ii),desPtSeq(jj),hypPara);
    end 
end 


noiseMat = sigmaSq*eye(valPtCounterRel);
for kk = 1:size_grid
%     trueMeanFunc(kk,1) = covRowVectorTrue(:,kk)'*inv(trueCovMatrix)*PerfmIndxLib(1:valPtCounterRel);
%     trueCovFunc(kk,1) = 1 - covRowVectorTrue(:,kk)'*inv(trueCovMatrix)*covRowVectorTrue(:,kk);
    trueMeanFunc(kk,1) = (covRowVectorTrue(:,kk)'/(trueCovMatrix + noiseMat))*PerfmIndxLib(1:valPtCounterRel);
    trueCovFunc(kk,1) = 1 - (covRowVectorTrue(:,kk)'/(trueCovMatrix + noiseMat))*covRowVectorTrue(:,kk);
end 
%--------------------------------------------------------------------------
% Generation of variance window for GP 
meanFuncFUpperLimit = meanFuncFupdate + predVar(:,end);
meanFuncFLowerLimit = meanFuncFupdate - predVar(:,end);
%--------------------------------------------------------------------------
% Debug plotting of mean function
figure
hold on
plot(basisDesignSpace*(180/pi),meanFuncFupdate(:,end))
scatter(desPtSeq*(180/pi),PerfmIndxLib,'r','filled')
plot(basisDesignSpace*(180/pi),trueMeanFunc,'--k','LineWidth',2)
plot(basisDesignSpace*(180/pi),meanFuncFUpperLimit,'-r','LineWidth',2)
plot(basisDesignSpace*(180/pi),meanFuncFLowerLimit,'-r','LineWidth',2)
xlabel('$\theta$')
ylabel('$J_{inst}$')
grid on
legend('RGP','data','GP')
% legend('RGP','data','GP','$RGP_{UL}$','$RGP_{LL}$')
hold off

% Debug plotting of covariance function
figure
hold on
plot(basisDesignSpace*(180/pi),trueCovFunc,'-k','LineWidth',2)
plot(basisDesignSpace*(180/pi),predVar(:,end),'--r','LineWidth',2)
xlabel('$\theta$')
ylabel('$\sigma ^2 (\theta)$')
legend('GP','RGP')
grid on
hold off

% % Debug plotting of covariance function
% figure
% hold on
% % plot(basisDesignSpace*(180/pi),predVar(:,end))
% plot(basisDesignSpace*(180/pi),trueCovFunc,'-k','LineWidth',2)
% plot(basisDesignSpace*(180/pi),predVarDebug(:,end),'-r','LineWidth',2)
% xlabel('$\theta$')
% ylabel('$\sigma ^2 (\theta)$')
% legend('RGP', 'GP','$RGP"$')
% grid on
% hold off

% figure
% hold on
% plot(basisDesignSpace*(180/pi),predVar(:,end)*sigmaSq + 1)
% plot(basisDesignSpace*(180/pi),trueCovFunc,'-k','LineWidth',2)
% xlabel('$\theta$')
% ylabel('$\sigma ^2 (\theta)$')
% legend('RGP', 'GP')
% grid on
% hold off   
%--------------------------------------------------------------------------    
function [covFuncValue] = covFuncEval(designPt1,designPt2,theta)
    covFuncValue = exp(-1/(2*theta^2)*(abs(designPt1-designPt2))^2);
end 
    