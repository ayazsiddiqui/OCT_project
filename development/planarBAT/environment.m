clear
clc
format compact
rng(2);

%% environment
heights = [0:100:1000]';
meanFlow = 10;

% fit a polynomial to dummy data
ht = [0:200:1000]';
ft2 =  normrnd(meanFlow,0.15*meanFlow,size(ht));
ft2 = [5;5.5;6.5;7;6;5.5];

pt = polyfit(ht,ft2,3);

Flows = polyval(pt,heights);

% simTime in minutes
simTime = 30;
FlowInt = 60*(0:5:simTime);
nS = length(FlowInt);

for ii = 1:nS
    
    for jj = 1:length(pt)
        rdNum = rand;
        if rdNum < 0.33
            mut = -1;
        elseif rdNum >= 0.33 && rdNum < 0.66
            mut = 0;
        else
            mut = 1;
        end
        pt(jj) = pt(jj)*(1 + mut*0.1);
    end
    
    Flows(:,:,ii) = polyval(pt,heights);
    
end

%% train GP
gp = gaussianProcess(2,'kernel','squaredExponential','acquisitionFunction','upperConfidenceBound');
if strcmpi(class(gp.acquisitionFunction),'acquisitionFunctions.upperConfidenceBound')
    gp.acquisitionFunction.explorationFactor = 1;
end
gp.kernel.noiseVariance = 1*0.05;

tVals = reshape(transpose(FlowInt'.*ones(length(FlowInt),length(heights))),[],1)';
trainDsgns = [repmat(heights',1,nS);tVals];
trainFval = Flows(:);

% step 1: optimize hyper parameters
initialGuess = rand(1+gp.noInputs,1);
trainOpHyp = gp.optimizeHyperParameters(trainDsgns,trainFval,initialGuess);

gp.kernel.covarianceAmp = trainOpHyp(1);
gp.kernel.lengthScale = trainOpHyp(2:end);

trainCovMat = gp.buildCovarianceMatrix(trainDsgns,trainDsgns);

%% posterior
% posterior
nPost = 10;
postAlts = ((heights(end)-heights(1)).*rand(1,nPost) + heights(1));
postAlts = reshape(heights,1,[]);
postTime = 2500*ones(size(postAlts));

postDsgns = [postAlts;postTime];

[predMean,predVar] = gp.calcPredictiveMeanAndVariance(postDsgns,trainDsgns,trainCovMat,trainFval);


%% plot
figure(1)
xlabel('Flow speed (m/s)')
ylabel('Height (m)')
hold on
grid on

for ii = 1:nS
    if ii>1
        delete(p1)
    end
    p1 = plot(Flows(:,:,ii),heights);
    title(sprintf('Time = %0.1f s',FlowInt(ii)))
%     waitforbuttonpress
end




