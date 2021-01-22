%
clear all
close all;
clc;


%parameters
timeStepId = 1;             %pick  one bewteeen 1 and 100
trainStepSize = 3;
testStepSize = 10;
testStartId = round(9*rand(1,1))+1;

%script outputs
trainX = [];
trainy = [];
trainFieldIds = [];
testX = [];
testy = [];
testFieldIds = [];

%load and visualize the data
load('field.mat')
figure;
imagesc(TempField(:,:,timeStepId));
title('temperature measurement of the east coast in the states');

%****************************
%%train set

for rIds=1:trainStepSize:size(LatScale,1)
    for cIds=1:trainStepSize:size(LonScale,1)
        if ~isnan(TempField(rIds, cIds, timeStepId))
            trainX = [trainX; LatScale(rIds), LonScale(cIds)];
            trainy = [trainy; TempField(rIds, cIds, timeStepId)];
            trainFieldIds = [trainFieldIds; rIds, cIds];
        end
    end
end
%%test set
for rIds=testStartId:testStepSize:size(LatScale,1)
    for cIds=testStartId:testStepSize:size(LonScale,1)
        if ~isnan(TempField(rIds, cIds, timeStepId))
            testX = [testX; LatScale(rIds), LonScale(cIds)];
            testy = [testy; TempField(rIds, cIds, timeStepId)];
            testFieldIds = [testFieldIds; rIds, cIds];
        end
    end
end

%dimension should be D times n
trainX = double(trainX'); 
trainFieldIds = trainFieldIds';
testX = double(testX');
testFieldIds = testFieldIds';

data.xTrain = trainX;
data.yTrain = trainy;
data.xTest= testX;
data.yTest = testy;
save('tempData.mat', '-struct', 'data')

%% visualize the train and test sets
figure;
subplot(2,1,1);
imagesc(TempField(trainFieldIds(1,:), trainFieldIds(2,:), timeStepId));
subplot(2,1,2);
imagesc(TempField(testFieldIds(1,:), testFieldIds(2,:), timeStepId));

%% save the data 
save('TempFieldDataSubset.mat');

