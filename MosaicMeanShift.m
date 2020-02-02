%% mosaic mean shift
function [clustCent_all,point2cluster_all,clustMembsCell_all] = MosaicMeanShift(testDat,bandwidth,borderPerc,MStype)
% load('gm_20151207.mat')
% testDat = pPKAFI_cell001(:,2:3);
% load('20151223_STORMpts.mat')
% testDat = noprebleach_HeLa_pPKAFLINC_120_cell001(:,2:3);
% testSubSet = testDat;
if MStype == 1
    MSfcn = @(dataP,BW,plotFl)MeanShiftClusterGauss(dataP,BW,plotFl);
else
    MSfcn = @(dataP,BW,plotFl)MeanShiftCluster(dataP,BW,plotFl);
end
maxX = max(testDat(:,1));
minX = min(testDat(:,1));
maxY = max(testDat(:,2));
minY = min(testDat(:,2));
rangeX = maxX-minX;
rangeY = maxY-minY;
midX = (maxX+minX)/2;
midY = (maxY+minY)/2;



subSetInd_Q1= find(testDat(:,1)>=midX&testDat(:,2)>=midY);
dataSubset_Q1 = testDat(subSetInd_Q1,:);

subSetInd_Q2= find(testDat(:,1)<midX&testDat(:,2)>=midY);
dataSubset_Q2 = testDat(subSetInd_Q2,:);

subSetInd_Q3= find(testDat(:,1)<midX&testDat(:,2)<midY);
dataSubset_Q3 = testDat(subSetInd_Q3,:);

subSetInd_Q4= find(testDat(:,1)>=midX&testDat(:,2)<midY);
dataSubset_Q4 = testDat(subSetInd_Q4,:);

% figure
% hold on
% plot(dataSubset_Q1(:,1),dataSubset_Q1(:,2),'.')
% plot(dataSubset_Q2(:,1),dataSubset_Q2(:,2),'.k')
% plot(dataSubset_Q3(:,1),dataSubset_Q3(:,2),'.r')
% plot(dataSubset_Q4(:,1),dataSubset_Q4(:,2),'.g')

% borderPerc = 0.05;

subSetInd_Q1_int= find(dataSubset_Q1(:,1)>=(midX+borderPerc*rangeX)&dataSubset_Q1(:,2)>=(midY+borderPerc*rangeY));
subSetInd_Q1_bord = setdiff([1:length(subSetInd_Q1(:,1))],subSetInd_Q1_int);
% dataSubset_Q1 = testDat(subSetInd_Q1,:);

subSetInd_Q2_int= find(dataSubset_Q2(:,1)<=(midX-borderPerc*rangeX)&dataSubset_Q2(:,2)>=(midY+borderPerc*rangeY));
subSetInd_Q2_bord = setdiff([1:length(subSetInd_Q2(:,1))],subSetInd_Q2_int);

subSetInd_Q3_int= find(dataSubset_Q3(:,1)<=(midX-borderPerc*rangeX)&dataSubset_Q3(:,2)<=(midY-borderPerc*rangeY));
subSetInd_Q3_bord = setdiff([1:length(subSetInd_Q3(:,1))],subSetInd_Q3_int);

subSetInd_Q4_int= find(dataSubset_Q4(:,1)>=(midX+borderPerc*rangeX)&dataSubset_Q4(:,2)<=(midY-borderPerc*rangeY));
subSetInd_Q4_bord = setdiff([1:length(subSetInd_Q4(:,1))],subSetInd_Q4_int);

% figure
% hold on
% plot(dataSubset_Q1(subSetInd_Q1_int,1),dataSubset_Q1(subSetInd_Q1_int,2),'.')
% % plot(dataSubset_Q1(subSetInd_Q1_bord,1),dataSubset_Q1(subSetInd_Q1_bord,2),'.r')
% 
% plot(dataSubset_Q2(subSetInd_Q2_int,1),dataSubset_Q2(subSetInd_Q2_int,2),'.k')
% plot(dataSubset_Q3(subSetInd_Q3_int,1),dataSubset_Q3(subSetInd_Q3_int,2),'.r')
% plot(dataSubset_Q4(subSetInd_Q4_int,1),dataSubset_Q4(subSetInd_Q4_int,2),'.g')

%% Do the mean shifts

% tic
[clustCent_Q1,point2cluster_Q1,clustMembsCell_Q1] = MSfcn(dataSubset_Q1',bandwidth,0);
fprintf('Finished Q1 ...')
% toc
% tic
[clustCent_Q2,point2cluster_Q2,clustMembsCell_Q2] = MSfcn(dataSubset_Q2',bandwidth,0);
fprintf(' Q2 ...')
% toc
% tic
[clustCent_Q3,point2cluster_Q3,clustMembsCell_Q3] = MSfcn(dataSubset_Q3',bandwidth,0);
fprintf(' Q3 ...')
% toc
% tic
[clustCent_Q4,point2cluster_Q4,clustMembsCell_Q4] = MSfcn(dataSubset_Q4',bandwidth,0);
fprintf(' Q4 ...')
% toc

%% pull out ones not in border
boardClust_Q1 = unique(point2cluster_Q1(subSetInd_Q1_bord));
boardClustInd_Q1 = ismember(point2cluster_Q1,boardClust_Q1);

boardClust_Q2 = unique(point2cluster_Q2(subSetInd_Q2_bord));
boardClustInd_Q2 = ismember(point2cluster_Q2,boardClust_Q2);

boardClust_Q3 = unique(point2cluster_Q3(subSetInd_Q3_bord));
boardClustInd_Q3 = ismember(point2cluster_Q3,boardClust_Q3);

boardClust_Q4 = unique(point2cluster_Q4(subSetInd_Q4_bord));
boardClustInd_Q4 = ismember(point2cluster_Q4,boardClust_Q4);

% figure; hold on
% plot(dataSubset_Q1(~boardClustInd_Q1,1),dataSubset_Q1(~boardClustInd_Q1,2),'.','Color',[0 0 .7])
% plot(dataSubset_Q1(boardClustInd_Q1,1),dataSubset_Q1(boardClustInd_Q1,2),'.','Color',[.2 .2 1])
% 
% plot(dataSubset_Q2(~boardClustInd_Q2,1),dataSubset_Q2(~boardClustInd_Q2,2),'.','Color',[0 .7 0])
% plot(dataSubset_Q2(boardClustInd_Q2,1),dataSubset_Q2(boardClustInd_Q2,2),'.','Color',[.2 1 .2])
% 
% plot(dataSubset_Q3(~boardClustInd_Q3,1),dataSubset_Q3(~boardClustInd_Q3,2),'.','Color',[.7 0 0])
% plot(dataSubset_Q3(boardClustInd_Q3,1),dataSubset_Q3(boardClustInd_Q3,2),'.','Color',[1 .2 .2])
% 
% plot(dataSubset_Q4(~boardClustInd_Q4,1),dataSubset_Q4(~boardClustInd_Q4,2),'.','Color',[.1 .1 .1])
% plot(dataSubset_Q4(boardClustInd_Q4,1),dataSubset_Q4(boardClustInd_Q4,2),'.','Color',[.7 .7 .7])


midPoints = [];
midPoints = [midPoints;subSetInd_Q1(boardClustInd_Q1)];
midPoints = [midPoints;subSetInd_Q2(boardClustInd_Q2)];
midPoints = [midPoints;subSetInd_Q3(boardClustInd_Q3)];
midPoints = [midPoints;subSetInd_Q4(boardClustInd_Q4)];

% figure
% plot(testDat(midPoints,1),testDat(midPoints,2),'.')

% tic
[clustCent_mid,point2cluster_mid,clustMembsCell_mid] =MSfcn(testDat(midPoints,:)',bandwidth,0);
fprintf('mid ... collecting vals ... ')
% toc

% tic
% [clustCent_all,point2cluster_all,clustMembsCell_all] = MeanShiftClusterGauss(testDat',150,0);
% toc

% figure
% hold on
% plot(testDat(:,1),testDat(:,2),'.')
% 
% plot(clustCent_Q1(1,:),clustCent_Q1(2,:),'*g')
% plot(clustCent_Q2(1,:),clustCent_Q2(2,:),'*k')
% plot(clustCent_Q3(1,:),clustCent_Q3(2,:),'*c')
% plot(clustCent_Q4(1,:),clustCent_Q4(2,:),'*m')
% plot(clustCent_mid(1,:),clustCent_mid(2,:),'*r')

%% Gather all the clusters
mapClustInd = [];
point2cluster_all = zeros(1,length(testDat));

nClust = 1:length(clustCent_Q1);
clustmatch = find(~ismember(nClust,boardClust_Q1));
mapClustInd = 1:length(clustmatch);
clustMap = [mapClustInd;clustmatch];
nOld = length(mapClustInd);
clustCent_all(:,mapClustInd) = clustCent_Q1(:,clustmatch);
clustMembsCell_all(mapClustInd,1) = cellfun(@(x)subSetInd_Q1(x),clustMembsCell_Q1(clustmatch), 'UniformOutput', false);
pts2map = point2cluster_Q1(~boardClustInd_Q1);
mapVals = arrayfun(@(x)clustMap(1,find(clustMap(2,:)==x,1)),pts2map);
point2cluster_all(subSetInd_Q1(~boardClustInd_Q1)) = mapVals;

nClust = 1:length(clustCent_Q2);
clustmatch = find(~ismember(nClust,boardClust_Q2));
mapClustInd = 1+nOld:length(clustmatch)+nOld;
clustMap = [mapClustInd;clustmatch];
nOld = length(clustmatch)+nOld;
clustCent_all(:,mapClustInd) = clustCent_Q2(:,clustmatch);
clustMembsCell_all(mapClustInd,1) = cellfun(@(x)subSetInd_Q2(x),clustMembsCell_Q2(clustmatch), 'UniformOutput', false);
pts2map = point2cluster_Q2(~boardClustInd_Q2);
mapVals = arrayfun(@(x)clustMap(1,find(clustMap(2,:)==x,1)),pts2map);
point2cluster_all(subSetInd_Q2(~boardClustInd_Q2)) = mapVals;

nClust = 1:length(clustCent_Q3);
clustmatch = find(~ismember(nClust,boardClust_Q3));
mapClustInd = 1+nOld:length(clustmatch)+nOld;
clustMap = [mapClustInd;clustmatch];
nOld = length(clustmatch)+nOld;
clustCent_all(:,mapClustInd) = clustCent_Q3(:,clustmatch);
clustMembsCell_all(mapClustInd,1) = cellfun(@(x)subSetInd_Q3(x),clustMembsCell_Q3(clustmatch), 'UniformOutput', false);
pts2map = point2cluster_Q3(~boardClustInd_Q3);
mapVals = arrayfun(@(x)clustMap(1,find(clustMap(2,:)==x,1)),pts2map);
point2cluster_all(subSetInd_Q3(~boardClustInd_Q3)) = mapVals;

nClust = 1:length(clustCent_Q4);
clustmatch = find(~ismember(nClust,boardClust_Q4));
mapClustInd = 1+nOld:length(clustmatch)+nOld;
clustMap = [mapClustInd;clustmatch];
nOld = length(clustmatch)+nOld;
clustCent_all(:,mapClustInd) = clustCent_Q4(:,clustmatch);
clustMembsCell_all(mapClustInd,1) = cellfun(@(x)subSetInd_Q4(x),clustMembsCell_Q4(clustmatch), 'UniformOutput', false);
pts2map = point2cluster_Q4(~boardClustInd_Q4);
mapVals = arrayfun(@(x)clustMap(1,find(clustMap(2,:)==x,1)),pts2map);
point2cluster_all(subSetInd_Q4(~boardClustInd_Q4)) = mapVals;

nClust = 1:length(clustCent_mid);
clustmatch = nClust;
mapClustInd = 1+nOld:length(clustmatch)+nOld;
clustMap = [mapClustInd;clustmatch];
nOld = length(clustmatch)+nOld;
clustCent_all(:,mapClustInd) = clustCent_mid(:,clustmatch);
clustMembsCell_all(mapClustInd,1) = cellfun(@(x)midPoints(x),clustMembsCell_mid(clustmatch), 'UniformOutput', false);;
pts2map = point2cluster_mid;
mapVals = arrayfun(@(x)clustMap(1,find(clustMap(2,:)==x,1)),pts2map);
point2cluster_all(midPoints) = mapVals;
fprintf(' done \n')
% QuadInd = zeros(length(testDat),1);


% figure
% hold on
% plot(testDat(:,1),testDat(:,2),'.')
% plot(clustCent_mid(1,:),clustCent_mid(2,:),'*r')
% plot(clustCent_all(1,:),clustCent_all(2,:),'*g')
% 
% figure
% hold on
% plot(testDat(:,1),testDat(:,2),'.')
% plot(clustCent_all(1,:),clustCent_all(2,:),'*g')
% plot(clustCent_mid(1,:),clustCent_mid(2,:),'*r')
% plot(clustCent_Q1(1,:),clustCent_Q1(2,:),'*r')
% plot(clustCent_Q2(1,:),clustCent_Q2(2,:),'*r')
% plot(clustCent_Q3(1,:),clustCent_Q3(2,:),'*r')
% plot(clustCent_Q4(1,:),clustCent_Q4(2,:),'*r')

