%% quantification of cluster parameters

function [clustTable,clusterMemberInd] = quantClustProp(dataIn,clustCent,point2cluster,clustMembsCell,filtCutoff)

% figure
% hold on
clusterMemberInd = [];
for k = 1:length(clustCent(1,:))
numInClust(k) = length(find(point2cluster==k));
% end

if numInClust(k)>= filtCutoff
% nonSmallClust = find(numInClust>=2);

% numClustFilt = length(nonSmallClust);
% ClustInd = nonSmallClust;
% ClustTable = table(ClustInd,'VariableNames',{'Index'});
% figure 
% plot(testSubSet(:,1),testSubSet(:,2),'.')
% hold on
% axis equal
% for k = 1:numClustFilt
%     ClustInd = nonSmallClust(k);
ptInds = clustMembsCell{k};
% ptInds = find(point2cluster==k);
CovMat = cov(dataIn(ptInds,:));
[V,D] = eig(CovMat);
majPhi = atan(V(2,2)/V(1,2));
minPhi = atan(V(2,1)/V(1,1));
majVec = [cos(majPhi);sin(majPhi)];
minorVec = [cos(minPhi);sin(minPhi)];
chiCrit = 5.991;
majAxLen = sqrt(chiCrit*D(2,2));
minAxLen = sqrt(chiCrit*D(1,1));
center = clustCent(:,k);
ellipArea = pi*majAxLen*minAxLen;

reCenter = mean(dataIn(ptInds,:));
reMeanX = reCenter(1);reMeanY=reCenter(2);

testPlotBool = 0;
if testPlotBool == 1
    
    axis equal
    plot(dataIn(ptInds,1)-center(1),dataIn(ptInds,2)-center(2),'.')
    hold on
    clusterMemberInd = [clusterMemberInd; dataIn(ptInds,1)-center(1),dataIn(ptInds,2)-center(2)];
    %plot(center(1),center(2),'r*')
%     t = 0:.1:2*pi;
% %     xEl = center(1)+majAxLen.*cos(t).*cos(majPhi)-minAxLen.*sin(t).*sin(majPhi);
% %     yEl = center(2)+majAxLen.*cos(t).*cos(majPhi)+minAxLen.*sin(t).*sin(majPhi);
%     xEl = reMeanX+majAxLen.*cos(t).*cos(majPhi)-minAxLen.*sin(t).*sin(majPhi);
%     yEl = reMeanY+majAxLen.*cos(t).*sin(majPhi)+minAxLen.*sin(t).*cos(majPhi);
%     plot(xEl,yEl,'g')
end

% Find Nearest Neigbor
[~,nnd] = knnsearch(clustCent',clustCent(:,k)','K',2);
meanNNk1 = mean(nnd(2:end));
[~,nnd] = knnsearch(clustCent',clustCent(:,k)','K',4);
meanNNk3 = mean(nnd(2:end));
[~,nnd] = knnsearch(clustCent',clustCent(:,k)','K',6);
meanNNk5 = mean(nnd(2:end));
[~,nnd] = knnsearch(clustCent',clustCent(:,k)','K',10);
meanNNk9 = mean(nnd(2:end));

[~,ConvHullV] = convhull(dataIn(ptInds,1),dataIn(ptInds,2));


clustTable(k,:) = [majAxLen,minAxLen,center(1),center(2),ellipArea,...
                    majPhi,numInClust(k),meanNNk1,meanNNk3,meanNNk5,...
                    meanNNk9,reMeanX,reMeanY,ConvHullV];
else
    clustTable(k,:) = [NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN];
end
end