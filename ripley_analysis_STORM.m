function [xout, yout, frameout, H, Rmax]=ripley_analysis_STORM(data,dims,rs)
    % data: Full STORM data set such that the first column is the frame number, the
    % second column is x coordinate data, and 3rd colunm is y coordinate
    % data
    
    % dims: Dimension of ROI that you would like to do the analysis on.
    % rs: An array of the different values of R for which you would like to do
    % the Ripley K calculation
    % xout: X-coordinates of points within ROI
    % yout: Y-coordinates of points within ROI
    % frameout: frames of data points within ROI
    % H: array of values of L(R)-R for all values of rs
    % Rmax: the value of R for which H is maximum 
    
    % by Brian Ross
    %frames = data(:,1);
    x = data(:,1);
    y = data(:,2);
    xmin = dims(1);
    xmax = dims(2);
    ymin = dims(3);
    ymax = dims(4);
    j=1;
    table2=[];                  
    for i=1:size(data,1);
        if data(i,1)>xmin && data(i,1)<xmax && data(i,2)>ymin && data(i,2)<ymax %& table(i,3)<maxstddev & table(i,4)<maxstddev;
            table2(j,1)=1;
            table2(j,2)=data(i,1);
            table2(j,3)=data(i,2);
            j=j+1;
        end
    end
    xout = table2(:,2);
    yout = table2(:,3);
    frameout = table2(:,1);
    K=kfunction(table2(:,2:3),rs,dims,1);
    L=sqrt(K/pi);
    H=L-rs';
%     %dHdr=diff(H)/step;
%     figure
%     scatter(table2(:,2),table2(:,3),5,table2(:,1),'.')
    %axis equal
    %figure
    %plot([0.005:0.01:1.995]',dHdr)
    %         randomPointsX=rand(numberOfMolecules,1);
    %         randomPointsY=rand(numberOfMolecules,1);
    %         simulatedMolecules = [randomPointsX*(xmax-xmin)+xmin,randomPointsY*(ymax-ymin)+ymin];
    %         figure(2)
    %         plot(simulatedMolecules(:,1), simulatedMolecules(:,2), '.', 'MarkerSize',5)
    %         K_sim=kfunction(simulatedMolecules, rs,[xmin,xmax,ymin,ymax]);
    %         L_sim=sqrt(K_sim/pi);
    %         H_sim=L_sim-rs';
    %figure(1)
    [rsmax,rsargmax]=max(H);
    Rmax=rs(rsargmax)
%     figure
%     plot(rs',H)
%     xlabel('r (nm)')
%     ylabel('L(r)-r (nm)')

end
    