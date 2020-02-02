%%Organizes min6 traces from excel file from ImageAnalyzer

function data = phase_shift(conglom_data,interest)

    num_cells = round((size(conglom_data,2)-1)/7); %num of cells - assuming channels in IA
    num_timepoints = size(conglom_data,1); %num of time points

    data = struct;
    data.time = conglom_data(:,1);
    data.cfp_direct = conglom_data(:,2:num_cells+1);
    data.yfpfret = conglom_data(:,num_cells+2:2*num_cells+1);
    data.yfp_direct = conglom_data(:,2*num_cells+2:3*num_cells+1);
    data.FRET = conglom_data(:,3*num_cells+2:4*num_cells+1);
    data.FRET_norm = conglom_data(:,4*num_cells+2:5*num_cells+1);
    data.rfp = conglom_data(:,5*num_cells+2:6*num_cells+1);
    data.rfp_norm = conglom_data(:,6*num_cells+2:end);
    
    phase = [];
    strength = [];
    regularity = [];

    %rfp_lpfilt = data.rfp_norm(:,interest) - filter((1/10)*ones(1,10),1,data.rfp_norm(:,interest));
    %rfp_corrected = data.rfp_norm - repmat(mean(data.rfp_norm(1:8,:)),size(data.rfp_norm,1),1);
    %fret_corrected = data.FRET_norm - repmat(mean(data.FRET_norm(1:8,:)),size(data.rfp_norm,1),1);
    
    %high-pass filter calcium and fret
    rfp_corrected = data.rfp_norm - reshape(smooth(data.rfp_norm,20),size(data.rfp_norm));
    fret_corrected = data.FRET_norm - reshape(smooth(data.FRET_norm,20),size(data.rfp_norm));
    
    %normalize calcium and fret so max = 1
    rfp_corrected = rfp_corrected./repmat(max(rfp_corrected),size(rfp_corrected,1),1);
    fret_corrected = fret_corrected./repmat(max(fret_corrected),size(fret_corrected,1),1);
    
    %just needed to put quest dialog
    figure
    fprintf('Move figure away from the center of the screen then press enter... \n')
    pause
        
    %loop through each cell
    for i=1:size(rfp_corrected,2)
        
        %calculate cross corr between calcium and fret
        [acor,lag] = xcorr(rfp_corrected(6:end,i),fret_corrected(6:end,i)); %calculate autocorrelation between rfp sig and fret sig
        
        %redefine lag and acor so that we only look at 40 time-frame window
        lag = lag(ceil(end/2)-20:ceil(end/2)+20);
        acor = acor(ceil(end/2)-20:ceil(end/2)+20);
        
        %peak pick - change prom if necessary
        [pks,locs] = findpeaks(acor/max(acor),lag,'MinPeakProminence',0.7); % find the peaks of the xcorr
        
        disp(i)
        plotyy(data.time,data.FRET_norm(:,i),data.time,data.rfp_norm(:,i));
        
%         figure;
%         subplot(1,3,1);
%         plotyy(data.time,data.FRET_norm(:,i),data.time,data.rfp_norm(:,i));
%         subplot(1,3,2);
%         plotyy(data.time,fret_corrected(:,i),data.time,rfp_corrected(:,i));
%         subplot(1,3,3);
%         plot(lag,acor);
%         findpeaks(acor/max(acor),lag,'MinPeakProminence',0.7);
        
        %locs
        if ~isempty(min(abs(locs)))
            phase(1,i) = min(abs(locs));
            strength(1,i) = max(acor);
        else
            phase(1,i) = NaN;
            strength(1,i) =NaN;
        end
        
        questStr = 'What is this cells phase relationship?';
        choice = questdlg(questStr, 'Choose One','In-phase','Out-of-phase','Inconcl','Inconcl');
    
        switch choice
            case 'In-phase'   
                phase(2,i) = 1;
            case 'Out-of-phase'
                phase(2,i) = 2;
            case 'Inconcl'
                phase(2,i) = 3;
        end
        
        questStr = 'Is this cell responding (calcium)?';
        choice = questdlg(questStr, 'Choose One','Regular','Irregular','No','No');
    
        switch choice
            case 'Regular'
                regularity(1,i) = 1;
            case 'Irregular'
                regularity(1,i) = 2;
            case 'No'
                regularity(1,i) = 3;
        end
        
%         figure;
%         plotyy(data.time,fret_corrected(:,i),data.time,rfp_corrected(:,i));
        
    end
    
    data.phase = phase;
    data.strength = strength;
    data.regularity = regularity;
    
    close all
end