%uses object mask to extract data from processed image stacks
function Experiment = pomansfacs;

    prompt = {'How many channels?','Which slices need deleting?',...
        'Which slice is normalization?','Which channel oscillates?','Inverse probe?','FRET Peak Prom?'};
    answer = inputdlg(prompt);
    
    bad_slices = str2num(answer{2});
    norm_slice = str2num(answer{3});
    osc_channel = str2num(answer{4});
    inv_probe = str2num(answer{5});
    fret_peak_prom = str2num(answer{6});

    %central data-containing structure
    clear Experiment;
    Experiment = struct;
    
    %default folder (hard drive)
    h = 'E:\BrianT\BT_UCSD_Imaging\';
    if ~exist(h,'dir')
        h = pwd;
    end

    %find and extract names of all images to analyze
    folder = uigetdir(h,'Select folder containing multi-channel stack');
    all_images = dir(folder);
    
    %step through all images in the folder, save them locally
    k = waitbar(0,'Importing images...');
    for j = 1:numel(all_images)
        
        if ~all_images(j).isdir
            
            info = regexp(all_images(j).name,expression,'names');
            if(~isfield(info,'channel')) info.channel='1'; end;
            
            Experiment.Images{str2num(info.channel),str2num(info.frame)} = imread(strcat(folder,'\',all_images(j).name));
            
        end
        
        waitbar(j/numel(all_images),k);
        
    end
    close(k)

    %extract time from .inf file
    [time_file,time_file_path] = uigetfile(strcat(h,'\*.inf'),'Select data file');
    fdat = fopen(strcat(time_file_path,'\',time_file));
    d = textscan(fdat, '%f %c %f %*[^\n]','CommentStyle','*');
    fclose(fdat);
    Experiment.time = d{1,3}(:,1)'/6000;
    
    %delete bad slices
    for i = 1:length(bad_slices)
        Experiment.time = [Experiment.time(1,1:bad_slices(i)-i) Experiment.time(1,(bad_slices(i)-i+2):end)];
    end
    
    %find object mask from Cell Profiler
    [mask_file,mask_path] = uigetfile(strcat(h,'\*.tif'),'Select object mask');
    Experiment.obj_mask = imread(strcat(mask_path,'\',mask_file));
    
    %step through each channel
    for channel = 1:size(Experiment.Images,1) 
        
        k = waitbar(0,strcat('Measuring objects in Channel ',int2str(channel),'...'));
    
        %use object mask to find avg intensity for every object through
        %time
        traces=[];
        for object = 1:max(max(Experiment.obj_mask))
            
            bin_mask = Experiment.obj_mask==object;
            for slice = 1:length(Experiment.Images)
                
                temp_cell = Experiment.Images{channel,slice}(bin_mask);
                traces(object,slice) = mean(temp_cell);
                
            end
            
            waitbar(double(object)/double(max(max(Experiment.obj_mask))),k);
        end
        
        close(k)
        Experiment.Data{channel,1} = traces;
    
    end
    
    %calculate fret and norm signals for each object
    for k = 1:max(max(Experiment.obj_mask))
        
        %"oscillating" channel (used as ref) - rfp if using rcamp
        Experiment.AnalyzeData.OscChannel(k,:) = Experiment.Data{osc_channel,1}(k,:)...
            ./repmat(Experiment.Data{osc_channel,1}(k,norm_slice),1,length(Experiment.Data{osc_channel,1}(k,:)));

        if str2num(answer{1})>1
            Experiment.AnalyzeData.yfpfret(k,:) = Experiment.Data{1,1}(k,:);
            Experiment.AnalyzeData.cfp(k,:) = Experiment.Data{2,1}(k,:);
        
            if ~inv_probe
                Experiment.AnalyzeData.nonnorm_fret(k,:) = Experiment.AnalyzeData.yfpfret(k,:)./Experiment.AnalyzeData.cfp(k,:);
            else
                Experiment.AnalyzeData.nonnorm_fret(k,:) = Experiment.AnalyzeData.cfp(k,:)./Experiment.AnalyzeData.yfpfret(k,:);
            end
 
            Experiment.AnalyzeData.norm_fret(k,:) = Experiment.AnalyzeData.nonnorm_fret(k,:)...
                ./repmat(Experiment.AnalyzeData.nonnorm_fret(k,norm_slice),1,length(Experiment.AnalyzeData.nonnorm_fret(k,:)));
        end

    end
  
    %apply periodicity filter to find TEA-responsive cells
    %calcium - raw intensity filtered
    resp_cells = periodicity_filter(Experiment.Data{osc_channel,1},Experiment.time,0);
    Experiment.resp_cells.OscChannel = resp_cells;
    
    %plot stuff
    figure
    plot(Experiment.time,Experiment.AnalyzeData.OscChannel(Experiment.resp_cells.OscChannel(:,1),:))
    
    if str2num(answer{1})>1
        
        resp_cells_other = periodicity_filter(Experiment.AnalyzeData.norm_fret,Experiment.time,fret_peak_prom);
        Experiment.resp_cells.Other = resp_cells_other;
        
        intersec = intersect(Experiment.resp_cells.OscChannel(:,1),Experiment.resp_cells.Other(:,1));
        Experiment.intersec = intersec;
        
        figure
        fprintf('Move figure away from the center of the screen then press enter... \n')
        pause
        filtered_responded = [];
        oscchannel_resp = [];
        fret_resp = [];
        filtered_nonresp = [];
        oscchannel_nonresp = [];
        fret_nonresp = [];
        
        for y = 1:length(intersec)
         
            plotyy(Experiment.time,Experiment.AnalyzeData.OscChannel(intersec(y),:),Experiment.time,Experiment.AnalyzeData.norm_fret(intersec(y),:));
            
            questStr = 'Is this cell responding (non-calcium)?';
            choice = questdlg(questStr, 'Yes','No');
        
            switch choice
                case 'Yes'
                    filtered_responded = [filtered_responded; intersec(y)];
                    oscchannel_resp = [oscchannel_resp; Experiment.AnalyzeData.OscChannel(intersec(y),:)];
                    fret_resp = [fret_resp; Experiment.AnalyzeData.norm_fret(intersec(y),:)];
%                     figure
%                     pause
                case 'No'
                    filtered_nonresp = [filtered_nonresp; intersec(y)];
                    oscchannel_nonresp = [oscchannel_nonresp; Experiment.AnalyzeData.OscChannel(intersec(y),:)];
                    fret_nonresp = [fret_nonresp; Experiment.AnalyzeData.norm_fret(intersec(y),:)];
                end
        end
        
        Experiment.filtered_resp.filtered_responded_cells = filtered_responded;
        Experiment.filtered_resp.oscchannel = oscchannel_resp;
        Experiment.filtered_resp.fret = fret_resp;
        Experiment.filtered_nonresp.filtered_nonresponded_cells = filtered_nonresp;
        Experiment.filtered_nonresp.oscchannel = oscchannel_nonresp;
        Experiment.filtered_nonresp.fret = fret_nonresp;
        
%         figure
%         plot(Experiment.time,Experiment.AnalyzeData.norm_fret(Experiment.resp_cells.OscChannel(:,1),:))
    end
    
%     corr_cells = correlation_filter(resp_cells,osc_signal,norm_fret);
%     Experiment.corr_cells = corr_cells;
%     
%     for j = 1:size(corr_cells,1)
%         figure
%         plot(Experiment.time,osc_signal(corr_cells(j,1),:))
%         hold on
%         plot(Experiment.time,norm_fret(corr_cells(j,1),:))
%     end
    
    Experiment.Images = [];
end

%periodicity filter to find TEA-responding cells 
function resp_cells = periodicity_filter(analyzed_mat,time,peakprom)

    %compute fft
    freq_dom = abs(fft(analyzed_mat,[],2));
    freq_dom = freq_dom(:,4:end-4); %trimming first four out on both ends of traces
    time = time(1,4:end-4);
    %time = 1:size(freq_dom,2);
    size(freq_dom)
    size(time)
    resp_cells = [];
    for obj = 1:size(freq_dom,1)

        %peak picking
        [pks,loc, width, prom] = findpeaks(freq_dom(obj,:),time(1,:),'MinPeakProminence',peakprom);
        
        %store responding cells
        if ~isempty(loc) && (min(loc)>1 && min(loc)<10)
            [m,i] = min(loc);
            resp_cells = [resp_cells; obj min(loc) prom(i)];
%             plot(time,analyzed_mat(obj,4:end-4))
%             hold on
        end
        
    end
    
end

function corr_cells = correlation_filter(resp_cells,osc_signal,norm_fret)

    corr_cells = [];
    for i = 1:size(osc_signal,1)
        
        q(i,:) = xcorr(osc_signal(i,:),norm_fret(i,:));
        
        if max(q(i,:))>100
            corr_cells = [corr_cells; i resp_cells(i,1)];
        end
        
    end
      
end



