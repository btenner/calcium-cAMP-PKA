%rs are all the values of r that you will be running the ripley analysis
%for
rs = 1:1:1000;

% strings={%'grouped_140_100_pPKAa647_AKAP79a568_sap_HeLa_FI_cell004_A647'
%     %'grouped_140_100_pPKAa647_AKAP79a568_sap_HeLa_FI_cell005_A647'}
%     %'grouped_140_100_pPKAa647_AKAP79a568_sap_HeLaHG_FI_cell007_A647'
%     %'grouped_140_100_pPKAa647_AKAP79a568_sap_HeLaHG_FI_cell008_A647'
%     %'grouped_140_100_pPKAa647_AKAP79a568_sap_HeLaHG_FI_cell009_A647'
%     %'grouped_45_100_pPKAa647_AKAP79a568_sap_HeLa_FI_cell004_A568' 
%     %'grouped_45_100_pPKAa647_AKAP79a568_sap_HeLa_FI_cell005_A568'
%     %'grouped_45_100_pPKAa647_AKAP79a568_sap_HeLaHG_FI_cell007_A568'
%     %'grouped_45_100_pPKAa647_AKAP79a568_sap_HeLaHG_FI_cell008_A568'
%     %'grouped_45_100_pPKAa647_AKAP79a568_sap_HeLaHG_FI_cell009_A568'};
% 
strings =  {
    'grouped_160_100_Min6_AC8NtermEGFP_AC8_647_001.txt'};

% dimsMatrix= [15e3, 25e3, 15e3, 25e3;
%     10e3, 15e3, 20e3, 25e3;
%     14e3, 19e3, 22e3, 27e3;
%     20e3, 25e3, 22e3, 27e3;
%     21e3, 26e3, 18e3, 23e3;
%     18e3, 23e3, 18e3, 23e3;
%     10e3, 15e3, 20e3, 25e3;
%     14e3, 19e3, 22e3, 27e3;
%     20e3, 25e3, 22e3, 27e3;
%     21e3, 26e3, 18e3, 23e3];
% 


% dimsMatrix= [18e3, 28e3, 18e3, 28e3;
%     10e3, 20e3, 20e3, 30e3;
%     14e3, 24e3, 22e3, 32e3;
%     20e3, 30e3, 22e3, 32e3;
%     21e3, 31e3, 18e3, 28e3;
%     18e3, 28e3, 18e3, 28e3;
%     10e3, 20e3, 20e3, 30e3;
%     14e3, 24e3, 22e3, 32e3;
%     20e3, 30e3, 22e3, 32e3;
%     21e3, 31e3, 18e3, 28e3];


% FI cell1 = [1.5e4 2.5e4 1.5e4 2.5e4]
% FI cell3 = [1.5e4 2.5e4 2e4 3e4]
% FI cell4 = [1.5e4 2.5e4 1.5e4 2.5e4]
% 
% H89 cell1 = [1.5e4 2.5e4 2e4 3e4]
% H89 cell2 = [1.5e4 2.5e4 1.5e4 2.5e4]
% H89 cell3 = [1e4 2e4 1.5e4 2.5e4]
% H89 cell4 = [1e4 2e4 2e4 3e4]
% 
% nopre cell1 = [1.5e4 2.5e4 1.5e4 2.5e4]
% nopre cell2 = [2e4 3e4 1.5e4 2.5e4]
% nopre cell3 = [1.5e4 2.5e4 1.5e4 2.5e4]

dimsMatrix = [
    2.2e4 2.7e4 1.2e4 1.7e4;
    
];

%[1.5e4 2.5e4 1.5e4 2.5e4;
%2e4 2.5e4 2e4 2.5e4;
%1.5e4 2.5e4 1.5e4 2.5e4;
%1.5e4 2.5e4 2e4 3e4;
%1.5e4 2.5e4 1.5e4 2.5e4;
%1e4 2e4 1.5e4 2.5e4;
%1e4 2e4 2e4 3e4;
%1.5e4 2.5e4 1.5e4 2.5e4;
%2e4 3e4 1.5e4 2.5e4;
%1.5e4 2.5e4 1.5e4 2.5e4];

shortnames = { 'AC8cell31'};
%shortnames = {'cell4_pPKA','cell5_pPKA','cell7_pPKA','cell8_pPKA','cell9_pPKA', 'cell4_AKAP79','cell5_AKAP79','cell7_AKAP79','cell8_AKAP79','cell9_AKAP79'};

graphTitles = {'AC8 cell 3'};
%graphTitles = {'Cell 4 pPKA', 'Cell 5 pPKA','Cell 7 pPKA', 'Cell 8 pPKA','Cell 9 pPKA', 'Cell 4 AKAP79', 'Cell 5 AKAP79', 'Cell 7 AKAP79', 'Cell 8 AKAP79','Cell 9 AKAP79' };


for i=1:length(strings)
    stringName = strings{i}
    varName = eval(stringName);
    dims = dimsMatrix(i,:);
    xstring = strcat('x',shortnames{i});
    ystring = strcat('y',shortnames{i});
    fstring = strcat('f',shortnames{i});
    Hstring = strcat('H',shortnames{i});
    Rstring = strcat('R',shortnames{i});
    %eval(strcat('[', xstring , ',', ystring, ',', fstring, ',', Hstring, ',' , Rstring, ']=ripley_analysis_STORM_lowmem(varName,dims,rs);'));
    eval(strcat('[', xstring , ',', ystring, ',', fstring, ',', Hstring, ',' , Rstring, ']=ripley_analysis_STORM(varName,dims,rs);'));

    figure
    eval(strcat('scatter(',xstring,',',ystring,',4,', fstring,',', char(39),'filled',char(39),')'))
    axis equal tight
    title(graphTitles{i});
    figure(i+10)
    eval(strcat('plot(rs,', Hstring, ')'))
    xlabel('r (nm)');
    ylabel('L(r)-r (nm)');
    title(graphTitles{i});
    axis square

    
end

% plot(rs,HH891,rs', HH892, rs, HH893,rs,HH894);
% title('Fsk/IBMX');
% xlabel('r (nm)');
% ylabel('L(r)-r (nm)');
% title('H89');
% figure
% plot(rs,Hnopre1,rs, Hnopre2, rs, Hnopre3)
% xlabel('r (nm)');
% ylabel('L(r)-r (nm)');
% title('no pretreatment');