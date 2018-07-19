%
% read NWB file of one new/old session
%
% this NWB file is created using no2nwb.m, which exports the released data from a proprietary format to NWB
%
% plots raster of a single cell to show category tuning during recognition (regardless of correct/incorrect)
%
%
%urut/060618

function no2nwb_readResult(fNameIn)

%% == import
nwb_in = nwbRead( fNameIn );

%% == get the data (here only the recog part is used)
sessionID = [ 'NWB-' nwb_in.general.session_id{1} ];   %which session was this experiment done with

% == loop over epochs, to create periods variable
periodsRecog=[];
for k=1:100
    epochData = nwb_in.epochs.(sprintf(['trial%03.0f_recog'], k));    
    periodsRecog(k,:) = [ k epochData.start_time*1e6 - 1e6 epochData.start_time*1e6 + 2*1e6];   % start 1s before stim onset, till 2s after stim onset 
    
    %parse attributes string
    tags = epochData.tags;    
    ind=strfind(tags{1},'Cat');    
    catOfTrial = str2num(tags{1}(ind+4));
    
    categoryOfStimulus(k) = catOfTrial;
    
    %Todo: add here to go from the epoch into nwbFile.acquisition.timeseries to determine what the response and RT was
end

%=== neurons
% which channel/cell to plot
channelNr = 2;
cellNr = 2;

clustering_ofChannel = nwb_in.processing.clustering.(sprintf('A%02.0f_cells', channelNr));

% create timestampsOfCell equivalent for this cell
timestampsOfCell = clustering_ofChannel.times( find( clustering_ofChannel.num==cellNr ) )' * 1e6;   % convert from sec to us
brainArea_ofChannel = clustering_ofChannel.description{1};

%% simplified plotting function (raster), showing category selectivity. Uses functions provided in the data release

stimOnset = 1000; %relative to begin of trial,when does the stimulus come on.
stimLength = 1000; 
delayPeriod = 500; %how long till the question comes on
trialLength = stimOnset + stimLength + delayPeriod;
countPeriod = [ stimOnset stimOnset+stimLength ];

indCat1 = find( categoryOfStimulus == 1 );
indCat2 = find( categoryOfStimulus == 2 );
indCat3 = find( categoryOfStimulus == 3 );
indCat4 = find( categoryOfStimulus == 4 );
indCat5 = find( categoryOfStimulus == 5 );

countStimulus_long = extractPeriodCountsSimple( timestampsOfCell, periodsRecog, stimOnset+200, stimOnset+1700, 1 );  % 200 to 1.7s after stim onset

DVs = {   categoryOfStimulus };
[pCategory,table,stats] = anovan( countStimulus_long, DVs,'alpha', 0.05,'display','off', 'model', 'interaction');

groupLabels = {'C1','C2','C3','C4','C5','all'};
binsizePlotting = 250;
alphaLim = 0.05;
normalize = 0;
figNr=1;

subplotSize = plotRasters(binsizePlotting, alphaLim, figNr, ['ID:' sessionID ' p:' num2str(pCategory) ' C' num2str(channelNr) '-' num2str(cellNr) ' ' brainArea_ofChannel], timestampsOfCell, trialLength+500, [stimOnset stimOnset+stimLength], countPeriod, normalize, groupLabels, ...
    periodsRecog(indCat1,:), periodsRecog(indCat2,:), periodsRecog(indCat3,:), periodsRecog(indCat4,:), periodsRecog(indCat5,:), periodsRecog );
