function nwbFile = no2nwb(session, dataDir, stimFilesDir)
% NO2NWB Convert new/old (i.e. no) session to NWB file.
%
%   NWBFILE = NO2NWB(SESSION,DATADIR,STIMFILESDIR) converts the new/old 
%   session SESSION into the NWB file NWBFILE by mapping the session data
%   to appropriate containers defined by the NWB format.
%
%   Example:
%     % Create the new/old session struct.
%     NOsession.session = 'P9HMH_032306'; 
%     NOsession.sessionID = 'P9S1';
%     NOsession.EXPERIMENTIDLearn = 80;
%     NOsession.EXPERIMENTIDRecog = 81; 
%     NOsession.taskDescr = 'NO';   
%     NOsession.variant = 1;            
%     NOsession.blockIDLearn = 1;       
%     NOsession.blockIDRecog = 2;       
%     NOsession.patient_nb = 1;         
%     NOsession.patientsession = 1;     
%     NOsession.diagnosisCode = 1;
%
%     % Define the location of the session data and stimulus files.
%     dataDir = '/home/user/RecogMemory_MTL_release_v2/Data';
%     stimFilesDir = '/home/user/RecogMemory_MTL_release_v2/Code/dataRelease/stimFiles';
%
%     % Execute the converter.
%     nwbFile = no2nwb(NOsession, dataDir, stimFilesDir);
%
%     % The NWB file can now be exported to disk.
%     nwbExport(nwbFile, 'my_session.nwb')

% Create the NWB file container. May want to add 'session_description'.
nwbFile = nwbfile( ...
    'file_create_date', {datestr(now, 30)}, ...
    'identifier', {session.session});

% Create the basic NWB data structure. This is an oddity of MatNWB at the
% moment. Hopefully in the future this is no longer necessary. All NWB
% files share this basic structure.
nwbFile.acquisition.timeseries = types.untyped.Group();
nwbFile.stimulus.presentation = types.untyped.Group();

% Load the session events file. The events in the events file serve as an
% experimental timeline for this session.
eventsFile = fullfile(dataDir, 'events', session.session, session.taskDescr, 'eventsRaw.mat');
data = load(eventsFile);
events = array2table(data.events(:,1:3), 'VariableNames', {'timestamp', 'ttl', 'experimentId'});

% Convert the event timestamps from microseconds to seconds. The NWB format
% says "all times are stored in seconds".
% See: http://nwb-schema.readthedocs.io/en/latest/format.html#storing-time-values
events.timestamp = events.timestamp / 1e6;

% Add experimental metadata to the NWB file.
addGeneralInfo(nwbFile, session);

% Add the stream of acquired TTL values, representing events, to the NWB file.
addAcquisitionEvents(nwbFile, events, [session.EXPERIMENTIDLearn session.EXPERIMENTIDRecog], {'learn','recog'} );

% Filter down the events to only those that are related to the experiments
% of this session. We could have done this before calling
% addAcquisitionEvents but we currently store all events in the file for 
% completeness. The epochs will serve to identify the relevant events.
tf = events.experimentId == session.EXPERIMENTIDRecog | events.experimentId == session.EXPERIMENTIDLearn;
experimentEvents = events(tf, :);

%--- process learning part first
mode=1; %learn
expStr='learn';
[order_learn,mapping,stimIDs_shown_learn] = loadStimInfo(session.variant, stimFilesDir,  session.blockIDLearn,  session.blockIDRecog, mode);
addPresentationImages(nwbFile, [], mapping, stimIDs_shown_learn, expStr, mode);
addEpochs(nwbFile, events(events.experimentId == session.EXPERIMENTIDLearn, :), expStr, mode, mapping, stimIDs_shown_learn, []);


%--- process recognition part second
mode=2; %recog
expStr='recog';
[order_recog,~,stimIDs_shown_recog] = loadStimInfo(session.variant, stimFilesDir,  session.blockIDLearn,  session.blockIDRecog, mode);
addPresentationImages(nwbFile, [], mapping, stimIDs_shown_recog, expStr, mode);
addEpochs(nwbFile, events(events.experimentId == session.EXPERIMENTIDRecog, :), expStr, mode, mapping, stimIDs_shown_learn, stimIDs_shown_recog);

% Add the filenames of images presented during the session, in the order they were presented, to the NWB file.
%addPresentationImages(nwbFile, experimentEvents, stimFilesDir, session.variant, session.blockIDLearn, session.blockIDRecog, 'learn', 1);
%addPresentationImages(nwbFile, experimentEvents, stimFilesDir, session.variant, session.blockIDLearn, session.blockIDRecog, 'recog', 2);

% Add epochs that represent each of the experimental trials performed
% during the experiment to the NWB file.
%addEpochs(nwbFile, events(events.experimentId == session.EXPERIMENTIDRecog, :), 'recog', 2, stimFilesDir, session.variant);

%--- process sorted clusters

% Add the cluster spike data and mean waveforms as intermediate analysis data (i.e. processing modules) to the NWB file.
addClusteringModules(nwbFile, dataDir, session.session, session.taskDescr);
end

function addGeneralInfo(nwbFile, session)
% Currently only adding the experimental metadata available to us through
% the session data. However, information from the published paper could be
% used to fill in many of the other general fields.
% See: http://nwb-schema.readthedocs.io/en/latest/format.html#groups-general
s.datasets.session_id = {session.sessionID};
%s.datasets.experiment_description = {'some experiment description here'};
%s.datasets.experimenter = {'experimenter here'};
nwbFile.general = types.untyped.Group(s);

% The NWB format has an ElectrodeTable container which would be an
% appropriate place to put the data in electrodePos.xlsx. However,
% ElectrodeTable is not yet supported by MatNWB.
% See: http://nwb-schema.readthedocs.io/en/latest/format.html#sec-electrodetable
%nwb.general.extracellular_ephys = types.untyped.Group();
%nwb.general.extracellular_ephys.electrodes = types.ElectrodeTable()
end

function addAcquisitionEvents(nwbFile, events, expIDs_toExport, expIDs_taskDescr)
% Most of the events table fits well into a basic TimeSeries container 
% under the NWB file 'acquisition' group. However, in the NOsession each 
% event is associated with an experimentId. We've lost the experimentIds
% here. I do not currently know of a good place to put them.


for k=1:length(expIDs_toExport)

    indsUse = find( events.experimentId == expIDs_toExport(k) );
    ts = types.TimeSeries( 'data', events.ttl(indsUse), 'timestamps', events.timestamp(indsUse), 'comment', {['EXPID=' num2str(expIDs_toExport(k))]} );

    nwbFile.acquisition.timeseries.(['Events_' expIDs_taskDescr{k}]) = ts;

end

% We're currently not mapping the experiment ids associated with each
% event. Where should these go?
end

function [order,mapping,stimIDs_shown] = loadStimInfo(variant, filesDir, learnBlock, recogBlock, mode)
% Load the two files that define the order of images shown during the
% session and the actual filenames of each of these images.
if variant == 1
    variant = [];
end

order = load(fullfile(filesDir, ['NewOldDelay' num2str(variant) '_v3.mat']));
mapping = load(fullfile(filesDir, ['newOldDelayStimuli' num2str(variant) '.mat']));

% Map from the order, which is defined in stimulusIds, to the actual
% filenames of each of the images.

if mode==1
    stimIDs_shown = order.experimentStimuli(learnBlock).stimuliLearn;
else
    stimIDs_shown = order.experimentStimuli(recogBlock).stimuliRecog;
end

end

%mode: 1 add learn, 2 add recog
function addPresentationImages(nwbFile, events, mapping, stimIDs_shown, expStr, mode)


files = mapping.fileMapping( stimIDs_shown );


%files = mapping.fileMapping;   % list of all filenames, in original order so they can be referenced.
%categoryID_forEachStim = mapping.categoryMapping(:,2);

% Extract the timestamps for just the "stimulus on" events. These
% timestamps indicate the time that each of these images were displayed.
%ttls = markerIds();
%tf = events.ttl == ttls.STIMULUS_ON;
%timestamps = events(tf, :).timestamp;

% Use an ImageSeries container to store timestamps and filenames of each
% image presented. We could also use two ImageSeries, one for the learning 
% phase and one for the recognization phase.
% See: http://nwb-schema.readthedocs.io/en/latest/format.html#sec-imageseries
is = types.ImageSeries( ...
    'external_file', files, ...
    'format', {'external'} );  
    
%'timestamps', timestamps);   %UR: seem redundant, not sure why to store the timestamps again

nwbFile.stimulus.presentation.(['images_' expStr]) = is;

% We're currently not mapping the data about each image's category (e.g.
% house, landscape, animal, etc.) or stimulus id. It would be nice if the
% ImageSeries allowed for this. Do we need to create a custom extension or 
% a separate TimeSeries for this data?

% Also, the files are currently defined by absolute paths on the
% experimentation machine. These should probably be changed to relative
% paths to the location of the NWB file so these files can actually be
% found.

nwbFile.stimulus.presentation.(['images_' expStr]) = is;
end

function addEpochs(nwbFile, events, expStr, mode, mapping, stimIDs_shown_learn, stimIDs_shown_recog)

categoryID_forEachStim = mapping.categoryMapping(:,2);

%== determine the subset of stimulus info relevant to the epochs being processed (learn or recog)
newold_ofTrial=[];
if mode==1
    %learn
    stimIDs_toUse = stimIDs_shown_learn;
    
    newold_ofTrial = zeros(length(stimIDs_toUse),1);  % they are all new
else
    %recog
    stimIDs_toUse = stimIDs_shown_recog;
        
    % determine which of these are old/new
    for k=1:length(stimIDs_shown_recog)       
        if ~isempty(find( stimIDs_shown_recog(k) == stimIDs_shown_learn ))
            newold_ofTrial(k) = 1;   % isOld=1
        else
            newold_ofTrial(k) = 0;   % isOld=0     
        end
    end
end


% For each interval between a "stimulus on" event and a "delay2 off" event, 
% we'll call this an experimental trial and create an epoch container.
ttls = markerIds();
indices = find(events.ttl == ttls.STIMULUS_ON);
for i = 1:numel(indices)
    ind = indices(i);
    
    % Lets do some sanity checks.
    assert(events.ttl(ind+1) == ttls.STIMULUS_OFF, ...
        'Expected STIMULUS_OFF TTL at event timestamp %f', events.timestamp(ind+1));
    assert(events.ttl(ind+2) == ttls.DELAY1_OFF, ...
        'Expected DELAY1_OFF TTL at event timestamp %f', events.timestamp(ind+2));
    % ind+3 should depend on the participant's response.
    assert(events.ttl(ind+4) == ttls.DELAY2_OFF, ...
        'Expected DELAY2_OFF TTL at event timestamp %f', events.timestamp(ind+3));
    
    % This EpochTimeSeries container will serve as a window into the events
    % TimeSeries container we already stored under the "acquisition" group.
    eventTs = types.EpochTimeSeries( ...
        'count', 5, ...
        'idx_start', ind-1, ...
        'timeseries', types.untyped.Link( ['/acquisition/timeseries/Events_' expStr]));
    
    % This EpochTimeSeries container will serve as a window into the
    % stimulus TimeSeries container we already stored under the "stimulus"
    % group.
    stimulusTs = types.EpochTimeSeries( ...
        'count', 1, ...
        'idx_start', i-1, ...
        'timeseries', types.untyped.Link( ['/stimulus/presentation/images_' expStr] ));
    
    
    stimID_ofTrial = stimIDs_toUse(i);
    catID_ofTrial = categoryID_forEachStim(stimID_ofTrial);
    
    % This is the actual Epoch container that we define to start at the
    % timestamp of the "stimulus on" event and stop at the timestamp of the
    % "delay2 off" event.
    % See: http://nwb-schema.readthedocs.io/en/latest/format.html#sec-epoch
    %
    %http://nwb-schema.readthedocs.io/en/latest/format.html#sec-eventtable
    epoch = types.Epoch( ...
        'start_time', events.timestamp(ind), ...
        'stop_time', events.timestamp(ind+4), ...     
        'tags', {    ['StimID=' sprintf('%03.0f', stimID_ofTrial) ';Cat=' num2str(catID_ofTrial) ';NO=' num2str(newold_ofTrial(i))]}, ...    
        'groups', struct('events', eventTs, 'stimulus', stimulusTs));       
    
    %UR: we store the category ID and the new/old status of each image here, because it is unclear to us how to store metadata in an ImageSeries    
    
    % The trials are labeled in sequential order with leading zeros so they are nicely sorted when viewed in HDFView.
    nwbFile.epochs.(sprintf(['trial%03.0f_' expStr], i)) = epoch;
end
end

function addClusteringModules(nwbFile, dataDir, sessionName, taskDescr)
% Load the file that defines all of the channel numbers for which at least
% one putative single neuron was identified.
brainAreaFile = fullfile(dataDir, 'events', sessionName, taskDescr, 'brainArea.mat');
data = load(brainAreaFile);
brainArea = array2table(data.brainArea(:,1:4), 'VariableNames', {'channelNr', 'clusterId', 'clusterIdOrig', 'brainAreaId'});

% Filter for cluster id 0. These appear to be invalid clusters?
tf = brainArea.clusterId ~= 0;
brainArea = brainArea(tf, :);

% For each of the channel numbers, place the clustered spike data and mean
% waveforms in appropriate containers under the "processing" group.
% See: http://nwb-schema.readthedocs.io/en/latest/format.html#groups-processing
channelNrs = unique(brainArea.channelNr);
for i = 1:numel(channelNrs)
    % Load the file containing the clustering data for this channel.
    chanNr = channelNrs(i);
    chanFile = fullfile(dataDir, 'sorted', sessionName, taskDescr, ['A' num2str(chanNr) '_cells.mat']);
    if ~exist(chanFile, 'file')
        warning('%s indicates there should be a channel file %s but it does not exist. Skipping this channel.', brainAreaFile, chanFile);
        continue;
    end
    chanData = load(chanFile);
    spikes = array2table(chanData.spikes, 'VariableNames', {'clusterId', 'clusterIdOrig', 'timestamp', 'unused', 'spikeIdOrig'});
    
    % Convert all timestamps from microseconds to seconds.
    spikes.timestamp = spikes.timestamp / 1e6;
    
    %brain area of this channel
    brainAreaCode = brainArea.brainAreaId( find( brainArea.channelNr==chanNr ) );
    
    % Store the clustered spike data in a Clustering container.
    % See: http://nwb-schema.readthedocs.io/en/latest/format.html#clustering
    cl = types.Clustering( ...
        'num', spikes.clusterId, ...
        'times', spikes.timestamp, ...
        'description', {['BrainArea=' num2str(brainAreaCode(1)) ]} );
    clustering.(sprintf('A%02.0f_cells', chanNr)) = cl;
    
    % Store the two mean waveforms for the spikes in separate
    % ClusterWaveforms containers. These containers link to the Clustering
    % container above to keep reference to the data that was used to create
    % them. More datasets ('waveform_sd', 'waveform_filtering') should be
    % added to these containers.
    % See: http://nwb-schema.readthedocs.io/en/latest/format.html#clusterwaveforms
    mw = chanData.meanWaveform_learn;
    mw = mw(~cellfun(@(a)isscalar(a), {mw.m_waveform})); % filter scalar waves
    cw = types.ClusterWaveforms( ...
        'waveform_mean', vertcat(mw.m_waveform)', ...
        'clustering_interface', types.untyped.Link(sprintf('/processing/clustering/A%02.0f_cells', chanNr)));
    clusterWaveform.(sprintf('A%02.0f_waves_learn', chanNr)) = cw;
    
    mw = chanData.meanWaveform_recog;
    mw = mw(~cellfun(@(a)isscalar(a), {mw.m_waveform})); % filter scalar waves
    cw = types.ClusterWaveforms( ...
        'waveform_mean', vertcat(mw.m_waveform)', ...
        'clustering_interface', types.untyped.Link(sprintf('/processing/clustering/A%02.0f_cells', chanNr)));
    clusterWaveform.(sprintf('A%02.0f_waves_recog', chanNr)) = cw;
    
    % The IsolDist and SNR (signal-to-noise) data is not currently being
    % mapped because I am not sure about the most appropriate place to put
    % these. Should they have their own ProcessingModules?
end
nwbFile.processing.clustering = types.ProcessingModule('groups', clustering);
nwbFile.processing.clustering_waveforms = types.ProcessingModule('groups', clusterWaveform);
end

function ttls = markerIds()
% These TTL identifiers fit better into an enumeration but a struct is used
% here just to keep all the no2nwb code in one file.
% See: https://www.mathworks.com/help/matlab/matlab_oop/enumerations.html
ttls.STIMULUS_ON = 1;
ttls.STIMULUS_OFF = 2;
ttls.DELAY1_OFF = 3;
ttls.DELAY2_OFF = 6;
ttls.EXPERIMENT_ON = 55;
ttls.EXPERIMENT_OFF = 66;
end