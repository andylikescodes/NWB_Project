
installMode = 0;   

if installMode
    %% Get MatNWB
    % Commit of MatNWB to download.
    commit = '824a82d88ff179b5fc6eedc62a860cd527bd8683';
    
    % Download the MatNWB code from GitHub as a zip archive.
    websave('matnwb-archive.zip', ['https://github.com/NeurodataWithoutBorders/matnwb/archive/' commit '.zip']);
    
    % Unzip the MatNWB code archive.
    unzip('matnwb-archive.zip', 'matnwb-archive');
    movefile(fullfile('matnwb-archive', ['matnwb-' commit]), 'matnwb');
    rmdir('matnwb-archive');
    
    % Add the MatNWB functions to the MATLAB path.
    addpath('matnwb');
    
    % Generate the core MatNWB types from the NWB schema.
    cd('matnwb');
    generateCore(fullfile('schema','core','nwb.namespace.yaml'));
    cd('..');
    
    %% Get new/old single-neuron MTL recordings data
    
    % Download the NO data from Dryad as a zip archive.
    websave('RecogMemory_MTL_release_v2.zip', 'https://datadryad.org/bitstream/handle/10255/dryad.163179/RecogMemory_MTL_release_v2.zip');
    
    % Unzip the NO data archive.
    unzip('RecogMemory_MTL_release_v2.zip', 'RecogMemory_MTL_release_v2');
    
    % Add dataRelease folder to path for defineNOsessions_release function.
    addpath(fullfile('RecogMemory_MTL_release_v2', 'Code', 'dataRelease'));
end

%% Run no2nwb
% Set data and stim files paths.


%adjust to where the data was downloaded to
dataPath = './RecogMemory_MTL_release_v2/Data';
stimFilesPath = './RecogMemory_MTL_release_v2/Code/dataRelease/stimFiles';

addpath(dataPath);
addpath(stimFilesPath);
addpath('./RecogMemory_MTL_release_v2/Code/dataRelease/helpers');
addpath('./RecogMemory_MTL_release_v2/Code/dataRelease/analysisBehavior');
addpath('./RecogMemory_MTL_release_v2/Code/dataRelease/newoldtask');
% Get session definitions.
[sessions, use] = defineNOsessions_release();

% Map the first usable session.
s1 = sessions(use(1));
file = no2nwb(s1, dataPath, stimFilesPath);

% Export the mapped NWBFile to disk.
fName_NWB = [s1.session '.h5'];
nwbExport(file, fName_NWB);

%% read back the info from the NWB file to plot
no2nwb_readResult( fName_NWB );