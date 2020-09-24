%% Plotting a dipole with moving orientation
%in this script, a dipole with fixed position is displayed with altered
%orientations
%Note: The Cortex display will take a while to load. You need to go to the
%GUI and unclick "show all time" on the right side. I have yet to figure
%out how to do that in the script. After that, you will have all 4 displays
%of the same activity. To make a movie out of them, arrange their windows
%the way you want. Then right-click one of them and go to 
%Snapshot --> Movie: all figures --> done!

%%%%%%%%%%%%SETTINGS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SubjectName = 'Subject26'; %specify the name of the subject that will be used in the brainstorm GUI
SubjectName = 'Subject03'; %specify the name of the subject that will be used in the brainstorm GUI
MAINPATH = 'C:\'; %path to brainstorm3
%PROTOCOLNAME = 'Ana'; %name of your current protocol in brainstorm
PROTOCOLNAME = 'Protocol01'; %name of your current protocol in brainstorm

%%%%%%%%%%%%SETTINGS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Copy existing folder
%this is just to avoid confusion, it gets messy with too many files in one
%folder, so we work in a copy of the original raw-folder (should only have the files from the first script)
listFolder = dir([MAINPATH,'brainstorm_db\',PROTOCOLNAME,'\data\',SubjectName,'\raw\results*']); %lists of the variables created in PrepareSource_cEEGrid.mat

% Input files
sFiles = {...
    [SubjectName,'/raw/',listFolder(1).name], ...
    [SubjectName,'/raw/',listFolder(2).name]};

% Start a new report
bst_report('Start', sFiles);

% Process: Duplicate folders: Add tag "_plots"
sFiles = bst_process('CallProcess', 'process_duplicate', sFiles, [], ...
    'target', 2, ...  % Duplicate folders
    'tag',    '_plots');

%% Load the channel file and the dipole file
listNewChannel = dir([MAINPATH,'brainstorm_db\',PROTOCOLNAME,'\data\',SubjectName,'\raw_plots*\channel*']);
export_matlab([listNewChannel(end).folder,'\',listNewChannel(end).name],'ReadyCap');

listSource = dir([MAINPATH,'brainstorm_db\',PROTOCOLNAME,'\data\',SubjectName,'\raw_plots*\results_MN*']); %lists of the variables created in PrepareSource_cEEGrid.mat
listChannel = dir([MAINPATH,'brainstorm_db\',PROTOCOLNAME,'\data\',SubjectName,'\raw_plots*\channel*']); %raw s.t. starts with R or r
listHeadmodel = dir([MAINPATH,'brainstorm_db\',PROTOCOLNAME,'\data\',SubjectName,'\raw_plots*\headmodel*']);

Source = load([listSource(1).folder,'\',listSource(1).name]); %load the source template that was created in PrepareSource_cEEGrid...
CM = load([listChannel(end).folder,'\',listChannel(end).name]); %...and the channelfile
HM = load([listHeadmodel(1).folder,'\',listHeadmodel(1).name]); %...and the headmodel

%% store new locations and orientations of the dipoles

NormOrient = [-0.0434816861376882,0.985049224918913,-0.103375097755849;0.993581456236887,0.0542738434228215,0.0992483740033164;-0.104270505035261,0.0992483740033164,0.989584469379643]; %dip.Dipole(1).Amplitude; %This is the orientation of the dipole

listDip = dir([MAINPATH,'brainstorm_db\',PROTOCOLNAME,'\data\',SubjectName,'\raw_plots*\dipoles*']);
export_matlab([listDip(end).folder,'\',listDip(end).name],'dip'); %load the dipole file from PrepareSourceDipole.m

dip.Dipole(1).Loc = HM.GridLoc(143,:)'; %3 dipoles with the same location and differ only in orientation
dip.Dipole(2).Loc = HM.GridLoc(143,:)'; % 143 is a random location from the GridLoc, it can be altered to be any position
dip.Dipole(3).Loc = HM.GridLoc(143,:)';

dip.Dipole(1).Amplitude = NormOrient(1,:)'; %here, the orientations for the loc are defined
dip.Dipole(2).Amplitude = NormOrient(2,:)';
dip.Dipole(3).Amplitude = NormOrient(3,:)';

dip.Dipole = dip.Dipole(1:3);

vecTime = linspace(0,1,length(dip.Dipole));
HM.GridLoc = round(HM.GridLoc(:,:),4); % this line is necessary because matlab rounds in a funny way

for w = 1:length(dip.Dipole)
    dip.Dipole(w).Loc = round(dip.Dipole(w).Loc,4); % this line is necessary because matlab rounds in a funny way
    dip.Dipole(w).Time = vecTime(w);
    dip.Dipole(w).Index = 1;
    dip.Dipole(w).Origin = [0,0,0];
    dip.Dipole(w).Goodness = 1;
end
dip.Time = [dip.Dipole.Time];

w=1;
for u = 1:length(dip.Dipole)
    [~,findHM(w,:)] = ismember(dip.Dipole(u).Loc',HM.GridLoc,'rows'); %this loop finds for each wanted dipole location the equivalent on the HM
    w=w+1; %only really useful when you have other dipole locations specified; in this example, it is always the same
end

%% Some settings for the cEEGrid channels

ProtocolInfo = bst_get('ProtocolInfo'); %needed for number of the current protocol. ALWAYS click the nedded protocol AND subject in the GUI
iChannels = 1:length(CM.Channel); %number of channels according to your channel file

smoothVec = linspace(0,90,50); %create a vector for 50 dipoles, with steps from 0 to 90 (for rotation of source direction of up to 90°)
timeVec = linspace(0,10,length(smoothVec)*2);

%% load in the raw signal file and extend it to be long enough for 50 dipoles to be shown for 100 samples each (=5000 sample points)

%listRaw = dir([MAINPATH,'brainstorm_db\',PROTOCOLNAME,'\data\',SubjectName,'\raw_plots*\data_block*']);
listRaw = dir([MAINPATH,'brainstorm_db\',PROTOCOLNAME,'\data\',SubjectName,'\raw_plots*\*data_sim*']);

export_matlab([listRaw(end).folder,'\',listRaw(end).name],'rawSim');
rawSim.F = rawSim.F(:,1);
rawSim.F = repmat(rawSim.F,1,length(timeVec)*50);
rawSim.Time = linspace(0,10,length(smoothVec)* length(timeVec));
rawSim.ChannelFlag(1:end) = 1;
rawSim.Comment = 'SimDip';

ProtocolInfo = bst_get('ProtocolInfo'); %reload the altered signal into bst
db_add(ProtocolInfo.iStudy, rawSim, 1);
%% get results from different orientations. The point here is to take the
%initial 3 dipole orientations and add dipoles with angles in between the
%orthogonals. This way, we get a "moving" dipole orientation that we can plot
  
dip_bst = dip;
dip_bst.Dipole(1) = dip.Dipole(1);
dip_bst.Time = timeVec;

vec1 = dip.Dipole(1).Amplitude; %1st vec
vec2 = dip.Dipole(3).Amplitude; % z vec

for t = 1 : length(smoothVec)
    angle(t) = smoothVec(t)*pi/180; % angle of rotation in rad
    v(:,t) = vec1 + tan(angle(t))*vec2; % length of vec1 is 1
    v(:,t) = v(:,t)/norm(v(:,t)); % v is a unit vector 5 degrees from vec1
    dip_bst.Dipole(t) = dip_bst.Dipole(1);
    dip_bst.Dipole(t).Amplitude = v(:,t);
    dip_bst.Dipole(t).Time = timeVec(t);
end

vec1 = dip.Dipole(3).Amplitude; % z vec
vec2 = dip.Dipole(2).Amplitude; % x vec

for t = length(smoothVec)+1: length(smoothVec)*2
    angle(t) = smoothVec(t-length(smoothVec))*pi/180;   % angle of rotation in rad
    v(:,t) = vec1 + tan(angle(t))*vec2;   % length of vec1 is 1
    v(:,t) = v(:,t)/norm(v(:,t));                % v is an unit vector 5 degrees from vec1
    dip_bst.Dipole(t) = dip_bst.Dipole(1);
    dip_bst.Dipole(t).Amplitude = v(:,t);
    dip_bst.Dipole(t).Time = timeVec(t);
end

dipReady = dip_bst;

u=0;
for bu = 1:100
    for st = 1: 50
        dipReady.Dipole(:,st+u) = dip_bst.Dipole(bu); % find the index of a dipole (findHM) in the IMAGEGRIDAMP.
    end
    u = u +50;
end

dipReady.Time = linspace(0,10,length(dip_bst.Dipole)*50);
for h = 1:length(dipReady.Time)
    dipReady.Dipole(h).Time = dipReady.Time(h); %adjust the time vector to 5000 sample points
end

IMAGEGRIDAMP = zeros(length(Source.ImageGridAmp),length(dipReady.Time));

for st = 1: length(dipReady.Time)
    IMAGEGRIDAMP(findHM(1)*3-2:findHM(1)*3,st) = dipReady.Dipole(st).Amplitude; % find the index of a dipole (findHM) in the IMAGEGRIDAMP.
end %here, we write the dipole signals into the source structure

%% From the dipoles, create simulations of the resulting topographies

Source.ImageGridAmp = IMAGEGRIDAMP;
Source.Time = dipReady.Time;

listRaw = dir([MAINPATH,'brainstorm_db\',PROTOCOLNAME,'\data\',SubjectName,'\raw_plots*\data_2*']);
Source.DataFile = [SubjectName,'/raw_plots/',listRaw(end).name]; % attach the source to the raw signal
dipReady.DataFile = [SubjectName,'/raw_plots/',listRaw(end).name]; % attach the dipole file to the raw signal
db_add(ProtocolInfo.iStudy, Source, 1); %load the source into bst
db_add(ProtocolInfo.iStudy, dipReady, 1); %the dipole file, too

%listResults = dir([MAINPATH,'brainstorm_db\',PROTOCOLNAME,'\data\',SubjectName,'\raw_plots*\results_20*']);
%Sim3DElec = bst_simulation([listResults(end).folder,'\',listResults(end).name]); %simulate from the new source for the 3DElectrodes (only cEEGrid)

%% Plot the results

listDipNew =  dir([MAINPATH,'brainstorm_db\',PROTOCOLNAME,'\data\',SubjectName,'\raw_plots*\dipoles_2*']);
listSim = dir([MAINPATH,'brainstorm_db\',PROTOCOLNAME,'\data\',SubjectName,'\raw_plots*\data_2*']);

%dipole orientations
%hFig = view_dipoles([listDipNew(end).folder,'\',listDipNew(end).name], 'cortex')
%set(gcf, 'Position', get(0, 'ScreenSize'),'Visible', 'on');

%cEEGrid elec-topo from the left side
%leftFig = view_topography([listSim(end).folder,'\',listSim(end).name], 'EEG', '3DElectrodes'); 
%figure_3d('SetStandardView', gcf, {'left'});

%Topoplots
%figTopo3D = view_topography([listSim(end).folder,'\',listSim(end).name], 'EEG', '3DSensorCap'); 
figTopo3D = view_topography([listSim(end).folder,'\',listSim(end).name], 'MEG', '3DSensorCap'); 
%figTopo2D = view_topography([listSim(end).folder,'\',listSim(end).name], 'EEG', '2DSensorCap'); 
figTopo2D = view_topography([listSim(end).folder,'\',listSim(end).name], 'MEG', '2DSensorCap')


%% clean the workspace in bst

% listDeleteAll = [MAINPATH,'brainstorm_db\',PROTOCOLNAME,'\data\',SubjectName,'\raw_plots'];
% isDeleted = file_delete( listDeleteAll, 0, 1);
% db_reload_studies(ProtocolInfo.iStudy) %reload to get rid of the deleted files in the GUI
