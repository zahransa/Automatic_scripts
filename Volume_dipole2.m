%% Here, we use the source- and dipole structure to create single dipoles with different orientations

%%%%%%%%%%%%SETTINGS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SubjectName = 'Subject03'; %specify the name of the subject that will be used in the brainstorm GUI
MAINPATH = 'C:\'; %path to brainstorm3
PROTOCOLNAME = 'Protocol01'; %name of your current protocol in brainstorm
%%%%%%%%%%%%SETTINGS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load source matrix, channel matrix, head model matrix and dipole matrix that were created in the previous script
listSource    = dir([MAINPATH,'brainstorm_db\',PROTOCOLNAME,'\data\',SubjectName,'\raw\results_MN*']); %lists of the variables created in the first script
listChannel   = dir([MAINPATH,'brainstorm_db\',PROTOCOLNAME,'\data\',SubjectName,'\raw\channel*']);
listHeadmodel = dir([MAINPATH,'brainstorm_db\',PROTOCOLNAME,'\data\',SubjectName,'\raw\headmodel*']);
listDip       = dir([MAINPATH,'brainstorm_db\',PROTOCOLNAME,'\data\',SubjectName,'\raw\dipoles*']);%load the dipole file

Source = load([listSource(1).folder,'\',listSource(1).name]); %load the source template
CM     = load([listChannel(end).folder,'\',listChannel(end).name]); %...and the channelfile
HM     = load([listHeadmodel(end).folder,'\',listHeadmodel(end).name]); %...and the head model
dip    = load([listDip(end).folder,'\',listDip(end).name]); %...amd the dipole matrix

%% Choose a dipole

NormOrient = [-0.0434816861376882,0.985049224918913,-0.103375097755849;0.993581456236887,0.0542738434228215,0.0992483740033164;-0.104270505035261,0.0992483740033164,0.989584469379643]; 
%These are 3 different orientations. This is the variable that needs to be
%altered when simulating different orientations! In this example, this
%creates 3 dipoles, with the same location and orthogonal orientations.

%% put the new source locations and orientations into the dipole file and load it to brainstorm

dip.Dipole(1).Loc = HM.GridLoc(14,:)'; %3 dipoles with the same location and differ only in orientation
dip.Dipole(2).Loc = HM.GridLoc(14,:)'; % 143 is a random location from the GridLoc, it can be altered to be any position
dip.Dipole(3).Loc = HM.GridLoc(14,:)'; % HM.GridLoc is the variable that determines the position of the dipole!

dip.Dipole(1).Amplitude = NormOrient(1,:)'; %here, the orientations for the loc are defined
dip.Dipole(2).Amplitude = NormOrient(2,:)';
dip.Dipole(3).Amplitude = NormOrient(3,:)';

dip.Dipole(4:end) = []; %delete the remaining values from the predefined dipole file. Only the once specified above remain
dip.Time(4:end) = []; %same with the time points. You can specify as many dipoles and orientations as you want for your own setup

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
    [~,findHM(w,:)] = ismember(dip.Dipole(u).Loc',HM.GridLoc,'rows'); %this loop finds for each wanted dipole location the equivalent on the HeadModel grid
    w=w+1; %only really useful when you have other dipole locations specified; in this example, it is always the same
end
%% Preapre an empty SourceGrid
%Now that we have the dipoles prepared, we enter them into the empty source
%grid and forward model them. This is to get the topographies of the signal

Source.ImageGridAmp = zeros(size(Source.ImageGridAmp)); %set source matrix to zero
Source.DataFile = []; % needed to avoid brainstorm confusion, otherwise it might try to reload things

Source.Time = 1;
Source.ImageGridAmp = Source.ImageGridAmp(:,1); %computationally less expensive

sFiles = [];
sFiles = {...
    [listSource(1).folder,'\',listSource(1).name]};

% Start a new report, get subject infos
bst_report('Start', sFiles);
sProcess = bst_process('GetInputStruct', sFiles);
[sSubject,iSubject] = bst_get('Subject', sProcess.SubjectFile); %gets the subject ID
ProtocolInfo = bst_get('ProtocolInfo'); %needed for number of the current protocol. ALWAYS click the nedded protocol AND subject in the GUI
iChannels = 1:length(CM.Channel); %number of channels according to your channel file

for d = 1: length(dip.Dipole) %loop through the dipoles with their 3 orientations
    SourceSimple = Source; %reset the source variable in every iteration
    dip_bst = dip; %reset the dipole variable in every iteration
    dip_bst.Dipole = dip_bst.Dipole(d); %take only the first dipole with the first orientation
    dip_bst.Time = 0;
    
    ProtocolInfo = bst_get('ProtocolInfo');
    dip_bst.Comment= ['dip_',num2str(d)];
    db_add(ProtocolInfo.iStudy, dip_bst, []); %add the dipole to brainstorm to have the plot of the dipole orientation
    
    IMAGEGRIDAMP = SourceSimple.ImageGridAmp; %This is the source structure
    IMAGEGRIDAMP(findHM(d)*3-2:findHM(d)*3,1) = dip_bst.Dipole.Amplitude; % find the position of a dipole (findHM) in the IMAGEGRIDAMP.
    %since per source there are 3 orientation, multiply the loc, which is
    %only one point (point 1 in this example), by 3. Replace these values with those ones of the 3
    %orientations from dip_bst.Dipole(r).Amplitude. The result is one
    %dipole source with one orientation
    
    SourceSimple.ImageGridAmp = IMAGEGRIDAMP;
    
    ProtocolInfo = bst_get('ProtocolInfo'); %we now have the simulated dipole source. We reload it into bst
    SourceSimple.Comment = ['dip_',num2str(d)];
    db_add(ProtocolInfo.iStudy, SourceSimple, []);
    db_reload_studies(sProcess.iStudy) %reload to get rid of the deleted files in the GUI
    
    %% simulate the recordings based on the dipoles
    listResults = dir([MAINPATH,'brainstorm_db\',PROTOCOLNAME,'\data\',SubjectName,'\raw\results_20*']);
    newDataFile = bst_simulation([listResults(end).folder,'\',listResults(end).name]);
    %we now forward model the source signal. This generates the topographies on the scalp surface. This step is the actual forward modeling of the source!
    
    listSimulation = dir([MAINPATH,'brainstorm_db\',PROTOCOLNAME,'\data\',SubjectName,'\raw\*simulation*.mat']);
    
    %rightFigTopo = view_topography([listSimulation(end).folder,'\',listSimulation(end).name], 'EEG', '2DSensorCap');
    rightFigTopo = view_topography([listSimulation(end).folder,'\',listSimulation(end).name], 'MEG', '2DSensorCap');
    set(gcf, 'Position', get(0, 'ScreenSize'),'Visible', 'on');
    %this displays the simulated signal as a topoplot. You can visualize
    %all recordings in the GUI as well, for example to display the dipoles
    %(right click the dipole symbol, "Display on MRI (3D)")
   
end
