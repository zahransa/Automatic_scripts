%% This script simulates a source and gives back a sensitivity map of the signals that can be recorded with the electrodes

%%%%%%%%%%%%SETTINGS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SubjectName = 'Subject03'; %the subject must be the one with the SURFACE HEADMODEL, not the volume one
MAINPATH = 'C:\'; %path to brainstorm3
PROTOCOLNAME = 'Protocol01'; %name of your current protocol in brainstorm
PATHOUT = 'C:\Users\Said\Desktop\alter_dipole\results\'; %store the results here
%%%%%%%%%%%%SETTINGS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% load source matrix, channel matrix, head model matrix and dipole matrix that were created in the previous script
ProtocolInfo = bst_get('ProtocolInfo'); %needed for number of the current protocol. ALWAYS click the nedded protocol AND subject in the GUI

listSource = dir([MAINPATH,'brainstorm_db\',PROTOCOLNAME,'\data\',SubjectName,'\*aw*\results_MN*']); %lists of the variables created in PrepareSource_cEEGrid.mat
listChannel = dir([MAINPATH,'brainstorm_db\',PROTOCOLNAME,'\data\',SubjectName,'\*aw*\channel*']); %raw s.t. starts with R or r
listHeadmodel = dir([MAINPATH,'brainstorm_db\',PROTOCOLNAME,'\data\',SubjectName,'\*aw*\headmodel*']);
listAnat = dir([MAINPATH,'brainstorm_db\',PROTOCOLNAME,'\anat\@default_subject\*cortex_pial_low*.mat']);

Source = load([listSource(1).folder,'\',listSource(1).name]); %load the source template that was created in a_PrepareSources.m...
CM = load([listChannel(end).folder,'\',listChannel(end).name]); %...and the channelfile
HM = load([listHeadmodel(1).folder,'\',listHeadmodel(1).name]); %...and the headmodel
HM.Gain = bst_gain_orient(HM.Gain, HM.GridOrient); %modify the gain-matrix, get constrained orientations (1 coordinate instead of xyz)
ANAT = load([listAnat(1).folder,'\',listAnat(1).name]);

Source.ImageGridAmp= zeros(size(Source.ImageGridAmp)); %set source matrix to zero
Source.DataFile=[]; % needed to avoid brainstorm confusion, otherwise it might try to reload things

Source.Time = 1;
Source.ImageGridAmp = Source.ImageGridAmp(:,1); %computationally less expensive, make sure to only have 1 col. More is unnecessary, the signal is always 1
SourceDist = Source; % source variable to work in

%calculate the sensitivity for a source (pointwise or for a cluster (ROI
%from atlas)) for the cEEGrids

nameResult = 'Cap';
nameVar = 'Full_';
iChannels = 1:length(CM.Channel); %number of channels according to your channel file



for  k = 1:size(Source.ImageGridAmp,1) %loop through every point on the brain mesh
    
    IMAGEGRIDAMP = Source.ImageGridAmp; %This is the source/scout
    IMAGEGRIDAMP(k,1) = 1; %put values into the source grid
    
    F = zeros(length(iChannels), 1); % simulation matrix of the signals on the electrode level
    F(iChannels,:) = HM.Gain(iChannels,:) * IMAGEGRIDAMP; %for every channel, calc the forward model of the simulated source with amp=1 (leadfield*activation)
    
    collect=[];
    index=1;
    
    for ch=1:length(iChannels)-1
        for b=ch+1:length(iChannels)
            collect(index)=abs((F(iChannels(ch),1))- (F(iChannels(b),1))); %calc every combination of electrodes. The result is a matrix of all bipolar electrode combinations
            index=index+1;
        end
    end
        
    SourceDist.ImageGridAmp(k,1) = max(collect); %get the maximum electrode, aka the one with the highest signal
    %these values indicate the highest signal that can be recorded with
    %the current electrode setup. The values are written into
    %SourceDist. When now displaying the source structure in bst, we
    %see for every vertex the highest signal that can be recorded from
    %that position.
    
    if mod(k,1000)==0
        disp('Almost there')
    end
    
end
[a,b] = maxk(SourceDist.ImageGridAmp(:,1),10); %find the 100 highest values
c = find(a > 200); %find unrealisticly high vlaues
SourceDist.ImageGridAmp(b(c),1)=0; % delete them

SourceDist.Comment= [nameResult,nameVar,'GridFull']; %set the name referring to the chosen grid (left/right/double)
SourceDist.ImageGridAmp = abs(SourceDist.ImageGridAmp); % use absolute values

ProtocolInfo = bst_get('ProtocolInfo'); %needed for number of the current protocol
db_add(ProtocolInfo.iStudy, SourceDist, []); %this adds the new source projection to the current protocol
