function [RTOTime, LTOTime, RHSTime, LHSTime, commSendTime, commSendFrame] = Speed_audioFeedback(velL,velR,FzThreshold,profilename,mode,signList,paramComputeFunc,paramCalibFunc)
%This function takes two vectors of speeds (one for each treadmill belt)
%and succesively updates the belt speed upon ipsilateral Toe-Off
%The function only updates the belts alternatively, i.e., a single belt
%speed cannot be updated twice without the other being updated
%The first value for velL and velR is the initial desired speed, and new
%speeds will be sent for the following N-1 steps, where N is the length of
%velL

%%
%close all
nirs = 1;
oxysoft_present = 1;
alphabetLetter = 'B'; %TODO: to be randomized for each participant 1st and last visit swap; maybe read from randomized file generated order
%         eventorder_all = [0,1,2; 2, 1, 0; 1,2,0; 1,0,2; 0,2,1; 2,0,1];%TODO: to be randomized for each participant
event_list = [0,1,2];
iterations = 6;

% Oxysoft = nan; %FIXME: for testing without nirs
Oxysoft = actxserver('OxySoft.OxyApplication'); %FIXME: uncomment

global feedbackFlag
if feedbackFlag==1  %&& size(get(0,'MonitorPositions'),1)>1
    ff=figure('Units','Normalized','Position',[1 0 1 1]);
    pp=gca;
    axis([0 length(velL)  0 2]);
    ccc1=animatedline('Parent',pp,'Marker','o','LineStyle','none','MarkerFaceColor',[0 0 1],'MarkerEdgeColor','none');
    ccc2=animatedline('Parent',pp,'Marker','o','LineStyle','none','MarkerFaceColor',[1 0 0],'MarkerEdgeColor','none');
    
end

% NET.addAssembly('System.Speech');
% obj = System.Speech.Synthesis.SpeechSynthesizer;
    
paramLHS=0;
paramRHS=0;
lastParamLHS=0;
lastParamRHS=0;
lastRstepCount=0;
lastLstepCount=0;
LANK=zeros(1,6);
RANK=zeros(1,6);
RHIP=zeros(1,6);
LHIP=zeros(1,6);
LANKvar=1e5*ones(6);
RANKvar=1e5*ones(6);
RHIPvar=1e5*ones(6);
LHIPvar=1e5*ones(6);
lastLANK=nan(1,6);
lastRANK=nan(1,6);
lastRHIP=nan(1,6);
lastLHIP=nan(1,6);
toneplayed=false;

if nargin<7 || isempty(paramComputeFunc)
    paramComputeFunc=@(w,x,y,z) (-w(2) + .5*(y(2)+z(2)))/1000; %Default
end
[d,~,~]=fileparts(which(mfilename));
if nargin<8 || isempty(paramCalibFunc)
    load([d '\..\calibrations\lastCalibration.mat'])
else
    save([d '\..\calibrations\lastCalibration.mat'],'paramCalibFunc','paramComputeFunc')
end




if nargin<5 || isempty(mode)
    error('Need to specify control mode')
end
if nargin<6
    %disp('no sign list')
    signList=[];
end
signCounter=1;

global PAUSE%pause button value
global STOP
global memory
global enableMemory
global firstPress
STOP = false;
fo=4000;
signS=0; %Initialize
endTone=3*sin(2*pi*[1:2048]*fo*.025/4096);  %.5 sec tone at 100Hz, even with noise cancelling this is heard
tone=sin(2*pi*[1:2048]*1.25*fo/4096); %.5 sec tone at 5kHz, even with noise cancelling this is heard

ghandle = guidata(AdaptationGUI);%get handle to the GUI so displayed data can be updated

%Initialize some graphical objects:
ppp1=animatedline('Parent',ghandle.profileaxes,'LineStyle','none','Marker','o','MarkerFaceColor',[0.68 .92 1],'MarkerEdgeColor','none');
ppp2=animatedline('Parent',ghandle.profileaxes,'LineStyle','none','Marker','o','MarkerFaceColor',[1 .6 .78],'MarkerEdgeColor','none');
ppp3=animatedline('Parent',ghandle.profileaxes,'LineStyle','none','Marker','.','Color',[0 0 1]);
ppp4=animatedline('Parent',ghandle.profileaxes,'LineStyle','none','Marker','.','Color',[1 0 0]);
ppv1=animatedline('Parent',ghandle.profileaxes,'LineStyle','none','Marker','.','Color',[0 .5 1]);
ppv2=animatedline('Parent',ghandle.profileaxes,'LineStyle','none','Marker','.','Color',[1 .5 0]);

yl=get(ghandle.profileaxes,'YLim');
%Add text with signlist:
%tt=text(1,yl(2)-.05*diff(yl),['Sign List for null trials= ' num2str(signList)],'Parent',ghandle.profileaxes);
%Add patches to show self-controlled strides:
aux1=diff(isnan(velL));
iL=find(aux1==1);
iL2=find(aux1==-1);
counter=1;
for i=1:length(iL)
    pp=patch([iL(i) iL2(i) iL2(i) iL(i)],[yl(1) yl(1) yl(2) yl(2)],[0 .3 .6],'FaceAlpha',.2,'EdgeColor','none');
    uistack(pp,'bottom')
    if (velL(iL(i))-velR(iL(i)))==0
        try
            auxT=signList(counter);
        catch
            if mode==0
                warning('Provided sign list seems to be shorter than # of null trials, adding default sign.')
            end
            signList(counter)=1;
            auxT=signList(counter);
        end
        counter=counter+1;
    else
        auxT=sign(velR(iL(i))-velL(iL(i)));
    end
    text(iL(i),yl(1)+.05*diff(yl),[num2str(auxT)],'Color',[0 0 1])
end
aux1=diff(isnan(velR));
iL=find(aux1==1);
iL2=find(aux1==-1);
counter=1;
for i=1:length(iL)
    pp=patch([iL(i) iL2(i) iL2(i) iL(i)],[yl(1) yl(1) yl(2) yl(2)],[.6 .3 0],'FaceAlpha',.2,'EdgeColor','none');
    uistack(pp,'bottom')
    if (velL(iL(i))-velR(iL(i)))==0
        try
            auxT=-signList(counter);
        catch
            if mode==0
                warning('Provided sign list seems to be shorter than # of null trials, adding default sign.')
            end
            signList(counter)=1;
            auxT=-signList(counter);
        end
        counter=counter+1;
    else
        auxT=sign(velL(iL(i))-velR(iL(i)));
    end
    text(iL(i),yl(1)+.1*diff(yl),[num2str(auxT)],'Color',[1 0 0])
end

% %initialize a data structure that saves information about the trial,
% %Shuqi TODO
datlog = struct();
datlog.buildtime = now;%timestamp
temp = datestr(datlog.buildtime,30); %Using ISO 8601 for date/time string format
a = regexp(temp,'-');
temp(a) = '_';
b = regexp(temp,':');
temp(b) = '_';
c = regexp(temp,' ');
temp(c) = '_';
[d,~,~]=fileparts(which(mfilename));
[~,n,~]=fileparts(profilename);
savename = [[d '\..\datlogs\'] temp '_' n];
set(ghandle.sessionnametxt,'String',[temp '_' n]);
datlog.session_name = savename;
datlog.errormsgs = {};
datlog.messages = {};
datlog.framenumbers.header = {'frame #','U Time','Relative Time'};
datlog.framenumbers.data = zeros(300*length(velR)+7200,2);
datlog.stepdata.header = {'Step#','U Time','frame #','Relative Time'};
datlog.stepdata.RHSdata = zeros(length(velR)+50,3);%empty cell theoretically big enough to house all the steps taken
datlog.stepdata.RTOdata = zeros(length(velR)+50,3);
datlog.stepdata.LHSdata = zeros(length(velL)+50,3);
datlog.stepdata.LTOdata = zeros(length(velL)+50,3);
datlog.stepdata.paramLHS = zeros(length(velR)+50,1);
datlog.stepdata.paramRHS = zeros(length(velL)+50,1);
% datlog.audioCues.start = nan(length(velL)+50,1); %Shuqi TODO
% datlog.audioCues.stop = nan(length(velL)+50,1);
datlog.audioCues.start = [];
datlog.audioCues.stop = [];
datlog.audioCues.audio_instruction_message = {};
histL=nan(1,50);
histR=nan(1,50);
histCount=1;
filtLHS=0;
filtRHS=0;
datlog.inclineang = [];
datlog.speedprofile.velL = velL;
datlog.speedprofile.velR = velR;
datlog.speedprofile.signListForNullTrials = signList;
datlog.TreadmillCommands.header = {'RBS','LBS','angle','U Time','Relative Time'};
datlog.TreadmillCommands.read = nan(300*length(velR)+7200,4);
datlog.TreadmillCommands.sent = nan(300*length(velR)+7200,4);%nan(length(velR)+50,4);
datlog.Markers.header = {'x','y','z'};
datlog.Markers.LHIP = nan(300*length(velR)+7200,3);
datlog.Markers.RHIP = nan(300*length(velR)+7200,3);
datlog.Markers.LANK = nan(300*length(velR)+7200,3);
datlog.Markers.RANK = nan(300*length(velR)+7200,3);
datlog.Markers.LHIPfilt = nan(300*length(velR)+7200,6);
datlog.Markers.RHIPfilt = nan(300*length(velR)+7200,6);
datlog.Markers.LANKfilt = nan(300*length(velR)+7200,6);
datlog.Markers.RANKfilt = nan(300*length(velR)+7200,6);
datlog.Markers.LHIPvar = nan(300*length(velR)+7200,36);
datlog.Markers.RHIPvar = nan(300*length(velR)+7200,36);
datlog.Markers.LANKvar = nan(300*length(velR)+7200,36);
datlog.Markers.RANKvar = nan(300*length(velR)+7200,36);
datlog.stepdata.paramComputeFunc=func2str(paramComputeFunc);
if mode==4
    datlog.stepdata.paramCalibFunc=func2str(paramCalibFunc);
end

%do initial save
try
    save(savename,'datlog');
catch ME
    disp(ME);
end
%Default threshold
if nargin<3
    FzThreshold=40; %Newtons (40 is minimum for noise not to be an issue)
elseif FzThreshold<40
    %     warning = ['Warning: Fz threshold too low to be robust to noise, using 30N instead'];
    datlog.messages{end+1} = 'Warning: Fz threshold too low to be robust to noise, using 40N instead';
    disp('Warning: Fz threshold too low to be robust to noise, using 40N instead');
    FzThreshold=40;
end


%Check that velL and velR are of equal length
N=length(velL)+1;
if length(velL)~=length(velR)
    disp('WARNING, velocity vectors of different length!');
    datlog.messages{end+1} = 'Velocity vectors of different length selected';
end
%Adding a fake step at the end to avoid out-of-bounds indexing, but those
%speeds should never get used in reality:
velL(end+1)=0;
velR(end+1)=0;

%Initialize nexus & treadmill communications
try
    % [MyClient] = openNexusIface();
    Client.LoadViconDataStreamSDK();
    MyClient = Client();
    Hostname = 'localhost:801';
    out = MyClient.Connect(Hostname);
    out = MyClient.EnableMarkerData();
    out = MyClient.EnableDeviceData();
    %MyClient.SetStreamMode(StreamMode.ServerPush);
    MyClient.SetStreamMode(StreamMode.ClientPullPreFetch);
    mn={'LHIP','RHIP','LANK','RANK'};
    altMn={'LGT','RGT','LANK','RANK'};
catch ME
    disp('Error in creating Nexus Client Object/communications see datlog for details');
    datlog.errormsgs{end+1} = 'Error in creating Nexus Client Object/communications';
    datlog.errormsgs{end+1} = ME;%store specific error
    disp(ME);
end
try
    t = openTreadmillComm();
catch ME
    disp('Error in creating TCP connection to Treadmill, see datlog for details...');
    datlog.errormsgs{end+1} = 'Error in creating TCP connection to Treadmill';
    datlog.errormsgs{end+1} = ME;
    disp(ME);
    %     log=['Error ocurred when opening communications with Nexus & Treadmill'];
    %     listbox{end+1}=log;
    %     disp(log);
end

try %So that if something fails, communications are closed properly
    % [FrameNo,TimeStamp,SubjectCount,LabeledMarkerCount,UnlabeledMarkerCount,DeviceCount,DeviceOutputCount] = NexusGetFrame(MyClient);
    MyClient.GetFrame();
    % listbox{end+1} = ['Nexus and Bertec Interfaces initialized: ' num2str(clock)];
    datlog.messages{end+1} = ['Nexus and Bertec Interfaces initialized: ' num2str(now)];
    % set(ghandle.listbox1,'String',listbox);
    
    %Initiate variables
    new_stanceL=false;
    new_stanceR=false;
    phase=0; %0= Double Support, 1 = single L support, 2= single R support
    LstepCount=1;
    RstepCount=1;
    % RTOTime(N)=TimeStamp;
    % LTOTime(N)=TimeStamp;
    % RHSTime(N)=TimeStamp;
    % LHSTime(N)=TimeStamp;
    RTOTime(N) = now;
    LTOTime(N) = now;
    RHSTime(N) = now;
    LHSTime(N) = now;
    commSendTime=zeros(2*N-1,6);
    commSendFrame=zeros(2*N-1,1);
    % stepFlag=0;
    
    fastest = 5.6; %7 meters/ 1.25 m/s
    slowest = 7; % 7 meters/ 1 m/s
    y_max = 4500;
    y_min = -2500;
 
    % TODO Shuqi NIRS
    inout1 = 0; %initialize both variables to 0
    inout2 = 0;
    
    LHS_time = 0;
    RHS_time = 0;
    LHS_pos = 0;
    RHS_pos = 0;
    HS_frame = 0;
    pad = 50;
    
    if ~nirs
        audioids = {'fast','good','slow'};
        instructions = containers.Map();
        for i = 1 : length(audioids)
            [audio_data,audio_fs]=audioread(strcat(audioids{i},'.mp3'));
            instructions(audioids{i}) = audioplayer(audio_data,audio_fs);
        end
    else
        %Connect to Oxysoft - Shuqi TODO
        restDuration = 20; %FIXME: time when rest, usually 20,
        disp('Initial Setup')
        eventorder_all = repmat(event_list, iterations, 1);
        trialIndex = 1;
        eventorder = eventorder_all(trialIndex,:); %0 = stand and alphabet, 1 = walk and alphabet, 2 = walk, 
        %set up audio players
        if alphabetLetter == 'A'
            audioFileNames = {'relax','walk','rest','standAndA','stopAndRest','turnAndStop','walkAndA'};
        else
            audioFileNames = {'relax','walk','rest','standAndB','stopAndRest','turnAndStop','walkAndB'};
        end
        audioids = {'relax','walk','rest','standAlphabet','stopAndRest','turnAndStop','walkAlphabet'};
        instructions = containers.Map();
        for i = 1 : length(audioFileNames)
            [audio_data,audio_fs]=audioread(strcat(audioFileNames{i},'.mp3'));
            instructions(audioids{i}) = audioplayer(audio_data,audio_fs);
        end
        
        % Write event I with description 'Instructions' to Oxysoft
        datlog = nirsEvent('', 'I', 'Connected', instructions, datlog, Oxysoft, oxysoft_present);
        
        %start with a rest block
        disp('Rest');
        datlog = nirsEvent('rest', 'R', 'Rest', instructions, datlog, Oxysoft, oxysoft_present);
        pause(restDuration);

        currentIndex = 1;
        if eventorder(currentIndex)== 0 %stand and alphabet
            datlog = nirsEvent('standAlphabet', 'A', ['Stand and Alphabet with ' alphabetLetter], instructions, datlog, Oxysoft, oxysoft_present);
            pause(restDuration)

            %complete follow a rest block
            datlog = nirsEvent('stopAndRest','R','Rest', instructions, datlog, Oxysoft, oxysoft_present);
            pause(restDuration);
            currentIndex = currentIndex + 1;
        end

        if eventorder(currentIndex)==1 %first event is walk and alphabet
            datlog = nirsEvent('walkAlphabet','B',['Walk and Alphabet' alphabetLetter], instructions, datlog, Oxysoft, oxysoft_present);
            currentIndex  = currentIndex + 1;
        elseif eventorder(currentIndex) == 2 %firstw event is walk
            datlog = nirsEvent('walk','W','Walk', instructions, datlog, Oxysoft, oxysoft_present);
            currentIndex = currentIndex + 1;
        end
    end
    
    %TODO: Shuqi
    if nirs
        y_max = 4250;%0 treadmill, max = 7 tile = 14 ft ~= 4.2672meter
        y_min = -2250;%0 treadmill, min = 2.3 tile = 4.6 ft ~= 1.4021 meter
%             y_max = 4500;
%     y_min = -2500;
    end

    
    %Send first speed command & store
    acc=3000;
%     [payload] = getPayload(velR(1),velL(1),acc,acc,cur_incl);
%     memoryR=velR(1);
%     memoryL=velL(1);
%     %sendTreadmillPacket(payload,t);
%     lastSent=now;
%     datlog.TreadmillCommands.firstSent = [velR(RstepCount),velL(LstepCount),acc,acc,cur_incl,lastSent];%record the command
    commSendTime(1,:)=clock;
%     datlog.TreadmillCommands.sent(1,:) = [velR(RstepCount),velL(LstepCount),cur_incl,now];%record the command
%     datlog.messages{end+1} = ['First speed command sent' num2str(now)];
%     datlog.messages{end+1} = ['Lspeed = ' num2str(velL(LstepCount)) ', Rspeed = ' num2str(velR(RstepCount))];
%     
    %% Main loop
    
    old_velR = libpointer('doublePtr',velR(1));
    old_velL = libpointer('doublePtr',velL(1));
    frameind = libpointer('doublePtr',1);
    framenum = libpointer('doublePtr',0);
    
    %Get first frame
    MyClient.GetFrame();
    framenum.Value = MyClient.GetFrameNumber().FrameNumber;
    datlog.framenumbers.data(frameind.Value,:) = [framenum.Value now];
    
    while ~STOP %only runs if stop button is not pressed
        pause(.001) %This ms pause is required to allow Matlab to process callbacks from the GUI.
        %It actually takes ~2ms, as Matlab seems to be incapable of pausing for less than that (overhead of the pause() function, I presume)
        while PAUSE %only runs if pause button is pressed
            pause(.2);
            datlog.messages{end+1} = ['Loop paused at ' num2str(now)];
            disp(['Paused at ' num2str(clock)]);
%             %bring treadmill to a stop and keep it there!...
%             [payload] = getPayload(0,0,500,500,cur_incl);
%             %cur_incl
%             sendTreadmillPacket(payload,t);
%             %do a quick save
            try
                save(savename,'datlog');
            catch ME
                disp(ME);
            end
            old_velR.Value = 1;%change the old values so that the treadmill knows to resume when the pause button is resumed
            old_velL.Value = 1;
        end
        %newSpeed
        %drawnow;
        %     lastFrameTime=curTime;
        %     curTime=clock;
        %     elapsedFrameTime=etime(curTime,lastFrameTime);
        old_stanceL=new_stanceL;
        old_stanceR=new_stanceR;
        
        
        
        %Read frame, update necessary structures4
        MyClient.GetFrame();
        framenum.Value = MyClient.GetFrameNumber().FrameNumber;
    
        if framenum.Value~= datlog.framenumbers.data(frameind.Value,1) %Frame num value changed, reading data
            frameind.Value = frameind.Value+1; %Frame counter
            datlog.framenumbers.data(frameind.Value,:) = [framenum.Value now];
            
            %Read treadmill, if enough time has elapsed since last read
%             aux=(datevec(now)-datevec(lastRead));
%             if aux(6)>.1 || any(aux(1:5)>0)  %Only read if enough time has elapsed
% %                 [RBS, LBS,read_theta] = readTreadmillPacket(t);%also read what the treadmill is doing
%                 lastRead=now;
%                 datlog.TreadmillCommands.read(frameind.Value,:) = [RBS,LBS,read_theta,lastRead];%record the read
%                 %set(ghandle.RBeltSpeed_textbox,'String',num2str(RBS/1000));
%                 %set(ghandle.LBeltSpeed_textbox,'String',num2str(LBS/1000));
%             end
            
            %Read markers:
            sn=MyClient.GetSubjectName(1).SubjectName;
            for l=1
                md=MyClient.GetMarkerGlobalTranslation(sn,mn{l});
                
                md_LHIP=MyClient.GetMarkerGlobalTranslation(sn,altMn{1});
                md_RHIP=MyClient.GetMarkerGlobalTranslation(sn,altMn{2});
                md_LANK=MyClient.GetMarkerGlobalTranslation(sn,mn{3});
                md_RANK=MyClient.GetMarkerGlobalTranslation(sn,mn{4});
                
                if md.Result.Value==2 %%Success getting marker
                    aux=double(md.Translation);
                    
                    
                    LHIP_pos = double(md_LHIP.Translation);
                    RHIP_pos = double(md_RHIP.Translation);
                    LANK_pos = double(md_LANK.Translation);
                    RANK_pos = double(md_RANK.Translation);
                else
                    md=MyClient.GetMarkerGlobalTranslation(sn,altMn{l});
                    if md.Result.Value==2
                        aux=double(md.Translation);
                        LHIP_pos = double(md_LHIP.Translation);
                        RHIP_pos = double(md_RHIP.Translation);
                        LANK_pos = double(md_LANK.Translation);
                        RANK_pos = double(md_RANK.Translation);
                    else
                        aux=[nan nan nan]';
                        LHIP_pos = [nan nan nan]';
                        RHIP_pos = [nan nan nan]';
                        LANK_pos = [nan nan nan]';
                        RANK_pos = [nan nan nan]';
                    end
                    if all(aux==0) || all(LHIP_pos==0) || all(RHIP_pos==0) || all(LANK_pos==0) || all(RANK_pos==0)
                        %Nexus returns 'Success' values even when marker is missing, and assigns it to origin!
                        %warning(['Nexus is returning origin for marker position of ' mn{l}]) %This shouldn't happen
                        aux=[nan nan nan]';
                        LHIP_pos = [nan nan nan]';
                        RHIP_pos = [nan nan nan]';
                        LANK_pos = [nan nan nan]';
                        RANK_pos = [nan nan nan]';
                    end
                end
                %aux
                %Here we should implement some Kalman-like filter:
                %eval(['v' mn(l) '=.8*v' mn(l) '+ (aux-' mn(l) ');']);
                eval(['lastEstimate=' mn{l} ''';']); %Saving previous estimate
                eval(['last' mn{l} '=lastEstimate'';']); %Saving previous estimate as row vec
                eval(['lastEstimateCovar=' mn{l} 'var;']); %Saving previous estimate's variance
                eval(['datlog.Markers.' mn{l} '(frameind.Value,:) = aux;']);%record the current read
                
                elapsedFramesSinceLastRead=framenum.Value-datlog.framenumbers.data(frameind.Value-1,1);
                %What follows assumes a 100Hz sampling rate
                n=elapsedFramesSinceLastRead;
                if n<0
                    error('inconsistent frame numbering')
                end
                
                
                %% added by Yashar to count steps OG and for verbal feedback action

                set(ghandle.text25,'String',num2str(framenum.Value/100));
                body_y_pos(frameind.Value) = (LHIP_pos(2)+RHIP_pos(2))/2;
                body_y_pos_diff(frameind.Value) = body_y_pos(frameind.Value) - body_y_pos(frameind.Value-1);
                Ank_diff(frameind.Value) = LANK_pos(2) - RANK_pos(2);
                Ank_diffdiff(frameind.Value) = Ank_diff(frameind.Value)-Ank_diff(frameind.Value-1);
                
                if body_y_pos(frameind.Value) > y_min && body_y_pos(frameind.Value) < y_max && framenum.Value-HS_frame>pad
                    
                    if body_y_pos_diff(frameind.Value) >= 0 && Ank_diffdiff(frameind.Value) >= 0 && sign(Ank_diffdiff(frameind.Value)) ~= sign(Ank_diffdiff(frameind.Value-1)) 
                        Left_HS(frameind.Value) = 1;
                        sumiL = sum(Left_HS);
                        HS_frame = framenum.Value;
                        LHS_time = [LHS_time framenum.Value/100];
                        LHS_pos = [LHS_pos body_y_pos(frameind.Value)/1000];
                        OG_speed_left = (LHS_pos(end)-LHS_pos(end-1))/(LHS_time(end)-LHS_time(end-1));
                        
                        set(ghandle.LBeltSpeed_textbox,'String',num2str(OG_speed_left));
                        set(ghandle.Left_step_textbox,'String',num2str(sumiL));
                        
                        addpoints(ppp3,sumiL,OG_speed_left)
                        %addpoints(ppv1,sumiL,sumiL/1000)%paramCalibFunc(paramLHS)/1000)
                        
                    elseif body_y_pos_diff(frameind.Value) >= 0 && Ank_diffdiff(frameind.Value) < 0 && sign(Ank_diffdiff(frameind.Value)) ~= sign(Ank_diffdiff(frameind.Value-1))
                        Right_HS(frameind.Value) = 1;
                        sumiR = sum(Right_HS);
                        HS_frame = framenum.Value;
                        RHS_time = [RHS_time framenum.Value/100];
                        RHS_pos = [RHS_pos body_y_pos(frameind.Value)/1000];
                        OG_speed_right = (RHS_pos(end)-RHS_pos(end-1))/(RHS_time(end)-RHS_time(end-1));
                        
                        set(ghandle.RBeltSpeed_textbox,'String',num2str(OG_speed_right));
                        set(ghandle.Right_step_textbox,'String',num2str(sumiR))
                        
                        addpoints(ppp4,sumiR,OG_speed_right);
                        %addpoints(ppv2,sumiR,sumiR/1000)%paramCalibFunc(paramRHS)/1000)
                        
                    elseif body_y_pos_diff(frameind.Value) < 0 && Ank_diffdiff(frameind.Value) >= 0 && sign(Ank_diffdiff(frameind.Value)) ~= sign(Ank_diffdiff(frameind.Value-1))
                        Right_HS(frameind.Value) = 1;
                        sumiR = sum(Right_HS);
                        HS_frame = framenum.Value;
                        RHS_time = [RHS_time framenum.Value/100];
                        RHS_pos = [RHS_pos body_y_pos(frameind.Value)/1000];
                        OG_speed_right = abs((RHS_pos(end)-RHS_pos(end-1))/(RHS_time(end)-RHS_time(end-1)));
                        
                        set(ghandle.RBeltSpeed_textbox,'String',num2str(OG_speed_right));
                        set(ghandle.Right_step_textbox,'String',num2str(sumiR))
                        
                        addpoints(ppp4,sumiR,OG_speed_right)
                        
                    elseif body_y_pos_diff(frameind.Value) < 0 && Ank_diffdiff(frameind.Value) < 0 && sign(Ank_diffdiff(frameind.Value)) ~= sign(Ank_diffdiff(frameind.Value-1))
                        Left_HS(frameind.Value) = 1;
                        sumiL = sum(Left_HS);
                        HS_frame = framenum.Value;
                        LHS_time = [LHS_time framenum.Value/100];
                        LHS_pos = [LHS_pos body_y_pos(frameind.Value)/1000];
                        OG_speed_left = abs((LHS_pos(end)-LHS_pos(end-1))/(LHS_time(end)-LHS_time(end-1)));
                        
                        set(ghandle.LBeltSpeed_textbox,'String',num2str(OG_speed_left));
                        set(ghandle.Left_step_textbox,'String',num2str(sumiL));
                        
                        addpoints(ppp3,sumiL,OG_speed_left)
                    end
                    
                    
                elseif body_y_pos(frameind.Value) <= y_min
                    %reach one end (computer side)
                    t1 = clock;
                    inout1 = 1;
%                     flagIn = True
                    
%                 elseif body_y_pos(frameind.Value) > y_min && flagIn = True
%                     %reach one end (computer side)
%                     flagIn = False
                    
                elseif body_y_pos(frameind.Value) >= y_max
                    %reach the door side
                    t2 = clock;
                    inout2 = 1;
                    if nirs
                        inout1 = 0;
                    end
                end
                
                if inout1 == 1 && inout2 == 1
                    
                    t_diff = t1-t2;
                    walk_duration = abs((t_diff(4)*3600)+(t_diff(5)*60)+t_diff(6));
                    
                    if (~nirs) %TODO: Shuqi
                        if walk_duration < fastest
                            play(instructions('fast'));
                        elseif walk_duration > slowest
                            play(instructions('slow'));
                        else
                            play(instructions('good'));
                        end
                        inout1 = 0;
                        inout2 = 0;
                    end
                end
                
                if nirs && inout1 == 1 && inout2 == 1 %TODO: Shuqi   
                    %always stop once, only stop at the computer side, then move on to next instruction.
                    inout1 = 0; %reset
                    inout2 = 0;
                    datlog = nirsEvent('turnAndStop','R','Rest', instructions, datlog, Oxysoft, oxysoft_present);
                    pause(restDuration);

                    if currentIndex > length (eventorder)
%                         Speak(obj, 'Trial End') %TODO: in actual test, pause or say stop and rest
                        datlog = nirsEvent('relax','E','Trial End', instructions, datlog, Oxysoft, oxysoft_present);
                        currentIndex = 1; %reset the current index
                        trialIndex = trialIndex + 1;
                        if (trialIndex >  size(eventorder_all,1))
                            STOP = 1;
%                             Speak(obj, 'Experiment End') %TODO: in actual test, pause or say stop and rest
                            datlog = nirsEvent('relax','O','Experiment End', instructions, datlog, Oxysoft, oxysoft_present);
                            %TODO: should automatically send signal to terminate the software
                        else  %Restart a trial
                            eventorder = eventorder_all(trialIndex,:); 
                            %0 = stand and alphabet, 1 = walk and alphabet, 2 = walk, 
                            %start the new trial with a rest
                            datlog = nirsEvent('rest','R','Rest', instructions, datlog, Oxysoft, oxysoft_present);
                            pause(restDuration);
                        end
                    end
                    disp('next current index');
                    disp(eventorder(currentIndex));
                    if STOP ~= 1&& eventorder(currentIndex)== 0 %stand and alphabet
                        datlog = nirsEvent('standAlphabet','A',['Stand and Alphabet' alphabetLetter], instructions, datlog, Oxysoft, oxysoft_present);
                        pause(restDuration)

                        %complete follow a rest block
                        datlog = nirsEvent('stopAndRest','R','Rest', instructions, datlog, Oxysoft, oxysoft_present);
                        pause(restDuration);
                        currentIndex = currentIndex + 1;

                        if currentIndex > length (eventorder)
                            datlog = nirsEvent('relax','E','Trial End', instructions, datlog, Oxysoft, oxysoft_present);
                            currentIndex = 1; %reset the current index
                            trialIndex = trialIndex + 1;
                            if (trialIndex >  size(eventorder_all,1))
                                STOP = 1;
                                datlog = nirsEvent('relax','O','Experiment End', instructions, datlog, Oxysoft, oxysoft_present); 
                            else
                                eventorder = eventorder_all(trialIndex,:); 
                                datlog = nirsEvent('rest','R','Rest', instructions, datlog, Oxysoft, oxysoft_present);
                                pause(restDuration);
                            end      
                        end 
                    end

                    if STOP ~= 1 && eventorder(currentIndex)==1 %first event is walk and alphabet
                        datlog = nirsEvent('walkAlphabet','B',['Walk and Alphabet' alphabetLetter], instructions, datlog, Oxysoft, oxysoft_present);
                        currentIndex  = currentIndex + 1;
                    elseif STOP ~= 1 && eventorder(currentIndex) == 2 %first event is walk
                        datlog = nirsEvent('walk','W','Walk', instructions, datlog, Oxysoft, oxysoft_present);
                        currentIndex = currentIndex + 1;
                    end
                end             
            end            
        end
    end %While, when STOP button is pressed
    if STOP
        datlog.messages{end+1} = ['Stop button pressed at: ' num2str(now) ' ,stopping... '];
        %     log=['Stop button pressed, stopping... ' num2str(clock)];
        %     listbox{end+1}=log;
        disp(['Stop button pressed, stopping... ' num2str(clock)]);
        set(ghandle.Status_textbox,'String','Stopping...');
        set(ghandle.Status_textbox,'BackgroundColor','red');
    else
    end
catch ME
    datlog.errormsgs{end+1} = 'Error ocurred during the control loop';
    datlog.errormsgs{end+1} = ME;
    disp(ME)
    %disp(ME.stack)
    ME.getReport
    %     log=['Error ocurred during the control loop'];%End try
    %     listbox{end+1}=log;
    disp('Error ocurred during the control loop, see datlog for details...');
end
%% Closing routine
%End communications
global addLog
try
    aux = cellfun(@(x) (x-datlog.framenumbers.data(1,2))*86400,addLog.keypress(:,2),'UniformOutput',false);
    addLog.keypress(:,2)=aux;
    addLog.keypress=addLog.keypress(cellfun(@(x) ~isempty(x),addLog.keypress(:,1)),:); %Eliminating empty entries
    datlog.addLog=addLog;
catch ME
    ME
end
try
    save(savename,'datlog');
catch ME
    disp(ME);
end

% try %stopping the treadmill
%     %see if the treadmill is supposed to stop at the end of the profile
%     if get(ghandle.StoptreadmillEND_checkbox,'Value')==1 && STOP ~=1
%         set(ghandle.Status_textbox,'String','Stopping...');
%         set(ghandle.Status_textbox,'BackgroundColor','red');
%         pause(1)%provide a little time to collect the last steps and so forth
%         smoothStop(t);
%         %see if the treadmill should be stopped when the STOP button is pressed
%     elseif get(ghandle.StoptreadmillSTOP_checkbox,'Value')==1 && STOP == 1
%         
%         set(ghandle.Status_textbox,'String','Stopping');
%         set(ghandle.Status_textbox,'BackgroundColor','red');
%         pause(1)
%         smoothStop(t);
%     end
%     
% catch ME
%     datlog.errormsgs{end+1} = 'Error stopping the treadmill';
%     datlog.errormsgs{end+1} = ME;
% end

% pause(1)
disp('closing comms');
try
%     closeNexusIface(MyClient);
%     closeTreadmillComm(t);
    %     keyboard
    disp('Closing Nexus Client') %Todo: Shuqi
    disp(clock);
    closeNexusIface(MyClient);
    closeTreadmillComm(t);
    disp('Done Closing'); %TODO: Shuqi
    disp(clock);
catch ME
    datlog.errormsgs{end+1} = ['Error ocurred when closing communications with Nexus & Treadmill at ' num2str(clock)];
    datlog.errormsgs{end+1} = ME;
    %     log=['Error ocurred when closing communications with Nexus & Treadmill (maybe they were not open?) ' num2str(clock)];
    %     listbox{end+1}=log;
    disp(['Error ocurred when closing communications with Nexus & Treadmill, see datlog for details ' num2str(clock)]);
    disp(ME);
end

disp('converting time in datlog...');
%convert time data into clock format then re-save
datlog.buildtime = datestr(datlog.buildtime);

%convert frame times & markers
temp = find(datlog.framenumbers.data(:,1)==0,1,'first');
datlog.framenumbers.data(temp:end,:) = [];
for z = 1:temp-1
    datlog.framenumbers.data(z,3) = etime(datevec(datlog.framenumbers.data(z,2)),datevec(datlog.framenumbers.data(1,2)));
end
%save marker positions for the SAME frames
%  for l=1:length(mn)
%      eval(['aux=datlog.Markers.' mn{l} ';']);
%      aux(temp:end,:)=[];
%      eval(['datlog.Markers.' mn{l} '=aux;']);
%      eval(['aux=datlog.Markers.' mn{l} 'filt;']);
%      aux(temp:end,:)=[];
%      eval(['datlog.Markers.' mn{l} 'filt=aux;']);
%      eval(['aux=datlog.Markers.' mn{l} 'var;']);
%      aux(temp:end,:)=[];
%      eval(['datlog.Markers.' mn{l} 'var=aux;']);
%  end

%convert RHS times
temp = find(datlog.stepdata.RHSdata(:,1) == 0,1,'first');
datlog.stepdata.RHSdata(temp:end,:) = [];
datlog.stepdata.paramRHS(temp:end,:) = [];
for z = 1:temp-1
    datlog.stepdata.RHSdata(z,4) = etime(datevec(datlog.stepdata.RHSdata(z,2)),datevec(datlog.framenumbers.data(1,2)));
end
%convert LHS times
temp = find(datlog.stepdata.LHSdata(:,1) == 0,1,'first');
datlog.stepdata.LHSdata(temp:end,:) = [];
datlog.stepdata.paramLHS(temp:end,:) = [];
for z = 1:temp-1
    datlog.stepdata.LHSdata(z,4) = etime(datevec(datlog.stepdata.LHSdata(z,2)),datevec(datlog.framenumbers.data(1,2)));
end
%convert RTO times
temp = find(datlog.stepdata.RTOdata(:,1) == 0,1,'first');
datlog.stepdata.RTOdata(temp:end,:) = [];
for z = 1:temp-1
    datlog.stepdata.RTOdata(z,4) = etime(datevec(datlog.stepdata.RTOdata(z,2)),datevec(datlog.framenumbers.data(1,2)));
end
%convert LTO times
temp = find(datlog.stepdata.LTOdata(:,1) == 0,1,'first');
datlog.stepdata.LTOdata(temp:end,:) = [];
for z = 1:temp-1
    datlog.stepdata.LTOdata(z,4) = etime(datevec(datlog.stepdata.LTOdata(z,2)),datevec(datlog.framenumbers.data(1,2)));
end

%convert command times
% temp = all(isnan(datlog.TreadmillCommands.read(:,4)),2);
% datlog.TreadmillCommands.read=datlog.TreadmillCommands.read(~temp,:);
% for z = 1:size(datlog.TreadmillCommands.read,1)
%     datlog.TreadmillCommands.read(z,4) = etime(datevec(datlog.TreadmillCommands.read(z,4)),datevec(datlog.framenumbers.data(1,2))); %This fails when no frames were received
% end

% temp = all(isnan(datlog.TreadmillCommands.sent(:,4)),2);
% datlog.TreadmillCommands.sent=datlog.TreadmillCommands.sent(~temp,:);
% for z = 1:size(datlog.TreadmillCommands.sent,1)
%     datlog.TreadmillCommands.sent(z,4) = etime(datevec(datlog.TreadmillCommands.sent(z,4)),datevec(datlog.framenumbers.data(1,2)));
% end

%convert audio times
datlog.audioCues.start = datlog.audioCues.start';
datlog.audioCues.audio_instruction_message = datlog.audioCues.audio_instruction_message';
datlog.audioCues.stop = datlog.audioCues.stop';
temp = isnan(datlog.audioCues.start);
disp('\nConverting datalog, current starts \n'); %TODo: Shuqi
disp(datlog.audioCues.start);
datlog.audioCues.start=datlog.audioCues.start(~temp);
datlog.audioCues.start = ((datlog.audioCues.start)-(datlog.framenumbers.data(1,2)))*86400;

temp = isnan(datlog.audioCues.stop);
datlog.audioCues.stop=datlog.audioCues.stop(~temp);
datlog.audioCues.stop = ((datlog.audioCues.stop) - (datlog.framenumbers.data(1,2)))*86400;


%Get rid of graphical objects we no longer need:
if exist('ff','var') && isvalid(ff)
    close(ff)
end
set(ghandle.profileaxes,'Color',[1,1,1])
delete(findobj(ghandle.profileaxes,'Type','Text'))
delete(findobj(ghandle.profileaxes,'Type','Patch'))

%Print some info:
NR=length(datlog.stepdata.paramRHS);
NL=length(datlog.stepdata.paramLHS);
for M=[20,50]
    i1=max([1 NR-M]);
    i2=max([1 NL-M]);
    disp(['Average params for last ' num2str(M) ' strides: '])
    disp(['R=' num2str(nanmean(datlog.stepdata.paramRHS(i1:end))) ', L=' num2str(nanmean(datlog.stepdata.paramLHS(i2:end)))])
end

disp('saving datlog...');
try
    save(savename,'datlog');
    [d,~,~]=fileparts(which(mfilename));
    save([d '\..\datlogs\lastDatlog.mat'],'datlog');
catch ME
    disp(ME);
end
