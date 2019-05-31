%% StimParamSweep_MatlabControl_amplitude
% Written: 4/17/17
% Edited: 10/19/17 to interleave staircases randomly rather than every
% other
% J Cronin

%% This script should be used with the TDT Project: StimParamSweep_v2
%% This code specifically adjusts the current amplitude
% This runs a double interleaved staircase with _-down _-up
% TDT will deliver stimulation, record the stim parameters used, and
% record 'yes' and 'no' depending on if subject felt the stim.
% Matlab code will get the 'yes' and 'no' info and set the next stim level
% based on the previous trials

% This MATLAB script will also save the stim parameters and the responses
% from the TDT and print figures right afterwards of the staircase results
% and the psychometric curve using Geoff Boynton's code.

% Experiment details (like up/down numbers) and number of trials are set in
% TDT, except for a few below
clc; clear all; close all
% dbstop if error

%% Experimental parameters set in MATLAB (others in TDT):
% Step size (uA): first value will be used for the first 4 steps, the 2nd
% value for the next 4 steps, and the last value for all of the subsequent
% steps until the alloted number of trials is reached
stepSize = [250, 100, 50]; % in uA
% stepSize = [125, 100, 50]; % in uA
% stepSize = [100, 100, 50]; % in uA
numForStepSize = [4, 8]; % Can adjust the number of each step sizes here, the 3rd step size will fill in the remainder of trials

%% Open connection with TDT and begin program
DA = actxcontrol('TDevAcc.X');
DA.ConnectServer('Local'); %initiates a connection with an OpenWorkbench server. The connection adds a client to the server
pause(1)

while DA.CheckServerConnection ~= 1
    disp('OpenWorkbench is not connected to server. Trying again...')
    close all
    DA = actxcontrol('TDevAcc.X');
    DA.ConnectServer('Local');
    pause(1) % seconds
end
clc
disp('Connected to server')
disp('Remember: if you want to end recording early, end TDT recording and Matlab will automatically save param sweep vars')

% If OpenWorkbench is not in Record mode, then this will set it to record
if DA.GetSysMode ~= 3
    DA.SetSysMode(3);
    while DA.GetSysMode ~= 3
        pause(.1)
    end
end

% Check the experiment type that was set on the TDT. If it's not train
% duration then throw error to ensure we run the intended Matlab script.
% amplitude = 1
% train duration = 0
experimentType = DA.GetTargetVal('RZ5D.expType');
if experimentType~=1
    warning('You are running the amplitude Matlab script, but Amp experiment is not selected on TDT, switching experiment type to Amp.')
    DA.SetTargetVal('RZ5D.expType',1);
else
    disp('Pulse amplitude script')
end

%% Initizlize values
% Check that you're ready to initialize values
respInitialize = questdlg('Are you ready for Matlab to read initialized values from TDT?');
while strcmp(respInitialize, 'No')
    % wait and then try again
    pause(10)
    respInitialize = questdlg('Are you ready for Matlab to read initialized values from TDT?');
end
if strcmp(respInitialize, 'Yes')
    tank = DA.GetTankName;
    noNum = DA.GetTargetVal('RZ5D.upNum'); % number of no's needed to increase intensity
    yesNum = DA.GetTargetVal('RZ5D.downNum'); % number of yes's needed to decrease intensity
    
    numTrials = DA.GetTargetVal('RZ5D.numTrials');
    numCatchTrials = DA.GetTargetVal('RZ5D.numCatchTrials'); % this is included in the total number of trials (i.e., if numTrials=60 and numCatchTrials=10, then normal trials=50)
    
    % Since this code is changing the amplitude, use the starting
    % train duration to set the PTD which will be used throughout
    % experiment:
    PTD_setting = DA.GetTargetVal('RZ5D.startingTD');
    DA.SetTargetVal('RZ5D.PTD', PTD_setting);
    
    % Read in a few more stim parameters from the TDT to save later
    PW = DA.GetTargetVal('RZ5D.cPW(us)');
    PulseFreq = DA.GetTargetVal('RZ5D.PulseFreqCalc');
    PTD = DA.GetTargetVal('RZ5D.PTD_to_use(ms)');
    
    % If the number of trials or catch trials isn't even, then make it even so
    % it's nicely divisible by 2 below
    if mod(numTrials, 2)~=0
        numTrials = numTrials + 1;
    end
    if mod(numCatchTrials, 2)~=0
        numCatchTrials = numCatchTrials + 1;
    end
    % Randomize the order of the stairs
    stairOrder = [ones(numTrials/2,1); 2*ones(numTrials/2,1)];
    stairOrder = stairOrder(randperm(numTrials));
    
    % Display the values that were just read from the TDT
    disp([num2str(yesNum), '-down, ', num2str(noNum), '-up; ',...
        num2str(numTrials), ' trials with ', num2str(numCatchTrials), ' catch trials'])
    
    % Randomly distribute the catch trials into all of the trials starting
    % after the 5th trial (want to get into a range where subject may be
    % feeling something before running catch trials
    catchTrial = zeros(numTrials/2-5, 2);
    catchTrial(1:numCatchTrials/2,:) = 1;
    catchTrial(:,1) = catchTrial(randperm(length(catchTrial)),1);
    catchTrial(:,2) = catchTrial(randperm(length(catchTrial)),2);
    catchTrial = [zeros(5,2); catchTrial]; % add those 5 non-catch trials back to the beginning of each staircase
    
    % Using the number of trials, establish a step size vector for each trial:
    stepSizeVec = zeros((numTrials/2)-1, 2); % Note - will pad the first row with 0's at the end
    stepSizeVec(1:numForStepSize(1),:) = stepSize(1);
    stepSizeVec(numForStepSize(1)+1:end,:) = stepSize(2);
    % Account for the location of the catch trials, so that you still use the
    % complete order of the step size progression (i.e., don't skip some of the
    % stepSize(2)'s just because there are catch trials in there
    stepSizeVec(logical(catchTrial(1:end-1,:))) = 0;
    ind = max(find(stepSizeVec(:,1)==stepSize(2),numForStepSize(2)));
    stepSizeVec(ind+1:end,1) = stepSize(3);
    ind = max(find(stepSizeVec(:,2)==stepSize(2),numForStepSize(2)));
    stepSizeVec(ind+1:end,2) = stepSize(3);
    stepSizeVec(logical(catchTrial(1:end-1,:))) = 0;
    stepSizeVec = [[0 0]; stepSizeVec]; % Adding 0's in the first row because the starting amp is set, and we'll never use the first row of step sizes

    
    % Initialize the staircases with zeros for now
    pulseAmp = zeros(numTrials, 1);
    startingPulseAmp = DA.GetTargetVal('RZ5D.startingPulseAmp');
    disp(['Starting pulse amplitude will be: ', num2str(startingPulseAmp), ' uA, with a train duration of ', num2str(PTD), ' ms'])
    pulseAmp(1) = startingPulseAmp;
    
    staircase.intensity = zeros(numTrials/2,2);
    staircase.intensity(1,:) = startingPulseAmp;
    
    intensity = zeros(numTrials, 1);
    intensity(1) = startingPulseAmp;
    
    % Responses are stored as 0: no, and 1: yes
    response = zeros(numTrials, 1);
    % We'll also store responses as a matrix separated by the two staircases
    staircase.response = zeros(numTrials/2,2);
    
elseif strcmp(respInitialize, 'Cancel')
    warning('Must initialize Matlab and TDT values')
    return
end

%% Before proceeding make sure that the system is armed:
if DA.GetTargetVal('RZ5D.IsArmed') == 0
    disp('System is not armed');
elseif DA.GetTargetVal('RZ5D.IsArmed') == 1
    disp('System armed');
end

while DA.GetTargetVal('RZ5D.IsArmed')~=1
    % waiting until the system is armed
end
disp(['amp = ', num2str(pulseAmp(1))])
% Once the system is armed, set the first TDT PulseAmp to the StartingPulseAmp value
pause(1) % make sure there's been enough time to get the first value
DA.SetTargetVal('RZ5D.PulseAmp', pulseAmp(1));
disp('System armed');

%% Run double interleaved staircase
% Run the staircase
nextIntensity_afterCatch = [0 0]; % this will save intensity values to use after a catch trial for the first and 2nd staricases respectively

yesInaRow = [0 0]; % this will track how many yes's there were in a row if you have more than a 1-down requirement for each staircase
noInaRow = [0 0];

stairOrderMatrix = [find(stairOrder==1), find(stairOrder==2)]; % this will be used to track which staircase each trial number is in and save responses and amps accordingly

% Turn the ready indicator on and wait for the stim button to be pressed
% for the first time
DA.SetTargetVal('RZ5D.readyIndicator', 1); % set ready indicator to 0 for computation
while DA.GetTargetVal('RZ5D.StimButton')~=1
    % waiting...
end

DA.SetTargetVal('RZ5D.readyIndicator', 0); % set ready indicator to 0 for computation
pause(1.2) % pause for 1.2 seconds while TDT delivers the stimulation (1 second for the audio cue pause, and 200 ms for approximate TD)
DA.SetTargetVal('RZ5D.readyForAnswerIndicator', 1); % set the indicator for answering to on, so experimenter knows that stim (or catch) has been delivered and system is waiting for response

for i=1:numTrials
    whichStair = stairOrder(i); % this will switch back and forth between the two staircases
    doubleMatPos = find(stairOrderMatrix==i); % This is the linear index within the matrices that we're using to track responses, intensity, catch trials, step sizes, etc.
    
    while DA.GetTargetVal('RZ5D.yesButton')~=1 && DA.GetTargetVal('RZ5D.noButton')~=1
        % waiting for a yes/no button to be pressed, and also checking that
        % the system is still in record mode
        % If the TDT has been taken out of record mode, then we'll end the
        % experiment early and save the variables:
        if DA.GetSysMode ~= 3
            disp('TDT Recording was ended early, ending and saving Matlab script now')
            Save_ParamSweepVariables
            return
        end
    end
    
    if DA.GetTargetVal('RZ5D.yesButton')==1
        resp = 'Yes';
        displayTrial = ['Trial #', num2str(i), ': Yes'];
        DA.SetTargetVal('RZ5D.yesButton', 0);
        DA.SetTargetVal('RZ5D.yesIndicator', 1); % set the indicator to 1
    elseif DA.GetTargetVal('RZ5D.noButton')==1
        resp = 'No';
        displayTrial = ['Trial #', num2str(i), ': No'];
        DA.SetTargetVal('RZ5D.noButton', 0);
        DA.SetTargetVal('RZ5D.noIndicator', 1); % set the indicator to 1
    else
        warning('Do not know which button was pressed')
    end
    
    % Now we assume that either the yesButton or the noButton was pressed
    % and use the 'yes' 'no' indicator targets in TDT to figure out which
    % one
    % Display what was pressed and what intensity was just delivered:
    if catchTrial(doubleMatPos) == 1
        disp([displayTrial, ', Intensity = ', num2str(pulseAmp(i)), ', Catch'])
    else
        disp([displayTrial, ', Intensity = ', num2str(pulseAmp(i))])
    end
    
    % First turn off the ready to answer indicator, since an answer was just
    % given:
    DA.SetTargetVal('RZ5D.readyForAnswerIndicator', 0);
    
    if strcmp(resp,'Yes') % yes was pressed, so step the intensity down
        % This pulls what was just answered/delivered into the final
        % response and intensity vectors
        response(i) = 1;
        staircase.response(doubleMatPos) = 1;
        intensity(i) = DA.GetTargetVal('RZ5D.PulseAmp'); % record the last pulse amp that was used
        
        if ~catchTrial(doubleMatPos)
            yesInaRow(whichStair) = yesInaRow(whichStair) + 1; % only increment the number of answers in a row when not a catch trial
            noInaRow(whichStair) = 0; % reset the no count
        end
        
        % set the next intensity value in the staircase matrix which tracks
        % both of the two staircases (one in each column)
        [row, col] = ind2sub(size(stairOrderMatrix), doubleMatPos);
        if row < numTrials/2
            if yesInaRow(whichStair) == yesNum
                staircase.intensity(doubleMatPos+1) = staircase.intensity(doubleMatPos) - stepSizeVec(doubleMatPos+1);
                yesInaRow(whichStair) = 0; % reset the count
            else % intensity doesn't change
                staircase.intensity(doubleMatPos+1) = staircase.intensity(doubleMatPos);
            end
        end
        
    elseif strcmp(resp,'No') % no was pressed, so step the intensity up
        response(i) = 0;
        staircase.response(doubleMatPos) = 0;
        intensity(i) = DA.GetTargetVal('RZ5D.PulseAmp'); % record the last pulse amp that was used
        
        if ~catchTrial(doubleMatPos)
            noInaRow(whichStair) = noInaRow(whichStair) + 1; % only increment the number of answers in a row when not a catch trial
            yesInaRow(whichStair) = 0; % reset the yes count
        end
        
        % set the next intensity value in the staircase matrix which tracks
        % both of the two staircases (one in each column)
        [row, col] = ind2sub(size(stairOrderMatrix), doubleMatPos);
        if row < numTrials/2
            if noInaRow(whichStair) == noNum
                staircase.intensity(doubleMatPos+1) = staircase.intensity(doubleMatPos) + stepSizeVec(doubleMatPos+1);       
                noInaRow(whichStair) = 0; % reset the count
            else % intensity doesn't change
                staircase.intensity(doubleMatPos+1) = staircase.intensity(doubleMatPos);
            end
        end 
        
    else
        warning('Something went wrong in going thru yes/no steps')
        return
    end
    
    if row < numTrials/2 && catchTrial(doubleMatPos+1)==1
        % This will just overwrite the next staircase value if it actually
        % needs to be a catch trial, but will first save the next value to
        % use in a varible that can be accessed after the catch trial;
        % however, should only save the next value when you're doing a
        % catch trial that was proceeded by a real trial. If you've got a
        % few catches in a row, then you should only save the next value to
        % use once
        if catchTrial(doubleMatPos)~=1
            nextIntensity_afterCatch(whichStair) = staircase.intensity(doubleMatPos+1);
        end
        staircase.intensity(doubleMatPos+1) = 0; % Set the next amp in the staircase to 0 if it's going to be a catch trial
        
    elseif row < numTrials/2 && catchTrial(doubleMatPos)==1
        % If the current trial was a catch trial, then rather then taking that 0
        % amplitude intensity and adding the step to it, we should take the
        % saved next value to use
        staircase.intensity(doubleMatPos+1) =  nextIntensity_afterCatch(whichStair);
    end
    
    % Now that new values are set, update the next pulse amp that will be
    % used
    if i<numTrials
        pulseAmp(i+1) = staircase.intensity(find(stairOrderMatrix==i+1));
        DA.SetTargetVal('RZ5D.PulseAmp', pulseAmp(i+1)); % send new pulseAmp to TDT

    end
          
    % Once all of the computation is done and new values are set, turn
    % the ready indicator on (just so the experimenter knows)
    DA.SetTargetVal('RZ5D.readyIndicator', 1); % set the indicator to 1 (LED on) so experimenter knows that the Stim button can be pressed again
    % and turn the catch trial off
    DA.SetTargetVal('RZ5D.catchTrialIndicator', 0); % Turn the catch trial indicator LED off since the catch trial should have just been delivered
    
    % Wait until stim button is pressed and then turn off the Yes and No
    % indicators to indicate to experimentor that they can press a Yes/No
    % button again
    % but only do this for every trial before the last trial
    if i<numTrials
        while DA.GetTargetVal('RZ5D.StimButton')~=1
            % waiting for a stim button to be pressed, and also checking that
            % the system is still in record mode
            % If the TDT has been taken out of record mode, then we'll end the
            % experiment early and save the variables:
            if DA.GetSysMode ~= 3
                disp('TDT Recording was ended early, ending and saving Matlab script now')
                Save_ParamSweepVariables
                return
            end
        end
        DA.SetTargetVal('RZ5D.noIndicator', 0); % set the indicator back to 0
        DA.SetTargetVal('RZ5D.yesIndicator', 0); % set the indicator back to 0
        if catchTrial(find(stairOrderMatrix==i+1))==1
            DA.SetTargetVal('RZ5D.catchTrialIndicator', 1); % Turn the catch indicator on if the next trial will be a catch trial
        end
        %     DA.SetTargetVal('RZ5D.catchTrialIndicator', 0); % Turn the catch trial indicator LED off since the catch trial should have just been delivered
        DA.SetTargetVal('RZ5D.readyIndicator', 0); % set ready indicator to 0 for computation
        pause(1.2) % pause for 1.2 seconds while TDT delivers the stimulation (1 second for the audio cue pause, and 200 ms for approximate TD)

        DA.SetTargetVal('RZ5D.readyForAnswerIndicator', 1); % set the indicator for answering to on, so experimenter knows that stim (or catch) has been delivered and system is waiting for response
    end 
    
end

% Now that trials are done, turn all indicators off:
DA.SetTargetVal('RZ5D.noIndicator', 0); 
DA.SetTargetVal('RZ5D.yesIndicator', 0); 
DA.SetTargetVal('RZ5D.catchTrialIndicator', 0);
DA.SetTargetVal('RZ5D.readyForAnswerIndicator', 0);
DA.SetTargetVal('RZ5D.readyIndicator', 0);

%% When run is ended, close the connection

% Disarm stim:
DA.SetTargetVal('RZ5D.ArmSystem', 0);

% Read the loaded circuit's name so that we can save this
circuitLoaded = DA.GetDeviceRCO('RZ5D');

% Close ActiveX connection:
DA.CloseConnection
if DA.CheckServerConnection == 0
    disp('Server was disconnected');
end

%% Save
Save_ParamSweepVariables

%% Run code from Geoff to plot psychophysics curves
% Put data into format needed for Geoff's code
results.intensity = intensity/1000;
results.response = response;

% Add path to other supporting files from Geoff
[pathstr,name,~] = fileparts(pwd); % check what folder your in (will work on different computers)
if ispc
    addpath([pathstr, '\', name, '\psychometric_curve_from_Geoff'])
elseif ismac
    addpath([pathstr, '/', name, '/psychometric_curve_from_Geoff'])
end

% Set the guess rate which is estimated as the number of false positives
% that the subject gave (i.e., the proportion of time that they responded
% 'yes' to a catch trial that didn't actually have stim
guessRate = sum(results.response(results.intensity==0))/sum(results.intensity==0);
p.g = guessRate;
FitECogData_noFileLoad
figure(1), title('Staircases Merged with Catch Trials')

%% Also plot the staircases separately and without the catch trials
stair1.intensity = staircase.intensity(staircase.intensity(:,1)~=0,1)/1000;
stair1.response = staircase.response(staircase.intensity(:,1)~=0,1);
stair2.intensity = staircase.intensity(staircase.intensity(:,2)~=0,2)/1000;
stair2.response = staircase.response(staircase.intensity(:,2)~=0,2);

randomlyInterleaved = 1; % set whether or not the files were randomly interleaved
% if 0, then assumes the staircases went every other trial
figure(3)
clf
if randomlyInterleaved == 0
    x = 1:2:length(stair1.intensity)*2-1;
elseif randomlyInterleaved == 1
    x = stairOrderMatrix(~catchTrial(:,1),1);
else
    error('randomlyInterleaved variable must be set to 0 or 1')
end
stairs(x, stair1.intensity,'k--', 'Linewidth', 1.5);
hold on
id = stair1.response == 1;
plot(x(id),stair1.intensity(id),'^', 'MarkerFaceColor', [0 191 255]/255,...
    'MarkerEdgeColor', [0 191 255]/255);
plot(x(~id),stair1.intensity(~id),'x', 'MarkerFaceColor', [0 0 205]/255,...
    'MarkerEdgeColor', [0 0 205]/255, 'Linewidth', 2, 'MarkerSize', 8);

if randomlyInterleaved == 0
    x = 2:2:length(stair1.intensity)*2;
elseif randomlyInterleaved == 1
    x = stairOrderMatrix(~catchTrial(:,2),2);
else
    error('randomlyInterleaved variable must be set to 0 or 1')
end
stairs(x, stair2.intensity, 'k:', 'Linewidth', 1.5);
id = stair2.response == 1;
plot(x(id),stair2.intensity(id),'v', 'MarkerFaceColor', [255 153 51]/255,...
    'MarkerEdgeColor', [255 153 51]/255);
plot(x(~id),stair2.intensity(~id), '+', 'MarkerFaceColor', [255 69 0]/255,...
    'MarkerEdgeColor', [255 69 0]/255, 'Linewidth', 2, 'MarkerSize', 8);

xlabel('Trial Number');
ylabel('Current Amplitude (mA)');
title('Double-interleaved Staircase');

legend('3-down, 1-up Staircase 1', 'Staircase 1: Yes', 'Staircase 1: No',...
    '3-down, 1-up Staircase 2', 'Staircase 2: Yes', 'Staircase 2: No')
legend('Location', 'southeast')
grid minor

%% Plot the staircases separately
% Uncomment the below if we want to plot staircases separately 
 
% % Compare staircase 1
% figure(4)
% subplot(2,1,1)
% stairs(staircase.intensity(:,1));
% hold on
% x = 1:length(staircase.intensity(:,1));
% id = staircase.response(:,1) == 1;
% plot(x(id),staircase.intensity(id,1),'ko','MarkerFaceColor','g');
% plot(x(~id),staircase.intensity(~id,1),'ko','MarkerFaceColor','r');
% xlabel('Trial Number');
% ylabel('Amplitude (mA)');
% title('Staircase 1 with catch trials');
% 
% subplot(2,1,2)
% stairs(stair1.intensity);
% hold on
% x = 1:length(stair1.intensity);
% id = stair1.response == 1;
% plot(x(id),stair1.intensity(id),'ko','MarkerFaceColor','g');
% plot(x(~id),stair1.intensity(~id),'ko','MarkerFaceColor','r');
% xlabel('Trial Number');
% ylabel('Amplitude (mA)');
% title('Without catch trials');
% 
% 
% % Compare staircase 2
% figure(5)
% subplot(2,1,1)
% stairs(staircase.intensity(:,2));
% hold on
% x = 1:length(staircase.intensity(:,2));
% id = staircase.response(:,2) == 1;
% plot(x(id),staircase.intensity(id,2),'ko','MarkerFaceColor','g');
% plot(x(~id),staircase.intensity(~id,2),'ko','MarkerFaceColor','r');
% xlabel('Trial Number');
% ylabel('Amplitude (mA)');
% title('Staircase 2 with catch trials');
% 
% subplot(2,1,2)
% stairs(stair2.intensity);
% hold on
% x = 1:length(stair2.intensity);
% id = stair2.response == 1;
% plot(x(id),stair2.intensity(id),'ko','MarkerFaceColor','g');
% plot(x(~id),stair2.intensity(~id),'ko','MarkerFaceColor','r');
% xlabel('Trial Number');
% ylabel('Amplitude (mA)');
% title('Without catch trials');