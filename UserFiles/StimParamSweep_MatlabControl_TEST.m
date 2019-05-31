%% StimParamSweep_MatlabControl
% Written: 4/17/17
% J Cronin

%% This script should be used with the TDT Project: StimParamSweep_v2 
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

close all; clear all; clc 
%% Experimental parameters set in MATLAB (others in TDT):
% Step size (uA): first value will be used for the first 4 steps, the 2nd
% value for the next 8 steps, and the last value for all of the subsequent
% steps until the alloted number of trials is reached
stepSize = [250, 100, 50]; % in uA 
numForStepSize = [4, 8]; % Can adjust the number of each step sizes here, the 3rd step size will fill in the remainder of trials

%% Initizlize values
noNum = 1; % number of no's needed to increase intensity
yesNum = 3; % number of yes's needed to decrease intensity

numTrials = 20;
numCatchTrials = 4; % this is included in the total number of trials (i.e., if numTrials=60 and numCatchTrials=10, then normal trials=50)

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
ind = max(find(stepSizeVec(:,1)==100,numForStepSize(2)));
stepSizeVec(ind+1:end,1) = stepSize(3);
ind = max(find(stepSizeVec(:,2)==100,numForStepSize(2)));
stepSizeVec(ind+1:end,2) = stepSize(3);
stepSizeVec(logical(catchTrial(1:end-1,:))) = 0;
stepSizeVec = [[0 0]; stepSizeVec]; % Adding 0's in the first row because the starting amp is set, and we'll never use the first row of step sizes

% Initialize the staircases with zeros for now
pulseAmp = zeros(numTrials, 1);
startingPulseAmp = 250;
pulseAmp(1) = startingPulseAmp;

staircase.intensity = ones(numTrials/2,2); %wouldn't actually use 1 uA, but this way can tell if something goes wrong
staircase.intensity(1,:) = startingPulseAmp;

intensity = zeros(numTrials, 1);
intensity(1) = startingPulseAmp;

% Responses are stored as 0: no, and 1: yes
response = zeros(numTrials, 1);
% We'll also store responses as a matrix separated by the two staircases
staircase.response = 2*ones(numTrials/2,2); % setting it to the number 2 so it's easy to tell if something goes wrong

%% Set up psychtoolbox
sca; % close any open psychtoolbox screens

PsychDefaultSetup(2);
screenNumber = max(Screen('Screens'));

% Define colors
black = BlackIndex(screenNumber);
white = WhiteIndex(screenNumber);

% % Probably need to remove this line!!!!!!!!!!!!!!!!!
% Screen('Preference', 'SkipSyncTests', 0);

% Open window screen and get size/coordinates
[window, windowRect] = PsychImaging('OpenWindow', screenNumber, black);
[xCenter, yCenter] = RectCenter(windowRect);
[screenXpixels, screenYpixels] = Screen('WindowSize', window);

% Query the frame duration
ifi = Screen('GetFlipInterval', window);

% Set up alpha-blending for smooth (anti-aliased) lines
Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

% Setup the text type for the window
Screen('TextFont', window, 'Ariel');

% The avaliable keys to press
escapeKey = KbName('ESCAPE');
leftKey = KbName('LeftArrow');
rightKey = KbName('RightArrow');

%% Start psychtoolbox screen
% Set fixation cross
fixCrossDimPix = 40;
% Now we set the coordinates (these are all relative to zero we will let
% the drawing routine center the cross in the center of our monitor for us)
xCoords = [-fixCrossDimPix fixCrossDimPix 0 0];
yCoords = [0 0 -fixCrossDimPix fixCrossDimPix];
fixCrossCoords = [xCoords; yCoords];
% Set the line width for our fixation cross
lineWidthPix = 4;
% Draw the fixation cross in white, set it to the center of our screen and
% set good quality antialiasing
Screen('DrawLines', window, fixCrossCoords, lineWidthPix, white, [xCenter yCenter], 2);

% Flip to the screen
vbl = Screen('Flip', window);

% Retreive the maximum priority number
topPriorityLevel = MaxPriority(window);

    % Length of time and number of frames we will use for each drawing test
numSecs = 2;
numFrames = round(numSecs / ifi);
waitframes = 1;



%% Run double interleaved staircase with PsychToolbox
% Run the staircase

% stairIter = repelem(1:numTrials/2,2); % this will count up by one every two trials (it's basically the row that should be accessed in the staircase matrix
nextIntensity_afterCatch = [0 0]; % this will save intensity values to use after a catch trial for the first and 2nd staricases respectively

yesInaRow = [0 0]; % this will track how many yes's there were in a row if you have more than a 1-down requirement for each staircase
noInaRow = [0 0];

stairOrderMatrix = [find(stairOrder==1), find(stairOrder==2)]; % this will be used to track which staircase each trial number is in and save responses and amps accordingly

for i=1:numTrials
    whichStair = stairOrder(i); % this will switch back and forth between the two staircases
    
    % The fixation cross was just displayed, so pause for stim delivery
    % before asking if subject felt the stim
    
    % Add text and yes/no buttons to the window
    Screen('TextSize', window, 56);
    DrawFormattedText(window, 'Did you feel the stimulation?', 'center', 'center', white);
    Screen('TextSize', window, 56);
    DrawFormattedText(window, 'Yes', screenXpixels*.25, screenYpixels*.75, white);
    DrawFormattedText(window, 'No', screenXpixels*.75, screenYpixels*.75, white);
    Screen('TextSize', window, 30);
    DrawFormattedText(window, 'press left arrow key', screenXpixels*.25-125, screenYpixels*.75+50, white);
    DrawFormattedText(window, 'press right arrow key', screenXpixels*.75-125, screenYpixels*.75+50, white);
    
    % draw left arrow
    % create a triangle
    head   = [screenXpixels*.25+25,  screenYpixels*.75+100]; % coordinates of head
    width  = 25;           % width of arrow head
    points = [ head-[0,width]         % left corner
        head-[width*1.5,0]         % right corner
        head+[0,width] ];      % vertex
    Screen('FillPoly', window, white, points);
    % create line for arrow
    lineCoords = [head', head'+[width*2,0]'];
    lineWidthPix = 4;
    Screen('DrawLines', window, lineCoords, lineWidthPix, white);
    
    % draw a right arrow
    % create a triangle
    head   = [screenXpixels*.75+35,  screenYpixels*.75+100]; % coordinates of head
    points = [ head-[0,width]         % left corner
        head+[width*1.5,0]         % right corner
        head+[0,width] ];      % vertex
    Screen('FillPoly', window, white, points);
    % create line for arrow
    lineCoords = [head', head'-[width*2,0]'];
    lineWidthPix = 4;
    Screen('DrawLines', window, lineCoords, lineWidthPix, white);
    
    vbl = Screen('Flip', window, vbl + numFrames*ifi);
    
    exitKeyboardLoop = 0;
    while exitKeyboardLoop == 0
        % Check the keyboard
        [keyIsDown, secs, keyCode] = KbCheck;
        if keyCode(escapeKey)
            % do something to save the trial as is, or double-check that you want
            % to stop
            exitKejyboardLoop = 1;
            sca;
            return
        elseif keyCode(leftKey) % subject indicated 'yes'
            resp = 1;
            exitKeyboardLoop = 1;
            vbl = Screen('Flip', window)
        elseif keyCode(rightKey) % subject indicated 'no'
            resp = 0;
            exitKeyboardLoop = 1;
        end
        clear keyCode keyIsDown secs
    end
    
    % Show yes/no gui
    %     resp = questdlg('Did you feel that?', ['Intensity = ', num2str(pulseAmp(i))]);
    doubleMatPos = find(stairOrderMatrix==i); % This is the linear index within the matrices that we're using to track responses, intensity, catch trials, step sizes, etc.
    
    % Now we assume that either the yesButton or the noButton was pressed
    % and use the 'yes' 'no' indicator targets in TDT to figure out which
    % one
    if isequal(resp,1) % yes was pressed, so step the intensity down
        % This pulls what was just answered/delivered into the final
        % response and intensity vectors
        response(i) = 1;
        staircase.response(doubleMatPos) = 1;
        intensity(i) = pulseAmp(i); % record the last pulse amp that was used
        
        if ~catchTrial(doubleMatPos)
            yesInaRow(whichStair) = yesInaRow(whichStair) + 1; % only increment the number of answers in a row when not a catch trial
            noInaRow(whichStair) = 0; % reset the no count
        end
        
        % set the next intensity value in the staircase matrix which tracks
        % both of the two staircases (one in each column)
        % but only want to do this if we're not in the last trial of that
        % staircase
        [row, col] = ind2sub(size(stairOrderMatrix), doubleMatPos);
        if row < numTrials/2
            if yesInaRow(whichStair) == yesNum
                staircase.intensity(doubleMatPos+1) = staircase.intensity(doubleMatPos) - stepSizeVec(doubleMatPos+1);
                yesInaRow(whichStair) = 0; % reset the count
            else % intensity doesn't change
                staircase.intensity(doubleMatPos+1) = staircase.intensity(doubleMatPos);
            end
        end
        
    elseif isequal (resp,0) % no was pressed, so step the intensity up
        response(i) = 0;
        staircase.response(doubleMatPos) = 0;
        intensity(i) = pulseAmp(i); % record the last pulse amp that was used
        
        if ~catchTrial(doubleMatPos)
            noInaRow(whichStair) = noInaRow(whichStair) + 1; % only increment the number of answers in a row when not a catch trial
            yesInaRow(whichStair) = 0; % reset the yes count
        end
        
        % set the next intensity value in the staircase matrix which tracks
        % both of the two staircases (one in each column)
        % but only want to do this if we're not in the last trial of that
        % staircase
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
        warning('Do not know which button was pressed')
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
    end
    
    pause(2)
    % Flip back to the fixation cross on subject's screen
    Screen('DrawLines', window, fixCrossCoords, lineWidthPix, white, [xCenter yCenter], 2);
    vbl = Screen('Flip', window);
    
end
% Clear the screen
sca;
%% Would want to end the connection to TDT here

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

% Also plot the staircases separately and without the catch trials
stair1.intensity = staircase.intensity(staircase.intensity(:,1)~=0,1)/1000;
stair1.response = staircase.response(staircase.intensity(:,1)~=0,1);
stair2.intensity = staircase.intensity(staircase.intensity(:,2)~=0,2)/1000;
stair2.response = staircase.response(staircase.intensity(:,2)~=0,2);

figure(1), title('Staircases Merged with Catch Trials')

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
% Compare staircase 1
figure(4)
subplot(2,1,1)
stairs(staircase.intensity(:,1));
hold on
x = 1:length(staircase.intensity(:,1));
id = staircase.response(:,1) == 1;
plot(x(id),staircase.intensity(id,1),'ko','MarkerFaceColor','g');
plot(x(~id),staircase.intensity(~id,1),'ko','MarkerFaceColor','r');
xlabel('Trial Number');
ylabel('Amplitude (mA)');
title('Staircase 1 with catch trials');

subplot(2,1,2)
stairs(stair1.intensity);
hold on
x = 1:length(stair1.intensity);
id = stair1.response == 1;
plot(x(id),stair1.intensity(id),'ko','MarkerFaceColor','g');
plot(x(~id),stair1.intensity(~id),'ko','MarkerFaceColor','r');
xlabel('Trial Number');
ylabel('Amplitude (mA)');
title('Without catch trials');


% Compare staircase 2
figure(5)
subplot(2,1,1)
stairs(staircase.intensity(:,2));
hold on
x = 1:length(staircase.intensity(:,2));
id = staircase.response(:,2) == 1;
plot(x(id),staircase.intensity(id,2),'ko','MarkerFaceColor','g');
plot(x(~id),staircase.intensity(~id,2),'ko','MarkerFaceColor','r');
xlabel('Trial Number');
ylabel('Amplitude (mA)');
title('Staircase 2 with catch trials');

subplot(2,1,2)
stairs(stair2.intensity);
hold on
x = 1:length(stair2.intensity);
id = stair2.response == 1;
plot(x(id),stair2.intensity(id),'ko','MarkerFaceColor','g');
plot(x(~id),stair2.intensity(~id),'ko','MarkerFaceColor','r');
xlabel('Trial Number');
ylabel('Amplitude (mA)');
title('Without catch trials');
