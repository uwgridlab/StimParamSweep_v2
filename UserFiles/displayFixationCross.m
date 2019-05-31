function [ vbl ] = displayFixationCross(window, white, xCenter, yCenter)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

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
end

