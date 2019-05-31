function [ vbl ] = displayYesNoQuestion(window, white, screenXpixels, screenYpixels, answer)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

yellow  = [1 0.75 0];

if ~exist('answer', 'var')
    % Add text and yes/no buttons to the window
    Screen('TextSize', window, 56);
    DrawFormattedText(window, 'Did you feel the stimulation?', 'center', 'center', white);
    Screen('TextSize', window, 56);
    DrawFormattedText(window, 'Yes', screenXpixels*.25, screenYpixels*.75, white);
    DrawFormattedText(window, 'No', screenXpixels*.75, screenYpixels*.75, white);
    Screen('TextSize', window, 30);
    DrawFormattedText(window, 'press left arrow key', screenXpixels*.25-125, screenYpixels*.75+50, white);
    DrawFormattedText(window, 'press right arrow key', screenXpixels*.75-125, screenYpixels*.75+50, white);
    
else % then we're color-coding the yes or no response text based on the answer given
    if strcmp(answer, 'yes') % color the yes answer yellow
        Screen('TextSize', window, 56);
        DrawFormattedText(window, 'Did you feel the stimulation?', 'center', 'center', white);
        Screen('TextSize', window, 56);
        DrawFormattedText(window, 'Yes', screenXpixels*.25, screenYpixels*.75, yellow);
        DrawFormattedText(window, 'No', screenXpixels*.75, screenYpixels*.75, white);
        Screen('TextSize', window, 30);
        DrawFormattedText(window, 'press left arrow key', screenXpixels*.25-125, screenYpixels*.75+50, yellow);
        DrawFormattedText(window, 'press right arrow key', screenXpixels*.75-125, screenYpixels*.75+50, white);
    elseif strcmp(answer, 'no') % color the no answer yellow
        Screen('TextSize', window, 56);
        DrawFormattedText(window, 'Did you feel the stimulation?', 'center', 'center', white);
        Screen('TextSize', window, 56);
        DrawFormattedText(window, 'Yes', screenXpixels*.25, screenYpixels*.75, white);
        DrawFormattedText(window, 'No', screenXpixels*.75, screenYpixels*.75, yellow);
        Screen('TextSize', window, 30);
        DrawFormattedText(window, 'press left arrow key', screenXpixels*.25-125, screenYpixels*.75+50, white);
        DrawFormattedText(window, 'press right arrow key', screenXpixels*.75-125, screenYpixels*.75+50, yellow);
    else
        warning('answer variable in displayYesNoQuestion must be a string of yes or no')
    end
end

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

vbl = Screen('Flip', window);
end

