% Clear variables we don't want to save
clearvars i ind nextIntensity_afterCatch noInaRow stairIter startingPulseAmp 
clearvars whichStair yesInaRow DA ans col doubleMatPos
clearvars resp respInitialize row

blockName = getLatestFile(tank);
if experimentType==1 % amplitude experiment
    clearvars PTD_setting
    save([tank, '\', blockName, '_Amp_Matlab']);
elseif experimentType==0 % TD experiment
    save([tank, '\', blockName, '_TD_Matlab']);
else
    warning('This is not a recognized experiment type, variables were saved as _unknownType_Matlab')
    save([tank, '\', blockName, '_unknownType_Matlab']);
end





