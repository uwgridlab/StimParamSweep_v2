function [filename] = getLatestFile(fileDir)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

a = dir(fileDir);

[~, index] = max([a.datenum]);


filename = a(index).name;

end

