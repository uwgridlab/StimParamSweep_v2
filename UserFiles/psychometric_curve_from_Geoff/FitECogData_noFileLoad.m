% FitECogData_noFileLoad.m

% Initial parameter values for the Weibull function
% p.g = 0; % guess rate (expected proportion correct at zero intensity)
p.b = 5; % slope 
p.t = 1; % threshold yeilding 79% performance

err = fitPsychometricFunction(p,results,'Weibull');

pBest = fit('fitPsychometricFunction',p,{'b','t'},results,'Weibull');

figure(1)
clf
stairs(results.intensity);
hold on

x = 1:length(results.intensity);

id = results.response ==1;
plot(x(id),results.intensity(id),'ko','MarkerFaceColor','g');
plot(x(~id),results.intensity(~id),'ko','MarkerFaceColor','r');
xlabel('Trial Number');
ylabel('Amplitude (mA)');
title('Staircase');

figure(2)
clf
plotPsycho(results,'Amplutide (mA)',pBest,'Weibull');
