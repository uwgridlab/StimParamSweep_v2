odd = 1:2:30;
even=2:2:30;
stair1.intensity = PAmp.data(odd);
stair2.intensity = PAmp.data(even);
figure, plot(Butt.data(:,2))

%%
resp = zeros(length(PAmp.data), 2);
yes = find(diff(Butt.data(:,2))==1);
resp(1:length(yes),1) = 1;
resp(1:length(yes),2) = yes;

no = find(diff(Butt.data(:,3))==1);
resp(length(yes)+1:end,1) = 0;
resp(length(yes)+1:end,2) = no;

resp = sortrows(resp,2);

%%
stair1.response = resp(odd,1);
stair2.response = resp(even,1);

staircase.intensity(:,1) = stair1.intensity;
staircase.intensity(:,2) = stair2.intensity;
staircase.response(:,1) = stair1.response;
staircase.response(:,2) = stair2.response;

%%
stair1.intensity = staircase.intensity(staircase.intensity(:,1)~=0,1)/1000;
stair1.response = staircase.response(staircase.intensity(:,1)~=0,1);
stair2.intensity = staircase.intensity(staircase.intensity(:,2)~=0,2)/1000;
stair2.response = staircase.response(staircase.intensity(:,2)~=0,2);


%%
results.intensity = PAmp.data/1000;
results.response = resp(:,1);

%% plot stuff
guessRate = sum(results.response(results.intensity==0))/sum(results.intensity==0);
p.g = guessRate;
% FitECogData_noFileLoad

figure
stairs(stair1.intensity);
hold on
x = 1:length(stair1.intensity);
id = stair1.response == 1;
plot(x(id),stair1.intensity(id),'ko','MarkerFaceColor','g');
plot(x(~id),stair1.intensity(~id),'ko','MarkerFaceColor','r');

stairs(stair2.intensity);
x = 1:length(stair2.intensity);
id = stair2.response == 1;
plot(x(id),stair2.intensity(id),'ko','MarkerFaceColor','g');
plot(x(~id),stair2.intensity(~id),'ko','MarkerFaceColor','r');

xlabel('Trial Number');
ylabel('Amplitude (mA)');
title('Two Staircases');

% Compare staircase 1
figure
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
figure
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



