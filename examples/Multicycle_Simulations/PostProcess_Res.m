load("LowRate.mat");
states1 = states;
load("NoReaction.mat");
states2 = states;
load('stepDates.mat');
load('schedule.mat') ;
[H2consCycles1, timeCycles1, H2consStep1, H2consTotal1, others1] = analyzeCycles_Res(states1, schedule, stepDates);
[H2consCycles2, timeCycles2, H2consStep2, H2consTotal2, others2] = analyzeCycles_Res(states2, schedule, stepDates);

% --- Case 1: With reaction ---
avgPresVals1 = others1.Average_Pres.Values;

% --- Case 2: No reaction ---
avgPresVals2 = others2.Average_Pres.Values;

% Reconstruct cycle boundaries from schedule
ctrl = schedule.step.control;
cycleEndIdx = find(ctrl(1:end-1) == 3 & ctrl(2:end) == 1); % 3->1 transition
cycleEndIdx = [cycleEndIdx; numel(avgPresVals1)];
cycleStartIdx = [1; cycleEndIdx(1:end-1)+1];
nCycles = numel(cycleStartIdx);

% Compute per-cycle averages
avgPresCycle1 = zeros(nCycles,1);
avgPresCycle2 = zeros(nCycles,1);
for c = 1:nCycles
    idxRange = cycleStartIdx(c):cycleEndIdx(c);
    avgPresCycle1(c) = mean(avgPresVals1(idxRange));
    avgPresCycle2(c) = mean(avgPresVals2(idxRange));
end

% --- Plot ---
figure;

% --- Left axis: pressure curves ---
yyaxis left
hold on;

% With Reaction
h1 = plot(1:nCycles, avgPresCycle1, '-ob','LineWidth',1.2, ...
    'MarkerSize',6,'MarkerFaceColor','r'); 

% No Reaction
h2 = plot(1:nCycles, avgPresCycle2, '--sk','LineWidth',1.2, ...
    'MarkerSize',6,'MarkerFaceColor','g'); 

ylabel('Average Reservoir Pressure [bar]');
xlabel('Cycle');
title('Cycle-Averaged Reservoir Pressure');
grid on;

% --- Right axis: difference bars ---
yyaxis right
diffCycle = avgPresCycle2 - avgPresCycle1;
h3 = bar(1:nCycles, diffCycle, 0.3, ...
    'FaceAlpha',0.4, 'FaceColor',[0.2 0.6 0.2]);
ylabel('Δ Pressure (No Reaction - Reaction) [bar]');

% Keep axes black
ax = gca;
ax.YColor = 'k';

% --- Legend with 3 entries only ---
legend([h1 h2 h3], {'With Reaction','No Reaction','Δ Pressure'}, 'Location','best');

%% Consumption
% Number of cycles
nCycles = numel(H2consCycles1);

% Preallocate
consMET  = zeros(nCycles,1);
consSRB  = zeros(nCycles,1);
consACE  = zeros(nCycles,1);
consTOT  = zeros(nCycles,1);

% Extract final cumulative consumption per cycle
for c = 1:nCycles
    consMET(c) = H2consCycles1{c}.MET(end);
    consSRB(c) = H2consCycles1{c}.SRB(end);
    consACE(c) = H2consCycles1{c}.ACE(end);
    consTOT(c) = H2consCycles1{c}.Total(end);
end

% Normalize by injected H2 (mol) and convert to percent
percMET = (consMET ./ injH2Cycles1) * 100;
percSRB = (consSRB ./ injH2Cycles1) * 100;
percACE = (consACE ./ injH2Cycles1) * 100;
percTOT = (consTOT ./ injH2Cycles1) * 100;

% --- Stacked bar data ---
Y = [percMET percSRB percACE]; % each column is a group

% --- Plot ---
figure; hold on;

% Stacked bars for MET, SRB, ACE
hBar = bar(1:nCycles, Y, 'stacked');
hBar(1).FaceColor = [0 0.4470 0.7410]; % MET (blue)
hBar(2).FaceColor = [0.8500 0.3250 0.0980]; % SRB (red)
hBar(3).FaceColor = [0.4660 0.6740 0.1880]; % ACE (green)

% Overlay total curve
plot(1:nCycles, percTOT, '-^k','LineWidth',1.5, ...
    'MarkerSize',6,'MarkerFaceColor','k');

xlabel('Cycle');
ylabel('H_2 Consumption [% of injected]');
title('Cycle-wise H_2 Consumption');
legend({'MET','SRB','ACE','Total'}, 'Location','best');
grid on;

