load('LowRate_wellSol.mat')
wellSol1 = wellSol;
load('NoReaction_wellSol.mat')
wellSol2 = wellSol;
load('stepDates.mat')
load('schedule.mat') 
[purityCycles1, recoveryCycles1, timeCycles1, H2SCycles1, avgH2S1, bhpCycles1, injH2Cycles1] = analyzeCycles_Well(wellSol1, schedule, stepDates);
[purityCycles2, recoveryCycles2, timeCycles2, H2SCycles2, avgH2S2, bhpCycles2, injH2Cycles2] = analyzeCycles_Well(wellSol2, schedule, stepDates);
nCycles = numel(purityCycles1);

%% 1. Recovery vs time (all cycles)
figure;
subplot(2,2,1)
hold on;
for c = 1:nCycles
    plot(timeCycles1{c}, 100*recoveryCycles1{c}, 'LineWidth', 1.2);
end
xlabel('Date'); ylabel('H_2 Recovery [%]');
title('H_2 Recovery vs Time');
%legend show;
grid on;

%% 2. Final Recovery vs cycle
finalRecovery1 = 100*cellfun(@(x) x(end), recoveryCycles1); % with reaction
finalRecovery2 = 100*cellfun(@(x) x(end), recoveryCycles2); % no reaction

subplot(2,2,2)

% --- Left axis: recovery curves ---
yyaxis left
plot(1:nCycles, finalRecovery1, '-ob','LineWidth',1.2, ...
    'MarkerSize',6,'MarkerFaceColor','b'); % blue line + markers
hold on
plot(1:nCycles, finalRecovery2, '-sr','LineWidth',1.2, ...
    'MarkerSize',6,'MarkerFaceColor','r'); % red line + markers
ylabel('Final Recovery [%]');
xlabel('Cycle');
title('Final H_2 Recovery per Cycle');
grid on;

ax = gca;
ax.YColor = 'k';  % keep axis black

% --- Right axis: difference bar chart ---
yyaxis right
diffRecovery = finalRecovery2 - finalRecovery1;
bar(1:nCycles, diffRecovery, 0.3, 'FaceAlpha',0.4, 'FaceColor',[0.2 0.6 0.2]);
ylabel('Δ Recovery (NoReaction - Reaction) [%]');

ax.YColor = 'k';  % keep axis black

legend({'With Microbial Reaction','No Microbial Reaction','Δ Purity'}, 'Location','best');


%% 3. Purity vs time (all cycles)
%figure;
subplot(2,2,3)
hold on;
for c = 1:nCycles
    plot(timeCycles1{c}, 100*purityCycles1{c}, 'LineWidth', 1.2);
end
xlabel('Date'); ylabel('H_2 Purity [%]');
title('H_2 Purity vs Time');
%legend show;
grid on;

%% 4. Average Purity vs cycle
avgPurity1 = 100*cellfun(@mean, purityCycles1); % with reaction
avgPurity2 = 100*cellfun(@mean, purityCycles2); % no reaction

subplot(2,2,4)

% --- Left axis: average purity curves ---
yyaxis left
plot(1:nCycles, avgPurity1, '-ob','LineWidth',1.2, ...
    'MarkerSize',6,'MarkerFaceColor','b'); % blue line + markers
hold on
plot(1:nCycles, avgPurity2, '-sr','LineWidth',1.2, ...
    'MarkerSize',6,'MarkerFaceColor','r'); % red line + markers

ylabel('Average Purity [%]');
xlabel('Cycle');
title('Average H_2 Purity per Cycle');
grid on;

% Force axis color to black
ax = gca;
ax.YColor = 'k'; 

% --- Right axis: difference bar chart ---
yyaxis right
diffPurity = avgPurity2 - avgPurity1;
bar(1:nCycles, diffPurity, 0.3, 'FaceAlpha',0.4, 'FaceColor',[0.2 0.6 0.2]);
ylabel('Δ Purity (No Reaction - Reaction) [%]');

% Force axis color to black again
ax.YColor = 'k';

legend({'With Microbial Reaction','No Microbial Reaction','Δ Purity'}, 'Location','best');


%% 5. H2S concentration vs time (all cycles)
figure; 
subplot(1,2,1)
hold on;
for c = 1:nCycles
    plot(timeCycles1{c}, H2SCycles1{c}, 'LineWidth', 1.2);
end
xlabel('Date'); ylabel('H_2S Concentration [ppm]');
title('H_2S Concentration vs Time');
%legend show; 
grid on;

%% 6. Average H2S per cycle
subplot(1,2,2); hold on;
plot(1:nCycles, avgH2S1, '-b','LineWidth',1.2);
plot(1:nCycles, avgH2S1, 'or','MarkerSize',6,'MarkerFaceColor','r');
xlabel('Cycle'); ylabel('Average H_2S [ppm]');
title('Average H_2S per Cycle');
grid on;

%% BHP
% Recompute injection start/end outside the function
ctrl = schedule.step.control;
dCtrl = diff([0; ctrl==1]);
injStart = find(dCtrl==1);
dCtrl = diff([ctrl==3; 0]);
prodEnd = find(dCtrl==-1);

nCycles = min([numel(injStart), numel(prodEnd), numel(bhpCycles1), numel(bhpCycles2)]);

% Make the figure wider
figure('Units','normalized','Position',[0.05 0.2 0.85 0.5]); 
hold on;

h1 = [];   % will hold first Reaction curve
h2 = [];   % will hold first No-Reaction curve

% Plot all cycles
for c = 1:nCycles
    range = injStart(c):prodEnd(c);

    % With reaction (blue)
    hR = plot(stepDates(range), bhpCycles1{c}, 'b-', 'LineWidth', 1.2);
    if isempty(h1), h1 = hR; end   % store only the FIRST blue line

    % No reaction (black)
    hN = plot(stepDates(range), bhpCycles2{c}, 'k--', 'LineWidth', 1.2);
    if isempty(h2), h2 = hN; end   % store only the FIRST black line
end

xlabel('Date');
ylabel('BHP [bar]');
title('Bottomhole Pressure vs Time (All Cycles)');
grid on;

% ---- Two-entry legend ----
legend([h1 h2], {'Moderate-Rate', 'No-Reaction'});
