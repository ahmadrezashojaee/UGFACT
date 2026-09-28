function [schedule, W, stepDates] = makeAnnualStorageSchedule525(G, rock, cells, ...
            startYear, endYear, injType, injVal, prodType, prodVal, dt)
% Annual storage schedule with 5-2-5 scheme, exact calendar alignment
% Injection: Apr–Aug
% Shut-in:   Sep–Oct
% Production: Nov–Mar (next year)

import ad_core.*

years = startYear:endYear-1;
dtVec = [];
ctrlSeq = [];
stepDates = datetime.empty(0,1);

for y = years
    % Period lengths in days
    injDays  = daysInMonths(y, 4:8);
    shutDays = daysInMonths(y, 9:10);
    prodDays = daysInMonths(y, 11:12) + daysInMonths(y+1, 1:3);

    % Convert to seconds
    injSecs  = injDays*day; 
    shutSecs = shutDays*day;
    prodSecs = prodDays*day;

    % --- Injection steps ---
    nInj = floor(injSecs/dt);
    remInj = injSecs - nInj*dt;
    injDt = [repmat(dt,nInj,1); remInj*(remInj>0)];
    dtVec = [dtVec; injDt];
    ctrlSeq = [ctrlSeq; ones(numel(injDt),1)*1]; % Injection control

    % --- Shut-in steps ---
    nShut = floor(shutSecs/dt);
    remShut = shutSecs - nShut*dt;
    shutDt = [repmat(dt,nShut,1); remShut*(remShut>0)];
    dtVec = [dtVec; shutDt];
    ctrlSeq = [ctrlSeq; ones(numel(shutDt),1)*2]; % Shut-in control

    % --- Production steps ---
    nProd = floor(prodSecs/dt);
    remProd = prodSecs - nProd*dt;
    prodDt = [repmat(dt,nProd,1); remProd*(remProd>0)];
    dtVec = [dtVec; prodDt];
    ctrlSeq = [ctrlSeq; ones(numel(prodDt),1)*3]; % Production control

    % Save dates (cumulative sum of dt for this year’s cycle)
    if isempty(stepDates)
        lastDate = datetime(y,4,1);
    else
        lastDate = stepDates(end);
    end
    cumT = cumsum([injDt; shutDt; prodDt])/day;
    stepDates = [stepDates; lastDate + days(cumT)];
end

% Build schedule
schedule = simpleSchedule(dtVec);
schedule.step.control = ctrlSeq;

% ---- Wells ----
W = [];
% Injector
W = addWell(W, G, rock, cells, 'Type', injType, 'Val', injVal, ...
            'Name','Injector','comp_i',[0 1],'sign',+1,'radius',0.1,'dir', 'z');
W(end).components = [0, 1, 0, 0, 0, 0];

% Shut-in
W = addWell(W, G, rock, cells, 'Type', 'grat', 'Val', 0, ...
            'Name','Injector','comp_i',[0 1],'sign',+1,'radius',0.1,'dir', 'z');
W(end).components = [0, 1, 0, 0, 0, 0]; 
W(end).status = false;

% Producer
W = addWell(W, G, rock, cells, 'Type', prodType, 'Val', prodVal, ...
            'Name','Producer','comp_i',[0.5 0.5],'sign',-1,'radius',0.1,'dir', 'z');
W(end).components = [0, 1, 0, 0, 0, 0]; 
%W(end).lims.bhp = 150*barsa;

% Assign wells to controls
schedule.control = repmat(struct('W',[],'bc',[],'src',[]),3,1);
schedule.control(1).W = W(1); % Injection
schedule.control(2).W = W(2); % Shut-in
schedule.control(3).W = W(3); % Production

end

%% Helper
function nDays = daysInMonths(y, months)
    nDays = 0;
    for m = months
        nDays = nDays + eomday(y,m);
    end
end
