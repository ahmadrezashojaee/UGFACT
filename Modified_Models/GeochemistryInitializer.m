function [model,state0] = GeochemistryInitializer(model,state0)

% Initializing Geochemistry variables

disp("****Start of Initializing Geochemistry****")
totTime = 1*day;
steps   = 1;
dt = totTime/steps;
dt = repmat(dt, steps, 1);
schedule = simpleSchedule(dt);
[~, states, ~] = simulateScheduleAD_Modified(state0, model, schedule);
%Solution
model.state0.Ca      = repmat(sum(states{steps,1}.Solution.Ca)/model.G.cells.num,model.G.cells.num,1);
model.state0.K       = repmat(sum(states{steps,1}.Solution.K)/model.G.cells.num,model.G.cells.num,1);
model.state0.Mg      = repmat(sum(states{steps,1}.Solution.Mg)/model.G.cells.num,model.G.cells.num,1);
model.state0.Na      = repmat(sum(states{steps,1}.Solution.Na)/model.G.cells.num,model.G.cells.num,1);
model.state0.C4      = repmat(sum(states{steps,1}.Solution.C4)/model.G.cells.num,model.G.cells.num,1);
model.state0.S6      = repmat(sum(states{steps,1}.Solution.S6)/model.G.cells.num,model.G.cells.num,1);
model.state0.Fe2     = repmat(sum(states{steps,1}.Solution.Fe2)/model.G.cells.num,model.G.cells.num,1);
model.state0.Fe3     = repmat(sum(states{steps,1}.Solution.Fe3)/model.G.cells.num,model.G.cells.num,1);
model.state0.Cl      = repmat(sum(states{steps,1}.Solution.Cl)/model.G.cells.num,model.G.cells.num,1);
model.state0.Si      = repmat(sum(states{steps,1}.Solution.Si)/model.G.cells.num,model.G.cells.num,1);
model.state0.S2      = repmat(sum(states{steps,1}.Solution.S2)/model.G.cells.num,model.G.cells.num,1);
model.state0.Acetate = repmat(sum(states{steps,1}.Solution.Acetate)/model.G.cells.num,model.G.cells.num,1);
model.state0.pH      = repmat(sum(states{steps,1}.Solution.pH)/model.G.cells.num,model.G.cells.num,1);

% Composition
state0.s = states{steps,1}.s;
state0.x = states{steps,1}.x;
state0.y = states{steps,1}.y;
state0.pressure = states{steps,1}.pressure;
state0.L = states{steps,1}.L;
state0.components = states{steps,1}.components;
state0.K = states{steps,1}.K;
state0.Z_V = states{steps,1}.Z_V;
state0.Z_L = states{steps,1}.Z_L;
model.state0.initGeoChem = false;
disp("****End of Initializing Geochemistry****")
end