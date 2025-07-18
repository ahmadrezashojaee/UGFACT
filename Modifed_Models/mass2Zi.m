function states = mass2Zi(states,model,t)

mol_Water = states{t,1}.FlowProps.ComponentTotalMass{1,1}./0.018015268; % Water mole
mol_H2    = states{t,1}.FlowProps.ComponentTotalMass{2,1}./0.00201588;  % H2 mole
mol_CO2   = states{t,1}.FlowProps.ComponentTotalMass{3,1}./0.0440098;  % CO2 mole
mol_CH4   = states{t,1}.FlowProps.ComponentTotalMass{4,1}./0.0160428;  % CH4 mole
mol_H2S   = states{t,1}.FlowProps.ComponentTotalMass{5,1}./0.03408088;  % H2S mole
mol_N2    = states{t,1}.FlowProps.ComponentTotalMass{6,1}./0.02801348;  % N2 mole
Z_Water   = mol_Water./(mol_Water + mol_H2 + mol_CO2 + mol_CH4 + mol_H2S + mol_N2);
Z_H2      = mol_H2./(mol_Water + mol_H2 + mol_CO2 + mol_CH4 + mol_H2S + mol_N2);
Z_CO2     = mol_CO2./(mol_Water + mol_H2 + mol_CO2 + mol_CH4 + mol_H2S + mol_N2);
Z_CH4     = mol_CH4./(mol_Water + mol_H2 + mol_CO2 + mol_CH4 + mol_H2S + mol_N2);
Z_H2S     = mol_H2S./(mol_Water + mol_H2 + mol_CO2 + mol_CH4 + mol_H2S + mol_N2);
Z_N2      = mol_N2./(mol_Water + mol_H2 + mol_CO2 + mol_CH4 + mol_H2S + mol_N2);
eos = EquationOfStateModel([], model.EOSModel.CompositionalMixture, 'Peng-Robinson');
Pressure = states{t,1}.pressure;
Temp     = states{t,1}.T;
X1       = zeros(model.G.cells.num,1);
X2       = zeros(model.G.cells.num,1);
X3       = zeros(model.G.cells.num,1);
X4       = zeros(model.G.cells.num,1);
X5       = zeros(model.G.cells.num,1);
X6       = zeros(model.G.cells.num,1);

Y1       = zeros(model.G.cells.num,1);
Y2       = zeros(model.G.cells.num,1);
Y3       = zeros(model.G.cells.num,1);
Y4       = zeros(model.G.cells.num,1);
Y5       = zeros(model.G.cells.num,1);
Y6       = zeros(model.G.cells.num,1);

L        = zeros(model.G.cells.num,1);
Z_V      = zeros(model.G.cells.num,1);
Z_L      = zeros(model.G.cells.num,1);


Z1       = zeros(model.G.cells.num,1);
Z2       = zeros(model.G.cells.num,1);
Z3       = zeros(model.G.cells.num,1);
Z4       = zeros(model.G.cells.num,1);
Z5       = zeros(model.G.cells.num,1);
Z6       = zeros(model.G.cells.num,1);

for i=1:model.G.cells.num
    [l, x, y, Z_l, Z_v] = standaloneFlash(Pressure(i), Temp(i),...
                           [Z_Water(i), Z_H2(i), Z_CO2(i), Z_CH4(i), Z_H2S(i), Z_N2(i)], eos);
    X1(i,1)  = x(1);
    X2(i,1)  = x(2);
    X3(i,1)  = x(3);
    X4(i,1)  = x(4);
    X5(i,1)  = x(5);
    X6(i,1)  = x(6);
    Y1(i,1)  = y(1);
    Y2(i,1)  = y(2);
    Y3(i,1)  = y(3);
    Y4(i,1)  = y(4);
    Y5(i,1)  = y(5);
    Y6(i,1)  = y(6);
    L(i)     = l;
    Z1(i,1) = x(1)* l + y(1) * (1-l); %Zi
    Z2(i,1) = x(2)* l + y(2) * (1-l); %Zi
    Z3(i,1) = x(3)* l + y(3) * (1-l); %Zi
    Z4(i,1) = x(4)* l + y(4) * (1-l); %Zi
    Z5(i,1) = x(5)* l + y(5) * (1-l); %Zi
    Z6(i,1) = x(6)* l + y(6) * (1-l); %Zi
    Z_V(i)    = Z_v;
    Z_L(i)    = Z_l;
end
states{t,1}.x   = [X1,X2,X3,X4,X5,X6];
states{t,1}.y   = [Y1,Y2,Y3,Y4,Y5,Y6];
states{t,1}.components = [Z1,Z2,Z3,Z4,Z5,Z6];
states{t,1}.L   = L;
states{t,1}.Z_V = Z_V;
states{t,1}.Z_L = Z_L;
states{t,1}.K   = states{t,1}.y./states{t,1}.x;
end