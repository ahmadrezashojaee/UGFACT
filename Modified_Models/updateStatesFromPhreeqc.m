function [states, reactionRates] = updateStatesFromPhreeqc(model, states, t, OUTPUT, phaseData, mineralData)
%UPDATESTATESFROMPHREEQC
% Use PHREEQC OUTPUT to update mineralogy, solution chemistry, pressure and
% component masses at time step t. Only cells flagged in phaseData.activeCells
% are updated by PHREEQC; others keep previous Solution/Mineralogy and FlowProps.

    G     = model.G;
    rock  = model.rock;
    nc    = G.cells.num;

    % ---------- Unpack phaseData ----------
    Sw               = phaseData.Sw;
    Sg               = phaseData.Sg;
    X                = phaseData.X;
    Y                = phaseData.Y;
    rho_molar_Gas    = phaseData.rho_molar_Gas;

    mass_H2O_Profile = phaseData.mass_H2O_Profile;
    mol_H2O_Profile  = phaseData.mol_H2O_Profile; %#ok<NASGU>

    MW_H2O           = phaseData.MW_H2O;
    MW_H2            = phaseData.MW_H2;
    MW_CO2           = phaseData.MW_CO2;
    MW_CH4           = phaseData.MW_CH4;
    MW_H2S           = phaseData.MW_H2S;
    MW_N2            = phaseData.MW_N2;

    pH               = phaseData.pH;
    K                = phaseData.K;
    Na               = phaseData.Na;
    Mg               = phaseData.Mg;
    Ca               = phaseData.Ca;
    Cl               = phaseData.Cl;
    C4               = phaseData.C4;
    S6               = phaseData.S6;
    S2               = phaseData.S2;
    Fe3              = phaseData.Fe3;
    Fe2              = phaseData.Fe2;
    Si               = phaseData.Si;
    Acetate          = phaseData.Acetate;

    SRB_Biomass      = phaseData.SRB_Biomass;
    MET_Biomass      = phaseData.MET_Biomass;
    ACE_Biomass      = phaseData.ACE_Biomass;
    FRB_Biomass      = phaseData.FRB_Biomass;

    pressure_atm     = phaseData.pressure_atm; %#ok<NASGU>  % might be useful if you want
    temp             = phaseData.temp;        % [K]
    Rgas             = phaseData.R;           % atm·L/(mol·K)

    idxKeepOldP      = phaseData.idxKeepOldPressure;
    activeCells      = phaseData.activeCells;
    steps            = phaseData.steps;

    % ---------- Unpack mineralData (initial moles per cell) ----------
    Calcite_i        = mineralData.Calcite_i;
    Dolomite_i       = mineralData.Dolomite_i;
    Quartz_i         = mineralData.Quartz_i;
    Anhydrite_i      = mineralData.Anhydrite_i;
    Goethite_i       = mineralData.Goethite_i;
    Brucite_i        = mineralData.Brucite_i;
    Portlandite_i    = mineralData.Portlandite_i;
    Pyrite_i         = mineralData.Pyrite_i;
    Gypsum_i         = mineralData.Gypsum_i;

    % ---------- Reaction rates (PHREEQC outputs) ----------
    MET_Rate   = zeros(nc,1);
    ACE_Rate   = zeros(nc,1);
    SRB_Rate   = zeros(nc,1);
    D_Calcite  = zeros(nc,1);
    D_Dolomite = zeros(nc,1);

    % ---------- Start from current state: only Solution & Mineralogy ----------
    % For t>1, we want Solution and Mineralogy to default to values from t-1.
    % FlowProps and other fields should already be in states{t,1} from MRST.

    if t > 1
        prevSol = states{t-1,1}.Solution;
        prevMin = states{t-1,1}.Mineralogy;

        solFields = {'pH','Ca','Cl','Na','Mg','K','C4','S6','S2','Fe3','Fe2','Si','Acetate',...
                     'SRB_Biomass','MET_Biomass','ACE_Biomass','FRB_Biomass'};
        minFields = {'Calcite','Dolomite','Anhydrite','Quartz','Goethite','Brucite','Portlandite','Pyrite','Gypsum'};

        for f = 1:numel(solFields)
            fld = solFields{f};
            states{t,1}.Solution.(fld) = prevSol.(fld);
        end
        for f = 1:numel(minFields)
            fld = minFields{f};
            states{t,1}.Mineralogy.(fld) = prevMin.(fld);
        end
    end

    % Now read these back into local arrays (for all cells)
    Calcite     = states{t,1}.Mineralogy.Calcite;
    Dolomite    = states{t,1}.Mineralogy.Dolomite;
    Anhydrite   = states{t,1}.Mineralogy.Anhydrite;
    Quartz      = states{t,1}.Mineralogy.Quartz;
    Goethite    = states{t,1}.Mineralogy.Goethite;
    Brucite     = states{t,1}.Mineralogy.Brucite;
    Portlandite = states{t,1}.Mineralogy.Portlandite;
    Pyrite      = states{t,1}.Mineralogy.Pyrite;
    Gypsum      = states{t,1}.Mineralogy.Gypsum;

    % Water in Solution: if not present yet, initialise with current mass_H2O_Profile
    if isfield(states{t,1}.Solution, 'Water')
        Water = states{t,1}.Solution.Water;
    else
        Water = mass_H2O_Profile;
    end

    % Total moles per gas component (only needed for active cells)
    H2  = zeros(nc,1);
    CO2 = zeros(nc,1);
    CH4 = zeros(nc,1);
    H2S = zeros(nc,1);
    N2  = zeros(nc,1);

    % Start pressure from MRST (Pa); active cells may be updated
    pressure = states{t,1}.pressure;

    % ---------- Main loop over cells ----------
    for i = 1:nc
        if activeCells(i) && ~isempty(OUTPUT{i,1}) && iscell(OUTPUT{i,1})
            output = OUTPUT{i,1};

            % ---- Integrate rates over time (columns 64,65,66) ----
            a = 0; b = 0; c = 0;
            for j = 2:steps+1
                if j == 2
                    a = output{j,64} .* output{j,1};
                    b = output{j,65} .* output{j,1};
                    c = output{j,66} .* output{j,1};
                else
                    dt_local = output{j,1} - output{j-1,1};
                    a = a + output{j,64} .* dt_local;
                    b = b + output{j,65} .* dt_local;
                    c = c + output{j,66} .* dt_local;
                end
            end

            MET_Rate(i) = a ./ output{steps+1,1} * 1000 * 86400; % mmol/day/kgw
            ACE_Rate(i) = b ./ output{steps+1,1} * 1000 * 86400; % mmol/day/kgw
            SRB_Rate(i) = c ./ output{steps+1,1} * 1000 * 86400; % mmol/day/kgw

            % ---- Calcite / Dolomite change (columns 24,30) ----
            d_calcite  = 0;
            d_dolomite = 0;
            for jj = 2:steps+1
                d_calcite  = d_calcite  + output{jj,24} * mass_H2O_Profile(i);
                d_dolomite = d_dolomite + output{jj,30} * mass_H2O_Profile(i);
            end
            D_Calcite(i)  = d_calcite;
            D_Dolomite(i) = d_dolomite;

            % ---- Updated minerals ----
            Calcite(i)     = Calcite_i(i)   + d_calcite;
            Dolomite(i)    = Dolomite_i(i)  + d_dolomite;
            Anhydrite(i)   = output{end,25} * mass_H2O_Profile(i);
            Quartz(i)      = output{end,39} * mass_H2O_Profile(i);
            Goethite(i)    = output{end,31} * mass_H2O_Profile(i);
            Brucite(i)     = output{end,35} * mass_H2O_Profile(i);
            Portlandite(i) = output{end,37} * mass_H2O_Profile(i);
            Pyrite(i)      = output{end,33} * mass_H2O_Profile(i);
            Gypsum(i)      = output{end,27} * mass_H2O_Profile(i);

            % ---- Updated aqueous chemistry (end row indices) ----
            pH(i)       = output{end,3};
            Ca(i)       = output{end,8};
            Cl(i)       = output{end,14};
            Na(i)       = output{end,13};
            Mg(i)       = output{end,9};
            K(i)        = output{end,12};
            % C4: new amendment vs dissolved CO2
            C4(i)       = output{end,6} - output{end,20};
            S6(i)       = output{end,7};
            S2(i)       = output{end,16};
            Fe3(i)      = output{end,11};
            Fe2(i)      = output{end,10};
            Si(i)       = output{end,17};
            Acetate(i)  = output{end,15};

            % ---- Pressure from PHREEQC and EOS correction ----
            % first: PHREEQC gas pressure (column 48) -> Pa
            pressure(i) = output{end,48} * 101325;

            % ---- Water mass (column 5) ----
            Water(i)    = output{end,5} * mass_H2O_Profile(i);

            % ---- Gas EOS correction (following your original logic) ----
            xH2_MRST  = states{t,1}.x(i,2) / states{t,1}.x(i,1) / MW_H2O;
            xCH4_MRST = states{t,1}.x(i,4) / states{t,1}.x(i,1) / MW_H2O;
            xCO2_MRST = states{t,1}.x(i,3) / states{t,1}.x(i,1) / MW_H2O;

            CH4moles = (output{end,21} - xCH4_MRST) * Water(i);
            H2moles  = (output{end,18} - xH2_MRST)  * Water(i);
            CO2moles = (output{end,71} - xCO2_MRST) * Water(i);

            Z_factor = (output{end,72} + states{t,1}.Z_V(i)) / 2;

            % output{end,49} total gas moles/kgw, output{end,50} volume/kgw
            pressure(i) = 101325 * Z_factor * ...
                (CO2moles + CH4moles + H2moles + output{end,49} * mass_H2O_Profile(i)) .* ...
                Rgas .* temp(i) ./ output{end,50} ./ mass_H2O_Profile(i);

            % keep old pressure if flagged
            if idxKeepOldP(i)
                pressure(i) = states{t,1}.pressure(i);
            end

            % ---- Biomass (columns 56,58,60,62) ----
            SRB_Biomass(i) = output{end,56} * mass_H2O_Profile(i);
            MET_Biomass(i) = output{end,58} * mass_H2O_Profile(i);
            ACE_Biomass(i) = output{end,60} * mass_H2O_Profile(i);
            FRB_Biomass(i) = output{end,62} * mass_H2O_Profile(i);

            % update water mass profile factor
            mass_H2O_Profile(i) = output{end,5} * mass_H2O_Profile(i);

            % ---- Total moles (aqueous + gas) from PHREEQC (columns 51–55,18–21,19) ----
            H2(i,1)  = output{end,51} * mass_H2O_Profile(i) + output{end,18} * Water(i);
            CO2(i,1) = output{end,52} * mass_H2O_Profile(i) + output{end,20} * Water(i);
            CH4(i,1) = output{end,53} * mass_H2O_Profile(i) + output{end,21} * Water(i);
            H2S(i,1) = output{end,54} * mass_H2O_Profile(i);
            N2(i,1)  = output{end,55} * mass_H2O_Profile(i) + output{end,19} * Water(i);

        else
            % inactive cell -> keep Solution, Mineralogy, FlowProps as they are
            continue;
        end
    end

    % ---------- Cumulative mineral changes ----------
    if t == 1
        states{t,1}.Mineralogy.Delta_Calcite  = D_Calcite;
        states{t,1}.Mineralogy.Delta_Dolomite = D_Dolomite;
    else
        states{t,1}.Mineralogy.Delta_Calcite  = D_Calcite  + states{t-1,1}.Mineralogy.Delta_Calcite;
        states{t,1}.Mineralogy.Delta_Dolomite = D_Dolomite + states{t-1,1}.Mineralogy.Delta_Dolomite;
    end

    % ---------- Write updated minerals ----------
    states{t,1}.Mineralogy.Calcite     = Calcite;
    states{t,1}.Mineralogy.Dolomite    = Dolomite;
    states{t,1}.Mineralogy.Anhydrite   = Anhydrite;
    states{t,1}.Mineralogy.Quartz      = Quartz;
    states{t,1}.Mineralogy.Goethite    = Goethite;
    states{t,1}.Mineralogy.Brucite     = Brucite;
    states{t,1}.Mineralogy.Portlandite = Portlandite;
    states{t,1}.Mineralogy.Pyrite      = Pyrite;
    states{t,1}.Mineralogy.Gypsum      = Gypsum;

    % ---------- Updated solution chemistry & biomass ----------
    states{t,1}.Solution.pH        = pH;
    states{t,1}.Solution.Ca        = Ca;
    states{t,1}.Solution.Cl        = Cl;
    states{t,1}.Solution.Na        = Na;
    states{t,1}.Solution.Mg        = Mg;
    states{t,1}.Solution.K         = K;
    states{t,1}.Solution.C4        = C4;
    states{t,1}.Solution.S6        = S6;
    states{t,1}.Solution.S2        = S2;
    states{t,1}.Solution.Fe3       = Fe3;
    states{t,1}.Solution.Fe2       = Fe2;
    states{t,1}.Solution.Si        = Si;
    states{t,1}.Solution.Acetate   = Acetate;
    states{t,1}.Solution.Water     = Water;

    states{t,1}.Solution.MET_Rate  = MET_Rate;
    states{t,1}.Solution.ACE_Rate  = ACE_Rate;
    states{t,1}.Solution.SRB_Rate  = SRB_Rate;

    states{t,1}.pressure           = pressure;

    states{t,1}.Solution.SRB_Biomass = SRB_Biomass;
    states{t,1}.Solution.MET_Biomass = MET_Biomass;
    states{t,1}.Solution.ACE_Biomass = ACE_Biomass;
    states{t,1}.Solution.FRB_Biomass = FRB_Biomass;

    % ---------- Recompute ComponentTotalMass only for active cells ----------
    for i = 1:nc
        if activeCells(i)
            % Water: liquid + water in gas phase
            states{t,1}.FlowProps.ComponentTotalMass{1,1}(i) = ...
                Water(i) + ...
                states{t,1}.s(i,2) * G.cells.volumes(i) * rock.poro(i) * ...
                rho_molar_Gas(i) * Y(i,1) * MW_H2O;

            states{t,1}.FlowProps.ComponentTotalMass{2,1}(i) = H2(i)  * MW_H2;
            states{t,1}.FlowProps.ComponentTotalMass{3,1}(i) = CO2(i) * MW_CO2;
            states{t,1}.FlowProps.ComponentTotalMass{4,1}(i) = CH4(i) * MW_CH4;
            states{t,1}.FlowProps.ComponentTotalMass{5,1}(i) = H2S(i) * MW_H2S;
            states{t,1}.FlowProps.ComponentTotalMass{6,1}(i) = N2(i)  * MW_N2;
        end
    end

    % ---------- Pack reaction rates ----------
    reactionRates = struct( ...
        'MET_Rate',   MET_Rate, ...
        'ACE_Rate',   ACE_Rate, ...
        'SRB_Rate',   SRB_Rate, ...
        'D_Calcite',  D_Calcite, ...
        'D_Dolomite', D_Dolomite ...
    );
end
