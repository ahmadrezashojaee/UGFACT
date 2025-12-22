function states = initBiomassAndGeochem(model, states, t)
%INITBIOMASSANDGEOCHEM
% Initialise / carry over biomass, aqueous chemistry and mineralogy.

    G  = model.G;
    nc = G.cells.num;

    % ---------- Biomass initialisation ----------
    % Kinetic parameters
    M0_SRB = model.Kinetic.M0_SRB;
    M0_MET = model.Kinetic.M0_MET;
    M0_ACE = model.Kinetic.M0_ACE;
    M0_FRB = model.Kinetic.M0_FRB;

    if t == 1
        % If not already set, initialise biomass uniformly
        states{1,1}.Solution.SRB_Biomass = repmat(M0_SRB, nc, 1);
        states{1,1}.Solution.MET_Biomass = repmat(M0_MET, nc, 1);
        states{1,1}.Solution.ACE_Biomass = repmat(M0_ACE, nc, 1);
        states{1,1}.Solution.FRB_Biomass = repmat(M0_FRB, nc, 1);
    else
        % Carry over biomass from previous step
        states{t,1}.Solution.SRB_Biomass = states{t-1,1}.Solution.SRB_Biomass;
        states{t,1}.Solution.MET_Biomass = states{t-1,1}.Solution.MET_Biomass;
        states{t,1}.Solution.ACE_Biomass = states{t-1,1}.Solution.ACE_Biomass;
        states{t,1}.Solution.FRB_Biomass = states{t-1,1}.Solution.FRB_Biomass;
    end

    % ---------- Solution & Mineralogy ----------
    if t == 1
        % Decide where to take the initial aqueous chemistry from
        useState0 = isfield(model, 'state0') ...
                    && isfield(model.state0, 'initGeoChem') ...
                    && ~model.state0.initGeoChem;

        if useState0
            solSrc = model.state0;      % user-specified initial chemistry
        else
            solSrc = model.Solution;    % PHREEQC-equilibrated or default chemistry
        end

        % Always available chemistry as fallback
        solDefault = model.Solution;

        % 1) Unit
        if isfield(solSrc, 'Unit')
            states{1,1}.Solution.Unit = solSrc.Unit;
        elseif isfield(solDefault, 'Unit')
            states{1,1}.Solution.Unit = solDefault.Unit;
        else
            states{1,1}.Solution.Unit = 'mol/kgw';
        end

        % 2) Minerals always from model.Mineralogy
        states{1,1}.Mineralogy.Calcite     = model.Mineralogy.Calcite;
        states{1,1}.Mineralogy.Dolomite    = model.Mineralogy.Dolomite;
        states{1,1}.Mineralogy.Anhydrite   = model.Mineralogy.Anhydrite;
        states{1,1}.Mineralogy.Quartz      = model.Mineralogy.Quartz;
        states{1,1}.Mineralogy.Goethite    = model.Mineralogy.Goethite;
        states{1,1}.Mineralogy.Brucite     = model.Mineralogy.Brucite;
        states{1,1}.Mineralogy.Portlandite = model.Mineralogy.Portlandite;
        states{1,1}.Mineralogy.Pyrite      = model.Mineralogy.Pyrite;
        states{1,1}.Mineralogy.Gypsum      = model.Mineralogy.Gypsum;

        % 3) Aqueous species: try model.state0.*, otherwise fall back to model.Solution.*
        fieldsSol = {'Ca','Cl','Na','K','C4','S6','S2','Fe3','Mg','Fe2','Si','Acetate','pH'};

        for k = 1:numel(fieldsSol)
            fld = fieldsSol{k};

            if isfield(solSrc, fld)
                states{1,1}.Solution.(fld) = solSrc.(fld);
            elseif isfield(solDefault, fld)
                % Fallback if state0 does not have this field
                states{1,1}.Solution.(fld) = solDefault.(fld);
            else
                error('initBiomassAndGeochem:MissingField', ...
                      'Neither model.state0 nor model.Solution defines field "%s".', fld);
            end
        end

    else
        % t > 1 -> just copy Solution and Mineralogy from previous step
        prevSol = states{t-1,1}.Solution;
        prevMin = states{t-1,1}.Mineralogy;

        % Copy all relevant solution fields (add more if needed)
        fieldsSol = {'Unit','pH','Ca','Cl','Na','K','C4','S6','S2','Fe3','Mg','Fe2','Si','Acetate',...
                     'SRB_Biomass','MET_Biomass','ACE_Biomass','FRB_Biomass'};
        for k = 1:numel(fieldsSol)
            fld = fieldsSol{k};
            if isfield(prevSol, fld)
                states{t,1}.Solution.(fld) = prevSol.(fld);
            end
        end

        % Copy mineralogy
        fieldsMin = {'Calcite','Dolomite','Anhydrite','Quartz','Goethite','Brucite','Portlandite','Pyrite','Gypsum'};
        for k = 1:numel(fieldsMin)
            fld = fieldsMin{k};
            if isfield(prevMin, fld)
                states{t,1}.Mineralogy.(fld) = prevMin.(fld);
            end
        end
    end
end
