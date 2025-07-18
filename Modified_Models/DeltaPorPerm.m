function [model,states] = DeltaPorPerm(model,states,t)

if t == 1
    model.initial_Phi  = model.rock.poro;
    model.initial_Perm = model.rock.perm;
end

V_b     = model.G.cells.volumes;

Delta_V_Minerals = -1e-5 .* ( ...
    states{t,1}.Mineralogy.Delta_Calcite    .* 3.693 + ...
    states{t,1}.Mineralogy.Delta_Dolomite   .* 6.440 + ...
    states{t,1}.Mineralogy.Delta_Quartz     .* 2.269 + ...
    states{t,1}.Mineralogy.Delta_Anhydrite  .* 4.627 + ...
    states{t,1}.Mineralogy.Delta_Gypsum     .* 7.469 + ...
    states{t,1}.Mineralogy.Delta_Pyrite     .* 2.390 + ...
    states{t,1}.Mineralogy.Delta_Brucite    .* 2.463 + ...
    states{t,1}.Mineralogy.Delta_Portlandite.* 3.306 + ...
    states{t,1}.Mineralogy.Delta_Goethite   .* 2.088);
rho_biomass = 10; %Kg/m^3
Delta_V_Biomass = -states{t,1}.Solution.SRB_Biomass.*model.Kinetic.MW_SRB./1000./rho_biomass;
model.rock.poro = model.initial_Phi + Delta_V_Minerals./V_b;% + Delta_V_Biomass./V_b;
model.rock.perm = model.initial_Perm.*(model.rock.poro./model.initial_Phi).^3.*((1-model.initial_Phi)./(1-model.rock.poro)).^2;
states{t,1}.Mineralogy.Porosity = model.rock.poro;
states{t,1}.Mineralogy.Perm     = convertTo(model.rock.perm,milli*darcy);
states{t,1}.Mineralogy.DeltaPorMinerals = Delta_V_Minerals./V_b;
states{t,1}.Mineralogy.DeltaPorBiomass  = Delta_V_Biomass./V_b;
model.operators = setupOperatorsTPFA(model.G,model.rock);

end