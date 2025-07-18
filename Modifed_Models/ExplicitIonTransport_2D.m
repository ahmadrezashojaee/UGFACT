function states = ExplicitIonTransport_2D_New(model,states,t,dt)

G = model.G;
Nx = G.cartDims(1);
Nz = G.cartDims(3);
Water_Density = states{t,1}.PVTProps.Density{1,1};
input_Water1    = zeros(G.cells.num,1);
input_Water2    = zeros(G.cells.num,1);
input_Water5    = zeros(G.cells.num,1);
input_Water6    = zeros(G.cells.num,1);
for i=1:G.cells.num %calculation for each cell in 2D system
    CurrentWater = states{t,1}.FlowProps.ComponentTotalMass{1,1};
    CurrentWater = CurrentWater(i);
    rhoW = Water_Density(i);
    [cf,~] = gridCellFaces(G,i);
    rem = mod(i,Nx);
    q   = floor(i/Nx);
    if (0 < rem) && (rem <= Nx-1)
        nx = rem;
    else
        nx = Nx;
    end
    if mod(nx,Nx) == 0
        nz = i/Nx;
    else
        nz = q + 1;
    end
    flux = states{t,1}.flux(cf,1);
    if (nx == 1) && (nz == 1) % Top Left Cell (1)
        if flux(2) < 0
            input_Water2(i,1) = flux(2)*rhoW*dt;
        end
        if flux(6) < 0
            input_Water6(i,1) = flux(6)*rhoW*dt;
        end
        states{t,1}.Solution.Ca(i) = max((input_Water2(i,1)*states{t,1}.Solution.Ca(i+1) + CurrentWater.*states{t,1}.Solution.Ca(i) + input_Water6(i,1)*states{t,1}.Solution.Ca(i+Nx))./(input_Water2(i,1)+CurrentWater+input_Water6(i,1)),0);
        states{t,1}.Solution.Na(i) = max((input_Water2(i,1)*states{t,1}.Solution.Na(i+1) + CurrentWater.*states{t,1}.Solution.Na(i) + input_Water6(i,1)*states{t,1}.Solution.Na(i+Nx))./(input_Water2(i,1)+CurrentWater+input_Water6(i,1)),0);
        states{t,1}.Solution.Mg(i) = max((input_Water2(i,1)*states{t,1}.Solution.Mg(i+1) + CurrentWater.*states{t,1}.Solution.Mg(i) + input_Water6(i,1)*states{t,1}.Solution.Mg(i+Nx))./(input_Water2(i,1)+CurrentWater+input_Water6(i,1)),0);
        states{t,1}.Solution.Cl(i) = max((input_Water2(i,1)*states{t,1}.Solution.Cl(i+1) + CurrentWater.*states{t,1}.Solution.Cl(i) + input_Water6(i,1)*states{t,1}.Solution.Cl(i+Nx))./(input_Water2(i,1)+CurrentWater+input_Water6(i,1)),0);
        states{t,1}.Solution.K(i)  = max((input_Water2(i,1)*states{t,1}.Solution.K(i+1)  + CurrentWater.*states{t,1}.Solution.K(i)  + input_Water6(i,1)*states{t,1}.Solution.K(i+Nx))./(input_Water2(i,1)+CurrentWater+input_Water6(i,1)),0);
        states{t,1}.Solution.C4(i) = max((input_Water2(i,1)*states{t,1}.Solution.C4(i+1) + CurrentWater.*states{t,1}.Solution.C4(i) + input_Water6(i,1)*states{t,1}.Solution.C4(i+Nx))./(input_Water2(i,1)+CurrentWater+input_Water6(i,1)),0);
        states{t,1}.Solution.S6(i) = max((input_Water2(i,1)*states{t,1}.Solution.S6(i+1) + CurrentWater.*states{t,1}.Solution.S6(i) + input_Water6(i,1)*states{t,1}.Solution.S6(i+Nx))./(input_Water2(i,1)+CurrentWater+input_Water6(i,1)),0);
        states{t,1}.Solution.S2(i) = max((input_Water2(i,1)*states{t,1}.Solution.S2(i+1) + CurrentWater.*states{t,1}.Solution.S2(i) + input_Water6(i,1)*states{t,1}.Solution.S2(i+Nx))./(input_Water2(i,1)+CurrentWater+input_Water6(i,1)),0);
        states{t,1}.Solution.Si(i) = max((input_Water2(i,1)*states{t,1}.Solution.Si(i+1) + CurrentWater.*states{t,1}.Solution.Si(i) + input_Water6(i,1)*states{t,1}.Solution.Si(i+Nx))./(input_Water2(i,1)+CurrentWater+input_Water6(i,1)),0);
        states{t,1}.Solution.Fe2(i) = max((input_Water2(i,1)*states{t,1}.Solution.Fe2(i+1) + CurrentWater.*states{t,1}.Solution.Fe2(i) + input_Water6(i,1)*states{t,1}.Solution.Fe2(i+Nx))./(input_Water2(i,1)+CurrentWater+input_Water6(i,1)),0);
        states{t,1}.Solution.Fe3(i) = max((input_Water2(i,1)*states{t,1}.Solution.Fe3(i+1) + CurrentWater.*states{t,1}.Solution.Fe3(i) + input_Water6(i,1)*states{t,1}.Solution.Fe3(i+Nx))./(input_Water2(i,1)+CurrentWater+input_Water6(i,1)),0);
        states{t,1}.Solution.Acetate(i) = max((input_Water2(i,1)*states{t,1}.Solution.Acetate(i+1) + CurrentWater.*states{t,1}.Solution.Acetate(i) + input_Water6(i,1)*states{t,1}.Solution.Acetate(i+Nx))./(input_Water2(i,1)+CurrentWater+input_Water6(i,1)),0);
    elseif (nx == Nx) && (nz == 1) % Top Right Cell (2)
        if flux(1) > 0
            input_Water1(i,1) = flux(1)*rhoW*dt;
        end
        if flux(6) < 0
            input_Water6(i,1) = flux(6)*rhoW*dt;
        end
        states{t,1}.Solution.Ca(i) = max((input_Water1(i,1)*states{t,1}.Solution.Ca(i-1) + CurrentWater.*states{t,1}.Solution.Ca(i) + input_Water6(i,1)*states{t,1}.Solution.Ca(i+Nx))./(input_Water1(i,1)+CurrentWater+input_Water6(i,1)),0);
        states{t,1}.Solution.Na(i) = max((input_Water1(i,1)*states{t,1}.Solution.Na(i-1) + CurrentWater.*states{t,1}.Solution.Na(i) + input_Water6(i,1)*states{t,1}.Solution.Na(i+Nx))./(input_Water1(i,1)+CurrentWater+input_Water6(i,1)),0);
        states{t,1}.Solution.Mg(i) = max((input_Water1(i,1)*states{t,1}.Solution.Mg(i-1) + CurrentWater.*states{t,1}.Solution.Mg(i) + input_Water6(i,1)*states{t,1}.Solution.Mg(i+Nx))./(input_Water1(i,1)+CurrentWater+input_Water6(i,1)),0);
        states{t,1}.Solution.Cl(i) = max((input_Water1(i,1)*states{t,1}.Solution.Cl(i-1) + CurrentWater.*states{t,1}.Solution.Cl(i) + input_Water6(i,1)*states{t,1}.Solution.Cl(i+Nx))./(input_Water1(i,1)+CurrentWater+input_Water6(i,1)),0);
        states{t,1}.Solution.K(i)  = max((input_Water1(i,1)*states{t,1}.Solution.K(i-1)  + CurrentWater.*states{t,1}.Solution.K(i)  + input_Water6(i,1)*states{t,1}.Solution.K(i+Nx))./(input_Water1(i,1)+CurrentWater+input_Water6(i,1)),0);
        states{t,1}.Solution.C4(i) = max((input_Water1(i,1)*states{t,1}.Solution.C4(i-1) + CurrentWater.*states{t,1}.Solution.C4(i) + input_Water6(i,1)*states{t,1}.Solution.C4(i+Nx))./(input_Water1(i,1)+CurrentWater+input_Water6(i,1)),0);
        states{t,1}.Solution.S6(i) = max((input_Water1(i,1)*states{t,1}.Solution.S6(i-1) + CurrentWater.*states{t,1}.Solution.S6(i) + input_Water6(i,1)*states{t,1}.Solution.S6(i+Nx))./(input_Water1(i,1)+CurrentWater+input_Water6(i,1)),0);
        states{t,1}.Solution.S2(i) = max((input_Water1(i,1)*states{t,1}.Solution.S2(i-1) + CurrentWater.*states{t,1}.Solution.S2(i) + input_Water6(i,1)*states{t,1}.Solution.S2(i+Nx))./(input_Water1(i,1)+CurrentWater+input_Water6(i,1)),0);
        states{t,1}.Solution.Si(i) = max((input_Water1(i,1)*states{t,1}.Solution.Si(i-1) + CurrentWater.*states{t,1}.Solution.Si(i) + input_Water6(i,1)*states{t,1}.Solution.Si(i+Nx))./(input_Water1(i,1)+CurrentWater+input_Water6(i,1)),0);
        states{t,1}.Solution.Fe2(i) = max((input_Water1(i,1)*states{t,1}.Solution.Fe2(i-1) + CurrentWater.*states{t,1}.Solution.Fe2(i) + input_Water6(i,1)*states{t,1}.Solution.Fe2(i+Nx))./(input_Water1(i,1)+CurrentWater+input_Water6(i,1)),0);
        states{t,1}.Solution.Fe3(i) = max((input_Water1(i,1)*states{t,1}.Solution.Fe3(i-1) + CurrentWater.*states{t,1}.Solution.Fe3(i) + input_Water6(i,1)*states{t,1}.Solution.Fe3(i+Nx))./(input_Water1(i,1)+CurrentWater+input_Water6(i,1)),0);
        states{t,1}.Solution.Acetate(i) = max((input_Water1(i,1)*states{t,1}.Solution.Acetate(i-1) + CurrentWater.*states{t,1}.Solution.Acetate(i) + input_Water6(i,1)*states{t,1}.Solution.Acetate(i+Nx))./(input_Water1(i,1)+CurrentWater+input_Water6(i,1)),0);
    elseif (nx == 1) && (nz == Nz) % Bottom Left Cell (3)
        if flux(2) < 0
            input_Water2(i,1) = flux(2)*rhoW*dt;
        end
        if flux(5) > 0
            input_Water5(i,1) = flux(5)*rhoW*dt;
        end
        states{t,1}.Solution.Ca(i) = max((input_Water2(i,1)*states{t,1}.Solution.Ca(i+1) + CurrentWater.*states{t,1}.Solution.Ca(i) + input_Water5(i,1)*states{t,1}.Solution.Ca(i-Nx))./(input_Water2(i,1)+CurrentWater+input_Water5(i,1)),0);
        states{t,1}.Solution.Na(i) = max((input_Water2(i,1)*states{t,1}.Solution.Na(i+1) + CurrentWater.*states{t,1}.Solution.Na(i) + input_Water5(i,1)*states{t,1}.Solution.Na(i-Nx))./(input_Water2(i,1)+CurrentWater+input_Water5(i,1)),0);
        states{t,1}.Solution.Mg(i) = max((input_Water2(i,1)*states{t,1}.Solution.Mg(i+1) + CurrentWater.*states{t,1}.Solution.Mg(i) + input_Water5(i,1)*states{t,1}.Solution.Mg(i-Nx))./(input_Water2(i,1)+CurrentWater+input_Water5(i,1)),0);
        states{t,1}.Solution.Cl(i) = max((input_Water2(i,1)*states{t,1}.Solution.Cl(i+1) + CurrentWater.*states{t,1}.Solution.Cl(i) + input_Water5(i,1)*states{t,1}.Solution.Cl(i-Nx))./(input_Water2(i,1)+CurrentWater+input_Water5(i,1)),0);
        states{t,1}.Solution.K(i)  = max((input_Water2(i,1)*states{t,1}.Solution.K(i+1)  + CurrentWater.*states{t,1}.Solution.K(i)  + input_Water5(i,1)*states{t,1}.Solution.K(i-Nx))./(input_Water2(i,1)+CurrentWater+input_Water5(i,1)),0);
        states{t,1}.Solution.C4(i) = max((input_Water2(i,1)*states{t,1}.Solution.C4(i+1) + CurrentWater.*states{t,1}.Solution.C4(i) + input_Water5(i,1)*states{t,1}.Solution.C4(i-Nx))./(input_Water2(i,1)+CurrentWater+input_Water5(i,1)),0);
        states{t,1}.Solution.S6(i) = max((input_Water2(i,1)*states{t,1}.Solution.S6(i+1) + CurrentWater.*states{t,1}.Solution.S6(i) + input_Water5(i,1)*states{t,1}.Solution.S6(i-Nx))./(input_Water2(i,1)+CurrentWater+input_Water5(i,1)),0);
        states{t,1}.Solution.S2(i) = max((input_Water2(i,1)*states{t,1}.Solution.S2(i+1) + CurrentWater.*states{t,1}.Solution.S2(i) + input_Water5(i,1)*states{t,1}.Solution.S2(i-Nx))./(input_Water2(i,1)+CurrentWater+input_Water5(i,1)),0);
        states{t,1}.Solution.Si(i) = max((input_Water2(i,1)*states{t,1}.Solution.Si(i+1) + CurrentWater.*states{t,1}.Solution.Si(i) + input_Water5(i,1)*states{t,1}.Solution.Si(i-Nx))./(input_Water2(i,1)+CurrentWater+input_Water5(i,1)),0);
        states{t,1}.Solution.Fe2(i) = max((input_Water2(i,1)*states{t,1}.Solution.Fe2(i+1) + CurrentWater.*states{t,1}.Solution.Fe2(i) + input_Water5(i,1)*states{t,1}.Solution.Fe2(i-Nx))./(input_Water2(i,1)+CurrentWater+input_Water5(i,1)),0);
        states{t,1}.Solution.Fe3(i) = max((input_Water2(i,1)*states{t,1}.Solution.Fe3(i+1) + CurrentWater.*states{t,1}.Solution.Fe3(i) + input_Water5(i,1)*states{t,1}.Solution.Fe3(i-Nx))./(input_Water2(i,1)+CurrentWater+input_Water5(i,1)),0);
        states{t,1}.Solution.Acetate(i) = max((input_Water2(i,1)*states{t,1}.Solution.Acetate(i+1) + CurrentWater.*states{t,1}.Solution.Acetate(i) + input_Water6(i,1)*states{t,1}.Solution.Acetate(i-Nx))./(input_Water2(i,1)+CurrentWater+input_Water5(i,1)),0);
    elseif (nx == Nx) && (nz == Nz) % Bottom Right Cell (4)
        if flux(1) > 0
            input_Water1(i,1) = flux(1)*rhoW*dt;
        end
        if flux(5) > 0
            input_Water5(i,1) = flux(5)*rhoW*dt;
        end
        states{t,1}.Solution.Ca(i) = max((input_Water1(i,1)*states{t,1}.Solution.Ca(i-1) + CurrentWater.*states{t,1}.Solution.Ca(i) + input_Water5(i,1)*states{t,1}.Solution.Ca(i-Nx))./(input_Water1(i,1)+CurrentWater+input_Water5(i,1)),0);
        states{t,1}.Solution.Na(i) = max((input_Water1(i,1)*states{t,1}.Solution.Na(i-1) + CurrentWater.*states{t,1}.Solution.Na(i) + input_Water5(i,1)*states{t,1}.Solution.Na(i-Nx))./(input_Water1(i,1)+CurrentWater+input_Water5(i,1)),0);
        states{t,1}.Solution.Mg(i) = max((input_Water1(i,1)*states{t,1}.Solution.Mg(i-1) + CurrentWater.*states{t,1}.Solution.Mg(i) + input_Water5(i,1)*states{t,1}.Solution.Mg(i-Nx))./(input_Water1(i,1)+CurrentWater+input_Water5(i,1)),0);
        states{t,1}.Solution.Cl(i) = max((input_Water1(i,1)*states{t,1}.Solution.Cl(i-1) + CurrentWater.*states{t,1}.Solution.Cl(i) + input_Water5(i,1)*states{t,1}.Solution.Cl(i-Nx))./(input_Water1(i,1)+CurrentWater+input_Water5(i,1)),0);
        states{t,1}.Solution.K(i)  = max((input_Water1(i,1)*states{t,1}.Solution.K(i-1)  + CurrentWater.*states{t,1}.Solution.K(i)  + input_Water5(i,1)*states{t,1}.Solution.K(i-Nx))./(input_Water1(i,1)+CurrentWater+input_Water5(i,1)),0);
        states{t,1}.Solution.C4(i) = max((input_Water1(i,1)*states{t,1}.Solution.C4(i-1) + CurrentWater.*states{t,1}.Solution.C4(i) + input_Water5(i,1)*states{t,1}.Solution.C4(i-Nx))./(input_Water1(i,1)+CurrentWater+input_Water5(i,1)),0);
        states{t,1}.Solution.S6(i) = max((input_Water1(i,1)*states{t,1}.Solution.S6(i-1) + CurrentWater.*states{t,1}.Solution.S6(i) + input_Water5(i,1)*states{t,1}.Solution.S6(i-Nx))./(input_Water1(i,1)+CurrentWater+input_Water5(i,1)),0);
        states{t,1}.Solution.S2(i) = max((input_Water1(i,1)*states{t,1}.Solution.S2(i-1) + CurrentWater.*states{t,1}.Solution.S2(i) + input_Water5(i,1)*states{t,1}.Solution.S2(i-Nx))./(input_Water1(i,1)+CurrentWater+input_Water5(i,1)),0);
        states{t,1}.Solution.Si(i) = max((input_Water1(i,1)*states{t,1}.Solution.Si(i-1) + CurrentWater.*states{t,1}.Solution.Si(i) + input_Water5(i,1)*states{t,1}.Solution.Si(i-Nx))./(input_Water1(i,1)+CurrentWater+input_Water5(i,1)),0);
        states{t,1}.Solution.Fe2(i) = max((input_Water1(i,1)*states{t,1}.Solution.Fe2(i-1) + CurrentWater.*states{t,1}.Solution.Fe2(i) + input_Water5(i,1)*states{t,1}.Solution.Fe2(i-Nx))./(input_Water1(i,1)+CurrentWater+input_Water5(i,1)),0);
        states{t,1}.Solution.Fe3(i) = max((input_Water1(i,1)*states{t,1}.Solution.Fe3(i-1) + CurrentWater.*states{t,1}.Solution.Fe3(i) + input_Water5(i,1)*states{t,1}.Solution.Fe3(i-Nx))./(input_Water1(i,1)+CurrentWater+input_Water5(i,1)),0);
        states{t,1}.Solution.Acetate(i) = max((input_Water2(i,1)*states{t,1}.Solution.Acetate(i-1) + CurrentWater.*states{t,1}.Solution.Acetate(i) + input_Water6(i,1)*states{t,1}.Solution.Acetate(i-Nx))./(input_Water1(i,1)+CurrentWater+input_Water5(i,1)),0);
    elseif (2 <= nx) && (nx <= Nx-1) && (nz == 1) % Top bar (5)
        if flux(1) > 0
            input_Water1(i,1) = flux(1)*rhoW*dt;
        end
        if flux(2) < 0
            input_Water2(i,1) = flux(2)*rhoW*dt;
        end
        if flux(6) < 0
            input_Water6(i,1) = flux(6)*rhoW*dt;
        end
        states{t,1}.Solution.Ca(i) = max((input_Water1(i,1)*states{t,1}.Solution.Ca(i-1) + input_Water2(i,1)*states{t,1}.Solution.Ca(i+1) + CurrentWater.*states{t,1}.Solution.Ca(i) + input_Water6(i,1)*states{t,1}.Solution.Ca(i+Nx))./(input_Water1(i,1) + input_Water2(i,1) + CurrentWater + input_Water6(i,1)),0);
        states{t,1}.Solution.Na(i) = max((input_Water1(i,1)*states{t,1}.Solution.Na(i-1) +  input_Water2(i,1)*states{t,1}.Solution.Na(i+1) + CurrentWater.*states{t,1}.Solution.Na(i) + input_Water6(i,1)*states{t,1}.Solution.Na(i+Nx))./(input_Water1(i,1) + input_Water2(i,1)+CurrentWater+input_Water6(i,1)),0);
        states{t,1}.Solution.Mg(i) = max((input_Water1(i,1)*states{t,1}.Solution.Mg(i-1) + input_Water2(i,1)*states{t,1}.Solution.Mg(i+1) + CurrentWater.*states{t,1}.Solution.Mg(i) + input_Water6(i,1)*states{t,1}.Solution.Mg(i+Nx))./(input_Water1(i,1) + input_Water2(i,1)+CurrentWater+input_Water6(i,1)),0);
        states{t,1}.Solution.Cl(i) = max((input_Water1(i,1)*states{t,1}.Solution.Cl(i-1) + input_Water2(i,1)*states{t,1}.Solution.Cl(i+1) + CurrentWater.*states{t,1}.Solution.Cl(i) + input_Water6(i,1)*states{t,1}.Solution.Cl(i+Nx))./(input_Water1(i,1) + input_Water2(i,1)+CurrentWater+input_Water6(i,1)),0);
        states{t,1}.Solution.K(i)  = max((input_Water1(i,1)*states{t,1}.Solution.K(i-1)  + input_Water2(i,1)*states{t,1}.Solution.K(i+1)  + CurrentWater.*states{t,1}.Solution.K(i)  + input_Water6(i,1)*states{t,1}.Solution.K(i+Nx))./(input_Water1(i,1) + input_Water2(i,1)+CurrentWater+input_Water6(i,1)),0);
        states{t,1}.Solution.C4(i) = max((input_Water1(i,1)*states{t,1}.Solution.C4(i-1) + input_Water2(i,1)*states{t,1}.Solution.C4(i+1) + CurrentWater.*states{t,1}.Solution.C4(i) + input_Water6(i,1)*states{t,1}.Solution.C4(i+Nx))./(input_Water1(i,1) + input_Water2(i,1)+CurrentWater+input_Water6(i,1)),0);
        states{t,1}.Solution.S6(i) = max((input_Water1(i,1)*states{t,1}.Solution.S6(i-1) + input_Water2(i,1)*states{t,1}.Solution.S6(i+1) + CurrentWater.*states{t,1}.Solution.S6(i) + input_Water6(i,1)*states{t,1}.Solution.S6(i+Nx))./(input_Water1(i,1) + input_Water2(i,1)+CurrentWater+input_Water6(i,1)),0);
        states{t,1}.Solution.S2(i) = max((input_Water1(i,1)*states{t,1}.Solution.S2(i-1) + input_Water2(i,1)*states{t,1}.Solution.S2(i+1) + CurrentWater.*states{t,1}.Solution.S2(i) + input_Water6(i,1)*states{t,1}.Solution.S2(i+Nx))./(input_Water1(i,1) + input_Water2(i,1)+CurrentWater+input_Water6(i,1)),0);
        states{t,1}.Solution.Si(i) = max((input_Water1(i,1)*states{t,1}.Solution.Si(i-1) + input_Water2(i,1)*states{t,1}.Solution.Si(i+1) + CurrentWater.*states{t,1}.Solution.Si(i) + input_Water6(i,1)*states{t,1}.Solution.Si(i+Nx))./(input_Water1(i,1) + input_Water2(i,1)+CurrentWater+input_Water6(i,1)),0);
        states{t,1}.Solution.Fe2(i) = max((input_Water1(i,1)*states{t,1}.Solution.Fe2(i-1) + input_Water2(i,1)*states{t,1}.Solution.Fe2(i+1) + CurrentWater.*states{t,1}.Solution.Fe2(i) + input_Water6(i,1)*states{t,1}.Solution.Fe2(i+Nx))./(input_Water1(i,1) + input_Water2(i,1)+CurrentWater+input_Water6(i,1)),0);
        states{t,1}.Solution.Fe3(i) = max((input_Water1(i,1)*states{t,1}.Solution.Fe3(i-1) + input_Water2(i,1)*states{t,1}.Solution.Fe3(i+1) + CurrentWater.*states{t,1}.Solution.Fe3(i) + input_Water6(i,1)*states{t,1}.Solution.Fe3(i+Nx))./(input_Water1(i,1) + input_Water2(i,1)+CurrentWater+input_Water6(i,1)),0);
        states{t,1}.Solution.Acetate(i) = max((input_Water1(i,1)*states{t,1}.Solution.Acetate(i-1) + input_Water2(i,1)*states{t,1}.Solution.Acetate(i+1) + CurrentWater.*states{t,1}.Solution.Acetate(i) + input_Water6(i,1)*states{t,1}.Solution.Acetate(i+Nx))./(input_Water1(i,1) + input_Water2(i,1)+CurrentWater+input_Water6(i,1)),0);
    elseif (2 <= nz) && (nz <= Nz - 1) && (nx == Nx) % Right bar (6)
        if flux(1) > 0
            input_Water1(i,1) = flux(1)*rhoW*dt;
        end
        if flux(5) > 0
            input_Water5(i,1) = flux(5)*rhoW*dt;
        end
        if flux(6) < 0
            input_Water6(i,1) = flux(6)*rhoW*dt;
        end
        states{t,1}.Solution.Ca(i) = max((input_Water1(i,1)*states{t,1}.Solution.Ca(i-1) + input_Water5(i,1)*states{t,1}.Solution.Ca(i-Nx) + CurrentWater.*states{t,1}.Solution.Ca(i) + input_Water6(i,1)*states{t,1}.Solution.Ca(i+Nx))./(input_Water1(i,1) + input_Water5(i,1)+CurrentWater+input_Water6(i,1)),0);
        states{t,1}.Solution.Na(i) = max((input_Water1(i,1)*states{t,1}.Solution.Na(i-1) +  input_Water5(i,1)*states{t,1}.Solution.Na(i-Nx) + CurrentWater.*states{t,1}.Solution.Na(i) + input_Water6(i,1)*states{t,1}.Solution.Na(i+Nx))./(input_Water1(i,1) + input_Water5(i,1)+CurrentWater+input_Water6(i,1)),0);
        states{t,1}.Solution.Mg(i) = max((input_Water1(i,1)*states{t,1}.Solution.Mg(i-1) + input_Water5(i,1)*states{t,1}.Solution.Mg(i-Nx) + CurrentWater.*states{t,1}.Solution.Mg(i) + input_Water6(i,1)*states{t,1}.Solution.Mg(i+Nx))./(input_Water1(i,1) + input_Water5(i,1)+CurrentWater+input_Water6(i,1)),0);
        states{t,1}.Solution.Cl(i) = max((input_Water1(i,1)*states{t,1}.Solution.Cl(i-1) + input_Water5(i,1)*states{t,1}.Solution.Cl(i-Nx) + CurrentWater.*states{t,1}.Solution.Cl(i) + input_Water6(i,1)*states{t,1}.Solution.Cl(i+Nx))./(input_Water1(i,1) + input_Water5(i,1)+CurrentWater+input_Water6(i,1)),0);
        states{t,1}.Solution.K(i)  = max((input_Water1(i,1)*states{t,1}.Solution.K(i-1)  + input_Water5(i,1)*states{t,1}.Solution.K(i-Nx)  + CurrentWater.*states{t,1}.Solution.K(i)  + input_Water6(i,1)*states{t,1}.Solution.K(i+Nx))./(input_Water1(i,1) + input_Water5(i,1)+CurrentWater+input_Water6(i,1)),0);
        states{t,1}.Solution.C4(i) = max((input_Water1(i,1)*states{t,1}.Solution.C4(i-1) + input_Water5(i,1)*states{t,1}.Solution.C4(i-Nx) + CurrentWater.*states{t,1}.Solution.C4(i) + input_Water6(i,1)*states{t,1}.Solution.C4(i+Nx))./(input_Water1(i,1) + input_Water5(i,1)+CurrentWater+input_Water6(i,1)),0);
        states{t,1}.Solution.S6(i) = max((input_Water1(i,1)*states{t,1}.Solution.S6(i-1) + input_Water5(i,1)*states{t,1}.Solution.S6(i-Nx) + CurrentWater.*states{t,1}.Solution.S6(i) + input_Water6(i,1)*states{t,1}.Solution.S6(i+Nx))./(input_Water1(i,1) + input_Water5(i,1)+CurrentWater+input_Water6(i,1)),0);
        states{t,1}.Solution.S2(i) = max((input_Water1(i,1)*states{t,1}.Solution.S2(i-1) + input_Water5(i,1)*states{t,1}.Solution.S2(i-Nx) + CurrentWater.*states{t,1}.Solution.S2(i) + input_Water6(i,1)*states{t,1}.Solution.S2(i+Nx))./(input_Water1(i,1) + input_Water5(i,1)+CurrentWater+input_Water6(i,1)),0);
        states{t,1}.Solution.Si(i) = max((input_Water1(i,1)*states{t,1}.Solution.Si(i-1) + input_Water5(i,1)*states{t,1}.Solution.Si(i-Nx) + CurrentWater.*states{t,1}.Solution.Si(i) + input_Water6(i,1)*states{t,1}.Solution.Si(i+Nx))./(input_Water1(i,1) + input_Water5(i,1)+CurrentWater+input_Water6(i,1)),0);
        states{t,1}.Solution.Fe2(i) = max((input_Water1(i,1)*states{t,1}.Solution.Fe2(i-1) + input_Water5(i,1)*states{t,1}.Solution.Fe2(i-Nx) + CurrentWater.*states{t,1}.Solution.Fe2(i) + input_Water6(i,1)*states{t,1}.Solution.Fe2(i+Nx))./(input_Water1(i,1) + input_Water5(i,1)+CurrentWater+input_Water6(i,1)),0);
        states{t,1}.Solution.Fe3(i) = max((input_Water1(i,1)*states{t,1}.Solution.Fe3(i-1) + input_Water5(i,1)*states{t,1}.Solution.Fe3(i-Nx) + CurrentWater.*states{t,1}.Solution.Fe3(i) + input_Water6(i,1)*states{t,1}.Solution.Fe3(i+Nx))./(input_Water1(i,1) + input_Water5(i,1)+CurrentWater+input_Water6(i,1)),0);
        states{t,1}.Solution.Acetate(i) = max((input_Water1(i,1)*states{t,1}.Solution.Acetate(i-1) + input_Water5(i,1)*states{t,1}.Solution.Acetate(i-Nx) + CurrentWater.*states{t,1}.Solution.Acetate(i) + input_Water6(i,1)*states{t,1}.Solution.Acetate(i+Nx))./(input_Water1(i,1) + input_Water5(i,1)+CurrentWater+input_Water6(i,1)),0);
    elseif (2 <= nx) && (nx <= Nx-1) && (nz == Nz) % Bottom bar (7)
        if flux(1) > 0
            input_Water1(i,1) = flux(1)*rhoW*dt;
        end
        if flux(2) < 0
            input_Water2(i,1) = flux(2)*rhoW*dt;
        end
        if flux(5) > 0
            input_Water5(i,1) = flux(5)*rhoW*dt;
        end
        states{t,1}.Solution.Ca(i) = max((input_Water1(i,1)*states{t,1}.Solution.Ca(i-1) + input_Water2(i,1)*states{t,1}.Solution.Ca(i+1) + CurrentWater.*states{t,1}.Solution.Ca(i) + input_Water5(i,1)*states{t,1}.Solution.Ca(i-Nx))./(input_Water1(i,1) + input_Water2(i,1)+CurrentWater+input_Water5(i,1)),0);
        states{t,1}.Solution.Na(i) = max((input_Water1(i,1)*states{t,1}.Solution.Na(i-1) +  input_Water2(i,1)*states{t,1}.Solution.Na(i+1) + CurrentWater.*states{t,1}.Solution.Na(i) + input_Water5(i,1)*states{t,1}.Solution.Na(i-Nx))./(input_Water1(i,1) + input_Water2(i,1)+CurrentWater+input_Water5(i,1)),0);
        states{t,1}.Solution.Mg(i) = max((input_Water1(i,1)*states{t,1}.Solution.Mg(i-1) + input_Water2(i,1)*states{t,1}.Solution.Mg(i+1) + CurrentWater.*states{t,1}.Solution.Mg(i) + input_Water5(i,1)*states{t,1}.Solution.Mg(i-Nx))./(input_Water1(i,1) + input_Water2(i,1)+CurrentWater+input_Water5(i,1)),0);
        states{t,1}.Solution.Cl(i) = max((input_Water1(i,1)*states{t,1}.Solution.Cl(i-1) + input_Water2(i,1)*states{t,1}.Solution.Cl(i+1) + CurrentWater.*states{t,1}.Solution.Cl(i) + input_Water5(i,1)*states{t,1}.Solution.Cl(i-Nx))./(input_Water1(i,1) + input_Water2(i,1)+CurrentWater+input_Water5(i,1)),0);
        states{t,1}.Solution.K(i)  = max((input_Water1(i,1)*states{t,1}.Solution.K(i-1)  + input_Water2(i,1)*states{t,1}.Solution.K(i+1)  + CurrentWater.*states{t,1}.Solution.K(i)  + input_Water5(i,1)*states{t,1}.Solution.K(i-Nx))./(input_Water1(i,1) + input_Water2(i,1)+CurrentWater+input_Water5(i,1)),0);
        states{t,1}.Solution.C4(i) = max((input_Water1(i,1)*states{t,1}.Solution.C4(i-1) + input_Water2(i,1)*states{t,1}.Solution.C4(i+1) + CurrentWater.*states{t,1}.Solution.C4(i) + input_Water5(i,1)*states{t,1}.Solution.C4(i-Nx))./(input_Water1(i,1) + input_Water2(i,1)+CurrentWater+input_Water5(i,1)),0);
        states{t,1}.Solution.S6(i) = max((input_Water1(i,1)*states{t,1}.Solution.S6(i-1) + input_Water2(i,1)*states{t,1}.Solution.S6(i+1) + CurrentWater.*states{t,1}.Solution.S6(i) + input_Water5(i,1)*states{t,1}.Solution.S6(i-Nx))./(input_Water1(i,1) + input_Water2(i,1)+CurrentWater+input_Water5(i,1)),0);
        states{t,1}.Solution.S2(i) = max((input_Water1(i,1)*states{t,1}.Solution.S2(i-1) + input_Water2(i,1)*states{t,1}.Solution.S2(i+1) + CurrentWater.*states{t,1}.Solution.S2(i) + input_Water5(i,1)*states{t,1}.Solution.S2(i-Nx))./(input_Water1(i,1) + input_Water2(i,1)+CurrentWater+input_Water5(i,1)),0);
        states{t,1}.Solution.Si(i) = max((input_Water1(i,1)*states{t,1}.Solution.Si(i-1) + input_Water2(i,1)*states{t,1}.Solution.Si(i+1) + CurrentWater.*states{t,1}.Solution.Si(i) + input_Water5(i,1)*states{t,1}.Solution.Si(i-Nx))./(input_Water1(i,1) + input_Water2(i,1)+CurrentWater+input_Water5(i,1)),0);
        states{t,1}.Solution.Fe2(i) = max((input_Water1(i,1)*states{t,1}.Solution.Fe2(i-1) + input_Water2(i,1)*states{t,1}.Solution.Fe2(i+1) + CurrentWater.*states{t,1}.Solution.Fe2(i) + input_Water5(i,1)*states{t,1}.Solution.Fe2(i-Nx))./(input_Water1(i,1) + input_Water2(i,1)+CurrentWater+input_Water5(i,1)),0);
        states{t,1}.Solution.Fe3(i) = max((input_Water1(i,1)*states{t,1}.Solution.Fe3(i-1) + input_Water2(i,1)*states{t,1}.Solution.Fe3(i+1) + CurrentWater.*states{t,1}.Solution.Fe3(i) + input_Water5(i,1)*states{t,1}.Solution.Fe3(i-Nx))./(input_Water1(i,1) + input_Water2(i,1)+CurrentWater+input_Water5(i,1)),0);
        states{t,1}.Solution.Acetate(i) = max((input_Water1(i,1)*states{t,1}.Solution.Acetate(i-1) + input_Water2(i,1)*states{t,1}.Solution.Acetate(i+1) + CurrentWater.*states{t,1}.Solution.Acetate(i) + input_Water5(i,1)*states{t,1}.Solution.Acetate(i-Nx))./(input_Water1(i,1) + input_Water2(i,1)+CurrentWater+input_Water5(i,1)),0);
    elseif (2 <= nz) && (nz <= Nz - 1) && (nx == 1) % Left bar (8)
        if flux(2) < 0
            input_Water2(i,1) = flux(2)*rhoW*dt;
        end
        if flux(5) > 0
            input_Water5(i,1) = flux(5)*rhoW*dt;
        end
        if flux(6) < 0
            input_Water6(i,1) = flux(6)*rhoW*dt;
        end
        states{t,1}.Solution.Ca(i) = max((input_Water2(i,1)*states{t,1}.Solution.Ca(i+1) + input_Water5(i,1)*states{t,1}.Solution.Ca(i-Nx) + CurrentWater.*states{t,1}.Solution.Ca(i) + input_Water6(i,1)*states{t,1}.Solution.Ca(i+Nx))./(input_Water2(i,1)+CurrentWater+input_Water5(i,1) + input_Water6(i,1)),0);
        states{t,1}.Solution.Na(i) = max((input_Water2(i,1)*states{t,1}.Solution.Na(i+1) +  input_Water5(i,1)*states{t,1}.Solution.Na(i-Nx) + CurrentWater.*states{t,1}.Solution.Na(i) + input_Water6(i,1)*states{t,1}.Solution.Na(i+Nx))./(input_Water2(i,1)+CurrentWater+input_Water5(i,1) + input_Water6(i,1)),0);
        states{t,1}.Solution.Mg(i) = max((input_Water2(i,1)*states{t,1}.Solution.Mg(i+1) + input_Water5(i,1)*states{t,1}.Solution.Mg(i-Nx) + CurrentWater.*states{t,1}.Solution.Mg(i) + input_Water6(i,1)*states{t,1}.Solution.Mg(i+Nx))./(input_Water2(i,1)+CurrentWater+input_Water5(i,1) + input_Water6(i,1)),0);
        states{t,1}.Solution.Cl(i) = max((input_Water2(i,1)*states{t,1}.Solution.Cl(i+1) + input_Water5(i,1)*states{t,1}.Solution.Cl(i-Nx) + CurrentWater.*states{t,1}.Solution.Cl(i) + input_Water6(i,1)*states{t,1}.Solution.Cl(i+Nx))./(input_Water2(i,1)+CurrentWater+input_Water5(i,1) + input_Water6(i,1)),0);
        states{t,1}.Solution.K(i)  = max((input_Water2(i,1)*states{t,1}.Solution.K(i+1)  + input_Water5(i,1)*states{t,1}.Solution.K(i-Nx)  + CurrentWater.*states{t,1}.Solution.K(i)  + input_Water6(i,1)*states{t,1}.Solution.K(i+Nx))./(input_Water2(i,1)+CurrentWater+input_Water5(i,1) + input_Water6(i,1)),0);
        states{t,1}.Solution.C4(i) = max((input_Water2(i,1)*states{t,1}.Solution.C4(i+1) + input_Water5(i,1)*states{t,1}.Solution.C4(i-Nx) + CurrentWater.*states{t,1}.Solution.C4(i) + input_Water6(i,1)*states{t,1}.Solution.C4(i+Nx))./(input_Water2(i,1)+CurrentWater+input_Water5(i,1) + input_Water6(i,1)),0);
        states{t,1}.Solution.S6(i) = max((input_Water2(i,1)*states{t,1}.Solution.S6(i+1) + input_Water5(i,1)*states{t,1}.Solution.S6(i-Nx) + CurrentWater.*states{t,1}.Solution.S6(i) + input_Water6(i,1)*states{t,1}.Solution.S6(i+Nx))./(input_Water2(i,1)+CurrentWater+input_Water5(i,1) + input_Water6(i,1)),0);
        states{t,1}.Solution.S2(i) = max((input_Water2(i,1)*states{t,1}.Solution.S2(i+1) + input_Water5(i,1)*states{t,1}.Solution.S2(i-Nx) + CurrentWater.*states{t,1}.Solution.S2(i) + input_Water6(i,1)*states{t,1}.Solution.S2(i+Nx))./(input_Water2(i,1)+CurrentWater+input_Water5(i,1) + input_Water6(i,1)),0);
        states{t,1}.Solution.Si(i) = max((input_Water2(i,1)*states{t,1}.Solution.Si(i+1) + input_Water5(i,1)*states{t,1}.Solution.Si(i-Nx) + CurrentWater.*states{t,1}.Solution.Si(i) + input_Water6(i,1)*states{t,1}.Solution.Si(i+Nx))./(input_Water2(i,1)+CurrentWater+input_Water5(i,1) + input_Water6(i,1)),0);
        states{t,1}.Solution.Fe2(i) = max((input_Water2(i,1)*states{t,1}.Solution.Fe2(i+1) + input_Water5(i,1)*states{t,1}.Solution.Fe2(i-Nx) + CurrentWater.*states{t,1}.Solution.Fe2(i) + input_Water6(i,1)*states{t,1}.Solution.Fe2(i+Nx))./(input_Water2(i,1)+CurrentWater+input_Water5(i,1) + input_Water6(i,1)),0);
        states{t,1}.Solution.Fe3(i) = max((input_Water2(i,1)*states{t,1}.Solution.Fe3(i+1) + input_Water5(i,1)*states{t,1}.Solution.Fe3(i-Nx) + CurrentWater.*states{t,1}.Solution.Fe3(i) + input_Water6(i,1)*states{t,1}.Solution.Fe3(i+Nx))./(input_Water2(i,1)+CurrentWater+input_Water5(i,1) + input_Water6(i,1)),0);
        states{t,1}.Solution.Acetate(i) = max((input_Water2(i,1)*states{t,1}.Solution.Acetate(i+1) + input_Water5(i,1)*states{t,1}.Solution.Acetate(i-Nx) + CurrentWater.*states{t,1}.Solution.Acetate(i) + input_Water6(i,1)*states{t,1}.Solution.Acetate(i+Nx))./(input_Water2(i,1)+CurrentWater+input_Water5(i,1) + input_Water6(i,1)),0);
    else
        if flux(1) > 0
            input_Water1(i,1) = flux(1)*rhoW*dt;
        end
        if flux(2) < 0
            input_Water2(i,1) = flux(2)*rhoW*dt;
        end
        if flux(5) > 0
            input_Water5(i,1) = flux(5)*rhoW*dt;
        end
        if flux(6) < 0
            input_Water6(i,1) = flux(6)*rhoW*dt;
        end
        states{t,1}.Solution.Ca(i) = max((input_Water1(i,1)*states{t,1}.Solution.Ca(i-1) + input_Water2(i,1)*states{t,1}.Solution.Ca(i+1) + input_Water5(i,1)*states{t,1}.Solution.Ca(i-Nx) + CurrentWater.*states{t,1}.Solution.Ca(i) + input_Water6(i,1)*states{t,1}.Solution.Ca(i+Nx))./(input_Water1(i,1) + input_Water2(i,1)+CurrentWater+input_Water5(i,1) + input_Water6(i,1)),0);
        states{t,1}.Solution.Na(i) = max((input_Water1(i,1)*states{t,1}.Solution.Na(i-1) + input_Water2(i,1)*states{t,1}.Solution.Na(i+1) +  input_Water5(i,1)*states{t,1}.Solution.Na(i-Nx) + CurrentWater.*states{t,1}.Solution.Na(i) + input_Water6(i,1)*states{t,1}.Solution.Na(i+Nx))./(input_Water1(i,1) + input_Water2(i,1)+CurrentWater+input_Water5(i,1) + input_Water6(i,1)),0);
        states{t,1}.Solution.Mg(i) = max((input_Water1(i,1)*states{t,1}.Solution.Mg(i-1) + input_Water2(i,1)*states{t,1}.Solution.Mg(i+1) + input_Water5(i,1)*states{t,1}.Solution.Mg(i-Nx) + CurrentWater.*states{t,1}.Solution.Mg(i) + input_Water6(i,1)*states{t,1}.Solution.Mg(i+Nx))./(input_Water1(i,1) + input_Water2(i,1)+CurrentWater+input_Water5(i,1) + input_Water6(i,1)),0);
        states{t,1}.Solution.Cl(i) = max((input_Water1(i,1)*states{t,1}.Solution.Cl(i-1) + input_Water2(i,1)*states{t,1}.Solution.Cl(i+1) + input_Water5(i,1)*states{t,1}.Solution.Cl(i-Nx) + CurrentWater.*states{t,1}.Solution.Cl(i) + input_Water6(i,1)*states{t,1}.Solution.Cl(i+Nx))./(input_Water1(i,1) + input_Water2(i,1)+CurrentWater+input_Water5(i,1) + input_Water6(i,1)),0);
        states{t,1}.Solution.K(i)  = max((input_Water1(i,1)*states{t,1}.Solution.K(i-1)  + input_Water2(i,1)*states{t,1}.Solution.K(i+1)  + input_Water5(i,1)*states{t,1}.Solution.K(i-Nx)  + CurrentWater.*states{t,1}.Solution.K(i)  + input_Water6(i,1)*states{t,1}.Solution.K(i+Nx))./(input_Water1(i,1) + input_Water2(i,1)+CurrentWater+input_Water5(i,1) + input_Water6(i,1)),0);
        states{t,1}.Solution.C4(i) = max((input_Water1(i,1)*states{t,1}.Solution.C4(i-1) + input_Water2(i,1)*states{t,1}.Solution.C4(i+1) + input_Water5(i,1)*states{t,1}.Solution.C4(i-Nx) + CurrentWater.*states{t,1}.Solution.C4(i) + input_Water6(i,1)*states{t,1}.Solution.C4(i+Nx))./(input_Water1(i,1) + input_Water2(i,1)+CurrentWater+input_Water5(i,1) + input_Water6(i,1)),0);
        states{t,1}.Solution.S6(i) = max((input_Water1(i,1)*states{t,1}.Solution.S6(i-1) + input_Water2(i,1)*states{t,1}.Solution.S6(i+1) + input_Water5(i,1)*states{t,1}.Solution.S6(i-Nx) + CurrentWater.*states{t,1}.Solution.S6(i) + input_Water6(i,1)*states{t,1}.Solution.S6(i+Nx))./(input_Water1(i,1) + input_Water2(i,1)+CurrentWater+input_Water5(i,1) + input_Water6(i,1)),0);
        states{t,1}.Solution.S2(i) = max((input_Water1(i,1)*states{t,1}.Solution.S2(i-1) + input_Water2(i,1)*states{t,1}.Solution.S2(i+1) + input_Water5(i,1)*states{t,1}.Solution.S2(i-Nx) + CurrentWater.*states{t,1}.Solution.S2(i) + input_Water6(i,1)*states{t,1}.Solution.S2(i+Nx))./(input_Water1(i,1) + input_Water2(i,1)+CurrentWater+input_Water5(i,1) + input_Water6(i,1)),0);
        states{t,1}.Solution.Si(i) = max((input_Water1(i,1)*states{t,1}.Solution.Si(i-1) + input_Water2(i,1)*states{t,1}.Solution.Si(i+1) + input_Water5(i,1)*states{t,1}.Solution.Si(i-Nx) + CurrentWater.*states{t,1}.Solution.Si(i) + input_Water6(i,1)*states{t,1}.Solution.Si(i+Nx))./(input_Water1(i,1) + input_Water2(i,1)+CurrentWater+input_Water5(i,1) + input_Water6(i,1)),0);
        states{t,1}.Solution.Fe2(i) = max((input_Water1(i,1)*states{t,1}.Solution.Fe2(i-1) + input_Water2(i,1)*states{t,1}.Solution.Fe2(i+1) + input_Water5(i,1)*states{t,1}.Solution.Fe2(i-Nx) + CurrentWater.*states{t,1}.Solution.Fe2(i) + input_Water6(i,1)*states{t,1}.Solution.Fe2(i+Nx))./(input_Water1(i,1) + input_Water2(i,1)+CurrentWater+input_Water5(i,1) + input_Water6(i,1)),0);
        states{t,1}.Solution.Fe3(i) = max((input_Water1(i,1)*states{t,1}.Solution.Fe3(i-1) + input_Water2(i,1)*states{t,1}.Solution.Fe3(i+1) + input_Water5(i,1)*states{t,1}.Solution.Fe3(i-Nx) + CurrentWater.*states{t,1}.Solution.Fe3(i) + input_Water6(i,1)*states{t,1}.Solution.Fe3(i+Nx))./(input_Water1(i,1) + input_Water2(i,1)+CurrentWater+input_Water5(i,1) + input_Water6(i,1)),0);
        states{t,1}.Solution.Acetate(i) = max((input_Water1(i,1)*states{t,1}.Solution.Acetate(i-1) + input_Water2(i,1)*states{t,1}.Solution.Acetate(i+1) + input_Water5(i,1)*states{t,1}.Solution.Acetate(i-Nx) + CurrentWater.*states{t,1}.Solution.Acetate(i) + input_Water6(i,1)*states{t,1}.Solution.Acetate(i+Nx))./(input_Water1(i,1) + input_Water2(i,1)+CurrentWater+input_Water5(i,1) + input_Water6(i,1)),0);
    end

end
end