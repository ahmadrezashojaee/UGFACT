function states = ExplicitIonTransport_1D_New(model,states,t,dt)

G = model.G;
Water_Density = states{t,1}.PVTProps.Density{1,1};
input_Water1    = zeros(G.cells.num,1);
input_Water2    = zeros(G.cells.num,1);
for i=1:G.cells.num %calculation for each cell in 1D system
    
    [cf,~] = gridCellFaces(G,i);
  
    for j=1:2 %calculation for 2 faces in a 1D system j=1 : left face, j==2 right face
        flux = states{t,1}.flux(cf(j),1);
        %Area = G.faces.areas(cf(j));
        rhoW = Water_Density(i);
        if j==1 %Left
            if i==1
                input_Water1(i,1) = 0;  
            elseif i==G.cells.num
                if flux>=0
                    input_Water1(i,1) = flux*rhoW*dt;
                end
            else
                if flux>=0
                    input_Water1(i,1) = flux*rhoW*dt;
                end
            end
        elseif j==2 %right face
            if i==1
                if flux<=0
                    input_Water2(i,1) = flux*rhoW*dt;
                end
            elseif i==G.cells.num
                input_Water2(i,1) = 0;  
            else
                if flux<=0
                    input_Water2(i,1) = flux*rhoW*dt;
                end
                
            end
        end
    end
end
for i=1:G.cells.num
    CurrentWater = states{t,1}.FlowProps.ComponentTotalMass{1,1};
    CurrentWater = CurrentWater(i);
    if i==1
        states{t,1}.Solution.Ca(i) = max((input_Water2(i,1)*states{t,1}.Solution.Ca(i+1) + CurrentWater.*states{t,1}.Solution.Ca(i))./(input_Water2(i,1)+CurrentWater),0);
        states{t,1}.IonMole.Ca(i) = max((input_Water2(i,1)*states{t,1}.IonMole.Ca(i+1) + CurrentWater.*states{t,1}.IonMole.Ca(i))./(input_Water2(i,1)+CurrentWater),0);
        states{t,1}.Solution.Mg(i) = max((input_Water2(i,1)*states{t,1}.Solution.Mg(i+1) + CurrentWater.*states{t,1}.Solution.Mg(i))./(input_Water2(i,1)+CurrentWater),0);
        states{t,1}.IonMole.Mg(i) = max((input_Water2(i,1)*states{t,1}.IonMole.Mg(i+1) + CurrentWater.*states{t,1}.IonMole.Mg(i))./(input_Water2(i,1)+CurrentWater),0);
        states{t,1}.Solution.Na(i) = max((input_Water2(i,1)*states{t,1}.Solution.Na(i+1) + CurrentWater.*states{t,1}.Solution.Na(i))./(input_Water2(i,1)+CurrentWater),0);
        states{t,1}.IonMole.Na(i) = max((input_Water2(i,1)*states{t,1}.IonMole.Na(i+1) + CurrentWater.*states{t,1}.IonMole.Na(i))./(input_Water2(i,1)+CurrentWater),0);
        states{t,1}.Solution.K(i)  = max((input_Water2(i,1)*states{t,1}.Solution.K(i+1) + CurrentWater.*states{t,1}.Solution.K(i))./(input_Water2(i,1)+CurrentWater),0);
        states{t,1}.IonMole.K(i)  = max((input_Water2(i,1)*states{t,1}.IonMole.K(i+1) + CurrentWater.*states{t,1}.IonMole.K(i))./(input_Water2(i,1)+CurrentWater),0);
        states{t,1}.Solution.Cl(i) = max((input_Water2(i,1)*states{t,1}.Solution.Cl(i+1) + CurrentWater.*states{t,1}.Solution.Cl(i))./(input_Water2(i,1)+CurrentWater),0);
        states{t,1}.IonMole.Cl(i) = max((input_Water2(i,1)*states{t,1}.IonMole.Cl(i+1) + CurrentWater.*states{t,1}.IonMole.Cl(i))./(input_Water2(i,1)+CurrentWater),0);
        states{t,1}.Solution.S6(i) = max((input_Water2(i,1)*states{t,1}.Solution.S6(i+1) + CurrentWater.*states{t,1}.Solution.S6(i))./(input_Water2(i,1)+CurrentWater),0);
        states{t,1}.IonMole.S6(i) = max((input_Water2(i,1)*states{t,1}.IonMole.S6(i+1) + CurrentWater.*states{t,1}.IonMole.S6(i))./(input_Water2(i,1)+CurrentWater),0);
        states{t,1}.Solution.S2(i) = max((input_Water2(i,1)*states{t,1}.Solution.S2(i+1) + CurrentWater.*states{t,1}.Solution.S2(i))./(input_Water2(i,1)+CurrentWater),0);
        states{t,1}.IonMole.S2(i) = max((input_Water2(i,1)*states{t,1}.IonMole.S2(i+1) + CurrentWater.*states{t,1}.IonMole.S2(i))./(input_Water2(i,1)+CurrentWater),0);
        states{t,1}.Solution.Fe3(i) = max((input_Water2(i,1)*states{t,1}.Solution.Fe3(i+1) + CurrentWater.*states{t,1}.Solution.Fe3(i))./(input_Water2(i,1)+CurrentWater),0);
        states{t,1}.IonMole.Fe3(i) = max((input_Water2(i,1)*states{t,1}.IonMole.Fe3(i+1) + CurrentWater.*states{t,1}.IonMole.Fe3(i))./(input_Water2(i,1)+CurrentWater),0);
        states{t,1}.Solution.Fe2(i) = max((input_Water2(i,1)*states{t,1}.Solution.Fe2(i+1) + CurrentWater.*states{t,1}.Solution.Fe2(i))./(input_Water2(i,1)+CurrentWater),0);
        states{t,1}.IonMole.Fe2(i) = max((input_Water2(i,1)*states{t,1}.IonMole.Fe2(i+1) + CurrentWater.*states{t,1}.IonMole.Fe2(i))./(input_Water2(i,1)+CurrentWater),0);
        states{t,1}.Solution.Si(i) = max((input_Water2(i,1)*states{t,1}.Solution.Si(i+1) + CurrentWater.*states{t,1}.Solution.Si(i))./(input_Water2(i,1)+CurrentWater),0);
        states{t,1}.IonMole.Si(i) = max((input_Water2(i,1)*states{t,1}.IonMole.Si(i+1) + CurrentWater.*states{t,1}.IonMole.Si(i))./(input_Water2(i,1)+CurrentWater),0);
        states{t,1}.Solution.Acetate(i) = max((input_Water2(i,1)*states{t,1}.Solution.Acetate(i+1) + CurrentWater.*states{t,1}.Solution.Acetate(i))./(input_Water2(i,1)+CurrentWater),0);
        states{t,1}.IonMole.Acetate(i) = max((input_Water2(i,1)*states{t,1}.IonMole.Acetate(i+1) + CurrentWater.*states{t,1}.IonMole.Acetate(i))./(input_Water2(i,1)+CurrentWater),0);
        states{t,1}.Solution.C4(i) = max((input_Water2(i,1)*states{t,1}.Solution.C4(i+1) + CurrentWater.*states{t,1}.Solution.C4(i))./(input_Water2(i,1)+CurrentWater),0);
        states{t,1}.IonMole.C4(i) = max((input_Water2(i,1)*states{t,1}.Solution.C4(i+1) + CurrentWater.*states{t,1}.IonMole.C4(i))./(input_Water2(i,1)+CurrentWater),0);
    elseif i==G.cells.num
        states{t,1}.Solution.Ca(i) = max((input_Water1(i,1)*states{t,1}.Solution.Ca(i-1) + CurrentWater.*states{t,1}.Solution.Ca(i))./(input_Water1(i,1)+CurrentWater),0);
        states{t,1}.Solution.Mg(i) = max((input_Water1(i,1)*states{t,1}.Solution.Mg(i-1) + CurrentWater.*states{t,1}.Solution.Mg(i))./(input_Water1(i,1)+CurrentWater),0);
        states{t,1}.Solution.Na(i) = max((input_Water1(i,1)*states{t,1}.Solution.Na(i-1) + CurrentWater.*states{t,1}.Solution.Na(i))./(input_Water1(i,1)+CurrentWater),0);
        states{t,1}.Solution.K(i) = max((input_Water1(i,1)*states{t,1}.Solution.K(i-1) + CurrentWater.*states{t,1}.Solution.K(i))./(input_Water1(i,1)+CurrentWater),0);
        states{t,1}.Solution.Cl(i) = max((input_Water1(i,1)*states{t,1}.Solution.Cl(i-1) + CurrentWater.*states{t,1}.Solution.Cl(i))./(input_Water1(i,1)+CurrentWater),0);
        states{t,1}.Solution.S6(i) = max((input_Water1(i,1)*states{t,1}.Solution.S6(i-1) + CurrentWater.*states{t,1}.Solution.S6(i))./(input_Water1(i,1)+CurrentWater),0);
        states{t,1}.Solution.S2(i) = max((input_Water1(i,1)*states{t,1}.Solution.S2(i-1) + CurrentWater.*states{t,1}.Solution.S2(i))./(input_Water1(i,1)+CurrentWater),0);
        states{t,1}.Solution.Fe3(i) = max((input_Water1(i,1)*states{t,1}.Solution.Fe3(i-1) + CurrentWater.*states{t,1}.Solution.Fe3(i))./(input_Water1(i,1)+CurrentWater),0);
        states{t,1}.Solution.Fe2(i) = max((input_Water1(i,1)*states{t,1}.Solution.Fe2(i-1) + CurrentWater.*states{t,1}.Solution.Fe2(i))./(input_Water1(i,1)+CurrentWater),0);
        states{t,1}.Solution.Si(i) = max((input_Water1(i,1)*states{t,1}.Solution.Si(i-1) + CurrentWater.*states{t,1}.Solution.Si(i))./(input_Water1(i,1)+CurrentWater),0);
        states{t,1}.Solution.Acetate(i) = max((input_Water1(i,1)*states{t,1}.Solution.Acetate(i-1) + CurrentWater.*states{t,1}.Solution.Acetate(i))./(input_Water1(i,1)+CurrentWater),0);
        states{t,1}.Solution.C4(i) = max((input_Water1(i,1)*states{t,1}.Solution.C4(i-1) + CurrentWater.*states{t,1}.Solution.C4(i))./(input_Water1(i,1)+CurrentWater),0);
        states{t,1}.IonMole.Ca(i) = max((input_Water1(i,1)*states{t,1}.IonMole.Ca(i-1) + CurrentWater.*states{t,1}.IonMole.Ca(i))./(input_Water1(i,1)+CurrentWater),0);
        states{t,1}.IonMole.Mg(i) = max((input_Water1(i,1)*states{t,1}.IonMole.Mg(i-1) + CurrentWater.*states{t,1}.IonMole.Mg(i))./(input_Water1(i,1)+CurrentWater),0);
        states{t,1}.IonMole.Na(i) = max((input_Water1(i,1)*states{t,1}.IonMole.Na(i-1) + CurrentWater.*states{t,1}.IonMole.Na(i))./(input_Water1(i,1)+CurrentWater),0);
        states{t,1}.IonMole.K(i) = max((input_Water1(i,1)*states{t,1}.IonMole.K(i-1) + CurrentWater.*states{t,1}.IonMole.K(i))./(input_Water1(i,1)+CurrentWater),0);
        states{t,1}.IonMole.Cl(i) = max((input_Water1(i,1)*states{t,1}.IonMole.Cl(i-1) + CurrentWater.*states{t,1}.IonMole.Cl(i))./(input_Water1(i,1)+CurrentWater),0);
        states{t,1}.IonMole.S6(i) = max((input_Water1(i,1)*states{t,1}.IonMole.S6(i-1) + CurrentWater.*states{t,1}.IonMole.S6(i))./(input_Water1(i,1)+CurrentWater),0);
        states{t,1}.IonMole.S2(i) = max((input_Water1(i,1)*states{t,1}.IonMole.S2(i-1) + CurrentWater.*states{t,1}.IonMole.S2(i))./(input_Water1(i,1)+CurrentWater),0);
        states{t,1}.IonMole.Fe3(i) = max((input_Water1(i,1)*states{t,1}.IonMole.Fe3(i-1) + CurrentWater.*states{t,1}.IonMole.Fe3(i))./(input_Water1(i,1)+CurrentWater),0);
        states{t,1}.IonMole.Fe2(i) = max((input_Water1(i,1)*states{t,1}.IonMole.Fe2(i-1) + CurrentWater.*states{t,1}.IonMole.Fe2(i))./(input_Water1(i,1)+CurrentWater),0);
        states{t,1}.IonMole.Si(i) = max((input_Water1(i,1)*states{t,1}.IonMole.Si(i-1) + CurrentWater.*states{t,1}.IonMole.Si(i))./(input_Water1(i,1)+CurrentWater),0);
        states{t,1}.IonMole.Acetate(i) = max((input_Water1(i,1)*states{t,1}.IonMole.Acetate(i-1) + CurrentWater.*states{t,1}.IonMole.Acetate(i))./(input_Water1(i,1)+CurrentWater),0);
        states{t,1}.IonMole.C4(i) = max((input_Water1(i,1)*states{t,1}.IonMole.C4(i-1) + CurrentWater.*states{t,1}.IonMole.C4(i))./(input_Water1(i,1)+CurrentWater),0);
    else
        states{t,1}.Solution.Ca(i) = max((input_Water1(i,1)*states{t,1}.Solution.Ca(i-1) + input_Water2(i,1)*states{t,1}.Solution.Ca(i+1) + CurrentWater.*states{t,1}.Solution.Ca(i))./(input_Water1(i,1) + input_Water2(i,1) + CurrentWater),0);
        states{t,1}.Solution.Mg(i) = max((input_Water1(i,1)*states{t,1}.Solution.Mg(i-1) + input_Water2(i,1)*states{t,1}.Solution.Mg(i+1) + CurrentWater.*states{t,1}.Solution.Mg(i))./(input_Water1(i,1) + input_Water2(i,1) + CurrentWater),0); 
        states{t,1}.Solution.Na(i) = max((input_Water1(i,1)*states{t,1}.Solution.Na(i-1) + input_Water2(i,1)*states{t,1}.Solution.Na(i+1) + CurrentWater.*states{t,1}.Solution.Na(i))./(input_Water1(i,1) + input_Water2(i,1) + CurrentWater),0); 
        states{t,1}.Solution.K(i) = max((input_Water1(i,1)*states{t,1}.Solution.K(i-1) + input_Water2(i,1)*states{t,1}.Solution.K(i+1) + CurrentWater.*states{t,1}.Solution.K(i))./(input_Water1(i,1) + input_Water2(i,1) + CurrentWater),0); 
        states{t,1}.Solution.Cl(i) = max((input_Water1(i,1)*states{t,1}.Solution.Cl(i-1) + input_Water2(i,1)*states{t,1}.Solution.Cl(i+1) + CurrentWater.*states{t,1}.Solution.Cl(i))./(input_Water1(i,1) + input_Water2(i,1) + CurrentWater),0); 
        states{t,1}.Solution.S6(i) = max((input_Water1(i,1)*states{t,1}.Solution.S6(i-1) + input_Water2(i,1)*states{t,1}.Solution.S6(i+1) + CurrentWater.*states{t,1}.Solution.S6(i))./(input_Water1(i,1) + input_Water2(i,1) + CurrentWater),0); 
        states{t,1}.Solution.S2(i) = max((input_Water1(i,1)*states{t,1}.Solution.S2(i-1) + input_Water2(i,1)*states{t,1}.Solution.S2(i+1) + CurrentWater.*states{t,1}.Solution.S2(i))./(input_Water1(i,1) + input_Water2(i,1) + CurrentWater),0); 
        states{t,1}.Solution.Fe3(i) = max((input_Water1(i,1)*states{t,1}.Solution.Fe3(i-1) + input_Water2(i,1)*states{t,1}.Solution.Fe3(i+1) + CurrentWater.*states{t,1}.Solution.Fe3(i))./(input_Water1(i,1) + input_Water2(i,1) + CurrentWater),0); 
        states{t,1}.Solution.Fe2(i) = max((input_Water1(i,1)*states{t,1}.Solution.Fe2(i-1) + input_Water2(i,1)*states{t,1}.Solution.Fe2(i+1) + CurrentWater.*states{t,1}.Solution.Fe2(i))./(input_Water1(i,1) + input_Water2(i,1) + CurrentWater),0); 
        states{t,1}.Solution.Si(i) = max((input_Water1(i,1)*states{t,1}.Solution.Si(i-1) + input_Water2(i,1)*states{t,1}.Solution.Si(i+1) + CurrentWater.*states{t,1}.Solution.Si(i))./(input_Water1(i,1) + input_Water2(i,1) + CurrentWater),0); 
        states{t,1}.Solution.Acetate(i) = max((input_Water1(i,1)*states{t,1}.Solution.Acetate(i-1) + input_Water2(i,1)*states{t,1}.Solution.Acetate(i+1) + CurrentWater.*states{t,1}.Solution.Acetate(i))./(input_Water1(i,1) + input_Water2(i,1) + CurrentWater),0); 
        states{t,1}.Solution.C4(i) = max((input_Water1(i,1)*states{t,1}.Solution.C4(i-1) + input_Water2(i,1)*states{t,1}.Solution.C4(i+1) + CurrentWater.*states{t,1}.Solution.C4(i))./(input_Water1(i,1) + input_Water2(i,1) + CurrentWater),0); 
        states{t,1}.IonMole.Ca(i) = max((input_Water1(i,1)*states{t,1}.IonMole.Ca(i-1) + input_Water2(i,1)*states{t,1}.IonMole.Ca(i+1) + CurrentWater.*states{t,1}.IonMole.Ca(i))./(input_Water1(i,1) + input_Water2(i,1) + CurrentWater),0);
        states{t,1}.IonMole.Mg(i) = max((input_Water1(i,1)*states{t,1}.IonMole.Mg(i-1) + input_Water2(i,1)*states{t,1}.IonMole.Mg(i+1) + CurrentWater.*states{t,1}.IonMole.Mg(i))./(input_Water1(i,1) + input_Water2(i,1) + CurrentWater),0); 
        states{t,1}.IonMole.Na(i) = max((input_Water1(i,1)*states{t,1}.IonMole.Na(i-1) + input_Water2(i,1)*states{t,1}.IonMole.Na(i+1) + CurrentWater.*states{t,1}.IonMole.Na(i))./(input_Water1(i,1) + input_Water2(i,1) + CurrentWater),0); 
        states{t,1}.IonMole.K(i) = max((input_Water1(i,1)*states{t,1}.IonMole.K(i-1) + input_Water2(i,1)*states{t,1}.IonMole.K(i+1) + CurrentWater.*states{t,1}.IonMole.K(i))./(input_Water1(i,1) + input_Water2(i,1) + CurrentWater),0); 
        states{t,1}.IonMole.Cl(i) = max((input_Water1(i,1)*states{t,1}.IonMole.Cl(i-1) + input_Water2(i,1)*states{t,1}.IonMole.Cl(i+1) + CurrentWater.*states{t,1}.IonMole.Cl(i))./(input_Water1(i,1) + input_Water2(i,1) + CurrentWater),0); 
        states{t,1}.IonMole.S6(i) = max((input_Water1(i,1)*states{t,1}.IonMole.S6(i-1) + input_Water2(i,1)*states{t,1}.IonMole.S6(i+1) + CurrentWater.*states{t,1}.IonMole.S6(i))./(input_Water1(i,1) + input_Water2(i,1) + CurrentWater),0); 
        states{t,1}.IonMole.S2(i) = max((input_Water1(i,1)*states{t,1}.IonMole.S2(i-1) + input_Water2(i,1)*states{t,1}.IonMole.S2(i+1) + CurrentWater.*states{t,1}.IonMole.S2(i))./(input_Water1(i,1) + input_Water2(i,1) + CurrentWater),0); 
        states{t,1}.IonMole.Fe3(i) = max((input_Water1(i,1)*states{t,1}.IonMole.Fe3(i-1) + input_Water2(i,1)*states{t,1}.IonMole.Fe3(i+1) + CurrentWater.*states{t,1}.IonMole.Fe3(i))./(input_Water1(i,1) + input_Water2(i,1) + CurrentWater),0); 
        states{t,1}.IonMole.Fe2(i) = max((input_Water1(i,1)*states{t,1}.IonMole.Fe2(i-1) + input_Water2(i,1)*states{t,1}.IonMole.Fe2(i+1) + CurrentWater.*states{t,1}.IonMole.Fe2(i))./(input_Water1(i,1) + input_Water2(i,1) + CurrentWater),0); 
        states{t,1}.IonMole.Si(i) = max((input_Water1(i,1)*states{t,1}.IonMole.Si(i-1) + input_Water2(i,1)*states{t,1}.IonMole.Si(i+1) + CurrentWater.*states{t,1}.IonMole.Si(i))./(input_Water1(i,1) + input_Water2(i,1) + CurrentWater),0); 
        states{t,1}.IonMole.Acetate(i) = max((input_Water1(i,1)*states{t,1}.IonMole.Acetate(i-1) + input_Water2(i,1)*states{t,1}.IonMole.Acetate(i+1) + CurrentWater.*states{t,1}.IonMole.Acetate(i))./(input_Water1(i,1) + input_Water2(i,1) + CurrentWater),0); 
        states{t,1}.IonMole.C4(i) = max((input_Water1(i,1)*states{t,1}.IonMole.C4(i-1) + input_Water2(i,1)*states{t,1}.IonMole.C4(i+1) + CurrentWater.*states{t,1}.IonMole.C4(i))./(input_Water1(i,1) + input_Water2(i,1) + CurrentWater),0); 
    end
        
end
            
            
            
                
    
    