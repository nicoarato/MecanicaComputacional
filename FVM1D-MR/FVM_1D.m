function [cells,K,F,Phi,model] = FVM_1D(nodos,centroides, model)

  %Resuelve sin adveccion -> sin velocidad
  % dT/dt = D2T/Dx2 -cT + G
  
  [cells] = fvm1d_inicialize1D(nodos, model);

  %calcular celdas interiores
  K = zeros(centroides,centroides);
  F= zeros(centroides,1);
  [K, F, cells] = fvm1d_inner_cells(K,F,cells,model);


  %calcular celdas en los bordes
  [K, F, cells] = fvm1d_boundary_cells(K,F,cells,model, centroides);
  
  if(model.ts==-1) %ESTACIONARIO
    Phi = K\F;
  end
  
  if( model.ts==0) %EXPLICITO
    dt= model.dt;
    model.phi_n = ones(centroides,1)*model.t0;
    [cells,K,F,Phi,model] = fvm1d_explicito(K,F,cells,model,dt, centroides);
    model.PHI_explicito= Phi;
    Phi = Phi(:,end);
  end

  if (model.ts == 1) %IMPLICITO
    dt= model.dt;
    model.phi_n = ones(centroides,1)*model.t0;
    [cells,K,F,Phi,model] = fvm1d_implicito(K,F,cells,model,dt, centroides);
    model.PHI_implicito= Phi;
    Phi = Phi(:,end);
  end
  

end