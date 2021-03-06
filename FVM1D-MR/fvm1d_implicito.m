function [cells,K,F,PHI_vec,model] = fvm1d_implicito(K,F,cells,model,dt, centroides)
  
  A = dt/(model.cp*model.rho);
  I = eye(centroides,centroides);
  
  x = model.dc; %coordenadas de los centroides
  T= model.TT;  %%solucion analitica evaluada en centroides
  
  PHI = model.phi_n;
  PHI_n = model.phi_n;
  PHI_vec = PHI;
  
  for n=2 : model.maxit  
    PHI = (I + A*K) \ (PHI_n + A*F);
    
##    [Q] = fdm2d_flux(PHI,neighb,xnode,model.k);
##     Q_vec = [Q_vec, Q];
      if( mod(n,10) == 0)  
        plot(x,T,'r--x-');
        hold on
        plot(x,PHI,'b--o-'); 
        legend('Solucion analitica', 'Backward Euler')
        title(n);
        pause(0.000000000001);
        hold off  
      end
     
      error = norm(PHI-PHI_n,2)/norm(PHI,2)
      
      PHI_n = PHI;
      PHI_vec = [PHI_vec PHI];
      
    
      if error < model.tol
          if (error == 0)
            disp ('La solucion ha llegado a un estado estacionario');
            fprintf('Los pasos realizados fueron: %d\n',n);
            fprintf('El tiempo transcurrido fue de: %f segundos\n',n*dt);
           else
            disp ('Se alcanz� la tolerancia especificada.');
            fprintf('Los pasos realizados fueron: %d\n',n);
            fprintf('El tiempo transcurrido fue de: %f segundos\n',n*dt);
           end
            return;
        end 
        
  end

##    PHI=PHI_vec;
##    Q=Q_vec;
  disp('Corte por alcanzar cantidad de iteraciones');
  
end