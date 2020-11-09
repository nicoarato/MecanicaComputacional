function [cells,K,F,PHI_vec,model] = fvm1d_explicito(K,F,cells,model,dt, centroides)
   
    A = dt/(model.rho*model.cp);
    n=length(cells);
    
    x = model.dc; %coordenadas de los centroides
    T= model.TT; %solucion analitica evaluada en centroides
    
    % Initialize phi solution vector.
    PHI = model.phi_n;
    PHI_n = model.phi_n;
    PHI_vec = PHI;
    err = 1;
##    Q_vec = zeros(model.ncells,2);

    % Identity matrix
    I = eye(centroides,centroides);

    for n = 1 : model.maxit
        PHI = A*F + (I - A*K)*PHI_n;
        
        plot(x,T,'r--x-');
        hold on
        plot(x,PHI,'b--o-'); 
        legend('Solucion analitica', 'Forward Euler')
        title(sprintf('nit: %d - error: %e',n,err));
        pause(0.000000000000001);
        hold off  
       
        err = norm(PHI-PHI_n,2)/norm(PHI,2);

        % actualizo phi(n+1) será phi(n) para el siguiente paso
        PHI_n = PHI;
        PHI_vec = [PHI_vec PHI];
        
##        [Q] = fvm2d_flux(PHI,cells,neighb);
##        Q_vec = [Q_vec, Q];
        
        if err < model.tol
          if (err == 0)
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
    
    disp('Corte por alcanzar cantidad de iteraciones');
    fprintf('El error es de: %d\n',err);
end
