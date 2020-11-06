function [K, F, cells] = fvm1d_boundary_cells(K,F,cells, model)

%---------------Borde izquierdo-------------------------------------------------
    if (model.bc(1,1) == 1) %DIRICHLET
      Rw = cells(1).kw*cells(1).Aw/cells(1).dw;
      Re = cells(1).ke*cells(1).Ae/cells(1).de;
      K(1,1) = cells(1).c*cells(1).v + Re + Rw;
      K(1,2) = -Re;
      
      phi   = model.bc(1,2); % vlaor de dirichlet en la izquierda
      F(1)  = cells(1).v*model.G(1)+ Rw*phi;
      
    elseif (model.bc(1,1) == 2) %NEUMANN
        
      Rw = cells(1).kw*cells(1).Aw/cells(1).dw;
      Re = cells(1).ke*cells(1).Ae/cells(1).de;
      K(1,1) = cells(1).c*cells(1).v + Re;
      K(1,2) = -Re;
      
      q   = model.bc(1,2); % vlaor de flujo en la izquierda
      F(1)  = cells(1).v*model.G(1)- q*cells(1).Aw;
      
    else  %(model.bc(1,1) == 'r') %ROBIN  en la izquierda
      
      h=model.bc(1,2);
      Tinf= model.bc(1,3);
         
      Rw = cells(1).kw*cells(1).Aw/cells(1).dw;
      Re = cells(1).ke*cells(1).Ae/cells(1).de;
      
      aw = h*Tinf/(cells(1).k+(h*cells(1).dw));
      bw = -h/(cells(1).k+(h*cells(1).dw));
      
      K(1,1) = cells(1).c*cells(1).v + Re - cells(1).kw*cells(1).Aw*bw;
      K(1,2) = -Re;
      
      r   = cells(1).kw*cells(1).Aw*aw; % robin en la izquierda
      F(1)  = cells(1).v*model.G(1) + r;
      
      
      
    end   %end if

%----------------borde derecho--------------------------------------------------
    N= length(cells);
    if (model.bc(2,1) == 1) %DIRICHLET
      Rw = cells(N).kw*cells(N).Aw/cells(N).dw;
      Re = cells(N).ke*cells(N).Ae/cells(N).de;
      K(N,N) = cells(N).c*cells(N).v + Re + Rw;
      K(N,N-1) = -Rw;
      phi   = model.bc(2,2); % vlaor de dirichlet en la derecha
      F(N)  = cells(N).v*model.G(N)+ Re*phi;
     
    elseif (model.bc(2,1) == 2) %NEUMANN borde derecho
      
      Rw = cells(N).kw*cells(N).Aw/cells(N).dw;

      K(N,N) = cells(N).c*cells(N).v + Rw;
      K(N,N-1) = -Rw;
      
      q   = model.bc(2,2); % vlaor de flujo en la derecha
      F(N)  = cells(N).v*model.G(N)- q*cells(N).Ae;
    
    else  %(model.bc(2,1) == 'r')     ROBIN borde derecho
      
      h=model.bc(2,2);
      Tinf= model.bc(2,3);
         
      Rw = cells(N).kw*cells(N).Aw/cells(N).dw;
      Re = cells(N).ke*cells(N).Ae/cells(N).de;
      
      ae = h*Tinf/(cells(N).k+(h*cells(N).de));
      be = -h/(cells(N).k+(h*cells(N).de));
      
      K(N,N) = cells(N).c*cells(N).v + Rw - cells(N).ke*cells(N).Ae*be;
      K(N,N-1) = -Rw;
      
      r   = cells(N).ke*cells(N).Ae*ae; % robin en la derecha
      F(N)  = cells(N).v*model.G(N) + r;
    end
    
end