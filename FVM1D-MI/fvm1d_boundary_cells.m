function [K, F, cells] = fvm1d_boundary_cells(K,F,cells, model)

%---------------Borde izquierdo-------------------------------------------------
    
    if (model.bc(1,1) == 1) %DIRICHLET
##      be=((cells(1).de)/2)/(cells(1).de);
##      bw=((cells(1).dw))/(cells(1).dw);
##      
##      Rw = cells(1).kw*cells(1).Aw/cells(1).dw;
##      Re = cells(1).ke*cells(1).Ae/cells(1).de;
##      K(1,1) = 1;
##      K(1,2) = 0;
##      
##      phi   = model.bc(1,2); % vlaor de dirichlet en la izquierda
##      F(1)  = phi;


      be=((cells(1).de)/2)/(cells(1).de);
      bw=1;
      Rw = cells(1).kw*cells(1).Aw/cells(1).dw; %si es borde lo contempla
      Re = cells(1).ke*cells(1).Ae/cells(1).de;
      K(1,1) = cells(1).c*cells(1).v + Re + Rw + model.u*(be);
      K(1,2) = -Re + model.u*be;
      phi   = model.bc(1,2); % vlaor de dirichlet en la izquierda
      F(1)  = cells(1).v*model.G(1) + Rw*phi + model.u*phi;
      
      if (model.schemeu==2) %upwind  
        Rw = cells(1).kw*cells(1).Aw/cells(1).dw;
        Re = cells(1).ke*cells(1).Ae/cells(1).de;
        K(1,1) = cells(1).c*cells(1).v + Re + Rw;
        K(1,2) = -Re;
        
        phi   = model.bc(1,2); % vlaor de dirichlet en la izquierda
        F(1)  = cells(1).v*model.G(1)+ Rw*phi;
        
        K(1,1) += model.u;
        F(1)  += model.u*phi;
      end
      
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
      
##      be=((cells(N).de))/(cells(N).de);
##      bw=((cells(N).dw)/2)/(cells(N).dw);
##      
##      Rw = cells(N).kw*cells(N).Aw/cells(N).dw;
##      Re = cells(N).ke*cells(N).Ae/cells(N).de;
##      K(N,N) = 1;
##      K(N,N-1) = 0;
##      phi   = model.bc(2,2); % vlaor de dirichlet en la derecha
##      F(N)  = phi;


      bw=((cells(N).dw)/2)/(cells(N).dw);
      Rw = cells(N).kw*cells(N).Aw/cells(N).dw;
      Re = cells(N).ke*cells(N).Ae/cells(N).de;
      K(N,N) = cells(N).c*cells(N).v + Re + Rw - model.u*(bw);
      K(N,N-1) = -Rw -model.u*bw;
      phi   = model.bc(2,2); % vlaor de dirichlet en la derecha
      F(N)  = cells(N).v*model.G(N)+ Re*phi - model.u*phi;      
      
      if (model.schemeu==2)   %upwind    
        Rw = cells(N).kw*cells(N).Aw/cells(N).dw;
        Re = cells(N).ke*cells(N).Ae/cells(N).de;
        K(N,N) = cells(N).c*cells(N).v + Re + Rw;
        K(N,N-1) = -Rw;
        phi   = model.bc(2,2); % vlaor de dirichlet en la derecha
        F(N)  = cells(N).v*model.G(N)+ Re*phi;
        
        K(N,N) += model.u;
        K(N,N-1) -= model.u;
      end
     
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