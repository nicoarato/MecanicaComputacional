function [K, F, cells] = fvm1d_inner_cells(K,F,cells,model)
G=model.G;

##coeficientes que se agregan al tener termino advectivo
##velocidad
ue=model.u;
uw=model.u;

N = length(cells);
  for i = 2 : N-1
    ##betas
    be=((cells(i).dx)/2)/(cells(i+1).dc-cells(i).dc);
    bw=((cells(i-1).dx)/2)/(cells(i).dc-cells(i-1).dc);
    
    Re = cells(i).ke*cells(i).Ae/cells(i).de;
    Rw = cells(i).kw*cells(i).Aw/cells(i).dw;
    Vp = cells(i).v;
    
    ##    sin adveccion
    if (model.u == 0)
      K(i,i)   =  Re + Rw + cells(i).c*Vp;
      K(i,i-1) = -Rw;
      K(i,i+1) = -Re;    

    ##  central difference
    elseif (model.u ~= 0 && model.schemeu == 1) 
      K(i,i+1) =  model.u*be - Re;
      K(i,i)   =  Re + Rw + model.u*(bw - be) + cells(i).c*Vp;
      K(i,i-1) = -Rw - model.u*bw;
      
      
    ##      upwind  
    elseif (model.u ~= 0 && model.schemeu == 2)
      if (model.u > 0) 
        K(i,i+1) = -Re;
        K(i,i)   =  Re + Rw + model.u + cells(i).c*Vp;
        K(i,i-1) = -Rw - model.u;
        
      elseif (model.u < 0) 
        K(i,i+1) = -Re + model.u;
        K(i,i)   =  Re + Rw - model.u + cells(i).c*Vp;
        K(i,i-1) = -Rw;
      end
    end
    
    F(i) = G(i)*Vp;
    
  end
end