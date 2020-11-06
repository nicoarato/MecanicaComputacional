function [K, F, cells] = fvm1d_inner_cells(K,F,cells,model)
G=model.G;
N = length(cells);
  for i = 2 : N-1
    
    Re = cells(i).ke*cells(i).Ae/cells(i).de;
    Rw = cells(i).kw*cells(i).Aw/cells(i).dw;
    Vp = cells(i).v;
    
    K(i,i)   =  Re + Rw + cells(i).c*Vp;
    K(i,i-1) = -Rw;
    K(i,i+1) = -Re;
    
    F(i) = G(i)*Vp;
    
  end
end