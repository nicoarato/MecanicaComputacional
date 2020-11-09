function [xnode] = genMalla(x0,xL,r,n)
## x0 -> Extremo izquierdo del dominio
## xL -> Extremo derecho del dominio
## r -> factor de comprension
## r < 1 -> malla refinada hacia el extremo derecho
## r > 1 -> malla refinada hacia el extremo izquierdo
## n -> cantidad de nodos de la malla
## xnode -> malla generada   
  L = xL-x0;
  if r==1  
    dx = L/(n-1);
    xnode = x0:dx:xL;
  else  
    S = (1-r^n)/(1-r);
    dx = L/S;
    xnode = [x0];
    x = x0;
    for i=0:n-1
      dxi = r^i*dx;
      x = x+dxi;
      xnode = [xnode x];
    endfor
  endif  
endfunction
