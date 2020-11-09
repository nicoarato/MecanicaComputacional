function [cells]= fvm1d_inicialize1D(nodos,model)


## inicializa las celdas con los valores 
## de distancia al centoide este
## dw distancia al centoide oeste
## ke valor de k al este
## kw valor de k al oeste
## t  valor de espesor
## dx valor de ancho de celda
## dy altura de la celda
## dc centro de centroide
## c el valor de constante c en la celda
##  v volumen de la celda
cells = struct('de', {}, 'dw', {},...
                   'Ae', {}, 'Aw', {},...
                   'ke', {},...
                   'kw', {},...
                   'dy', {}, 't', {},...
                   'dc', {},'dx', {}, 'c', {},'v',{});
N = length(nodos)-1; %cantidad de celdas
K = zeros(N);
F= zeros(N,1);   
G = model.G;   

  for i = 1 : N
    dx= (nodos(i+1)-nodos(i));%ancho de la celda
    cells(i).dx = dx;             %ancho de la celda
    cells(i).dc = nodos(i) + dx/2; %centro de la celda
    
    if (i == 1)
      cells(i).dw = dx/2;         %distancia al borde izquerdo
   
    else
      cells(i).dw = abs(cells(i-1).dc - cells(i).dc); %distancia entre centroides celdas internas
    end    %endif
    
    if (i == N)
      cells(i).de =  dx/2;     %distancia al borde derecho
    else
      cells(i).de = dx/2 +(nodos(i+2) - nodos(i+1))/2; %distancia al centroide derecho
    end  %endif
    
    cells(i).dy = model.dy;    %altura de la celda
    cells(i).t = model.s(i);   %espesor de la celda       
    cells(i).Aw = cells(i).dy*cells(i).t; %area a la izquierda
    cells(i).Ae = cells(i).dy*cells(i).t; %area a la derecha
    cells(i).c = model.c(i);        %valor de c de la ecuacion
    cells(i).v = cells(i).dy*cells(i).dx*cells(i).t; %volumen 
    cells(i).k = model.k(i);    %k componente de difusion en el centro de la celda
    
    if(i == 1)    %difusion en las caras west
      cells(i).kw = cells(i).k;
    else
      cells(i).kw = cells(i-1).k;
    end %end if
    
    if(i == N)      %difusion en las caras este
      cells(i).ke = cells(i).k;
    else 
      cells(i).ke = model.k(i+1);
    end %end if
    
    
  end %endfor

        
end                  