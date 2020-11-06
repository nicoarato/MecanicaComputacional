function [xnode] = generar_malla(x, y)
  nx = length(x);
  ny = length(y);
  xnode = [];
  for i = 1 : nx
    for j = 1 : ny
      xnode= [xnode; x(i) y(j)];
    end
  end
  
  
end