function [icone]= generar_icone(xnode,x ,y)
  n= length(x);
  m= length(y);
  icone = [];
  nxnode=length(xnode);
  for i=1 : nxnode-n
      
      if mod(i,n)~=0
        icone = [icone; i n+i i+n+1 i+1]
      
      end
      
  
  end
  
  
end