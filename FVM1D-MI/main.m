clear all;
close all;
clc;
nnodos = 41;      %cantidad de nodos
centroides = nnodos-1; %cantidad de centroides
intervalo = [0,1]; % intervalor de X

dx = (intervalo(2)-intervalo(1))/(centroides); % paso y ANCHO DE LA CELDA

c=0;      %COEFICIENTE REACTIVO
k=0.1;      %COEFICIENTE DIFUSIVO
s=1;      %ESPESOR DE LA PLACA
rho= 0;
cp= 0;

##condiciones para el termino advectivo
model.u=5;

%condiciones de borde 1 dirichlet; 2 neumann ; 3 robin
model.bc = [ 1 , 0, -1;
             1 , 0,  -1];
             
             
             
##nodos = intervalo(1):dx:intervalo(2); %nodos

## r -> factor de comprension
## r < 1 -> malla refinada hacia el extremo derecho
## r > 1 -> malla refinada hacia el extremo izquierdo
nodos = genMalla(0,1,0.8,nnodos); %nodos con estrechamiento

for i=1:length(nodos)-1
  xx(i) = nodos(i)+((nodos(i+1)-nodos(i))/2); 
end

model.nodos = nodos;
##model.dc = intervalo(1)+dx/2 : dx: intervalo(2)-dx/2; %posicion de centroides
##
##xx =model.dc ; %vector con x de centroides

cantidad_nodos = length(nodos); 
##GG =xx.^3;                  %FUENTE VARIABLE
##model.G = GG;

G= 10;                         %FUENTE CONSTANTE
model.G= ones(cantidad_nodos-1,1)*G;

model.c= ones(cantidad_nodos-1,1)*c;
model.k= ones(cantidad_nodos-1,1)*k;
model.dx = dx;
model.dy = 1;
model.t = 0;
model.s = ones(cantidad_nodos-1,1)*s;  %ESPESOR DE  LA PLACA

%---------------Parametros para elegir el modelo -------------------------------
## que modelo utilizar
## ts = -1 -> estacionario
## ts = 0 -> expl?¿½cito
## ts = 1 -> impl?¿½cito
model.ts = -1;

% Parametros para esquemas temporales
model.rho=rho;
model.cp=cp;
model.maxit = 2000;
model.tol = 0.00001;
model.t0=0;

model.dt = 0.01;

%---------------------------Metodo FVM -----------------------------------------

##CD
model.schemeu=1;
[cells,K,F,Phi,model] = FVM_1D(nodos, centroides, model); 

##Upwind
model.schemeu=2;
[cellsU,KU,FU,PhiU,model] = FVM_1D(nodos, centroides, model); 

##anal?tica
T=ones(length(xx),1);
for i=1:centroides
    xx(i)=cells(i).dc;
    T(i)= (model.G(i)/model.u)*(xx(i)-(1-exp(model.u*xx(i)/model.k(i)))/(1-exp(model.u/model.k(i))));
end 

%-------------Comparativa T - Phi en los centros de celdas----------------------
%-----------------------------Figuras-------------------------------------------
error = norm(T-Phi,2)/norm(T,2);
errorU = norm(T-PhiU,2)/norm(T,2);

Phi=[0; Phi; 0];
PhiU=[0; PhiU; 0];
xx=[0 xx 1];
T=[0; T; 0];

##CD
plot(xx,Phi,'b--x-')
hold on

##Upwind
plot(xx,PhiU,'g-x-')
hold on

Pe = model.u*dx/k;
##anal?tica
plot(xx,T,'r-o-')
title('Caso 3')
text(0.2,1.8,strcat('G=10[W/m]; v=5[m/s]; k=',num2str(k),'[m^2/s]'))
text(0.2,1.5,strcat("Pe_x=",num2str(Pe) ))
legend('PHI CD- FVM', 'PHI UP- FVM', 'Analitica', "location", "northeastoutside")

##title(strcat(strcat("Error CD=", num2str(error)),strcat("-Error UP=", num2str(errorU))));
text(0.1,1.2,strcat("Error CD=", num2str(error)));
text(0.1,1,strcat("Error UP=", num2str(errorU)));


hold off