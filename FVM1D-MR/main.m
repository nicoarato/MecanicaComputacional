clear all;
close all;
clc;
nnodos = 10;      %cantidad de nodos
centroides = nnodos-1; %cantidad de centroides
intervalo = [0,1]; % intervalor de X

dx = (intervalo(2)-intervalo(1))/(centroides); % paso y ANCHO DE LA CELDA


c=-2;      %COEFICIENTE REACTIVO
k=2;      %COEFICIENTE DIFUSIVO
s=1;      %ESPESOR DE LA PLACA
rho= 1;
cp= 1;

%condiciones de borde 1 dirichlet; 2 neumann ; 3 robin
model.bc = [ 1 , 50, -1;
             2 , 5, -1];

nodos = intervalo(1):dx:intervalo(2); %nodos
model.nodos = nodos;
model.dc = intervalo(1)+dx/2 : dx: intervalo(2)-dx/2; %posicion de centroides


xx =model.dc ; %vector con x de centroides

cantidad_nodos = length(nodos); 
##GG =xx.^3;                  %FUENTE VARIABLE
##model.G = GG;

G= 0;                         %FUENTE CONSTANTE
model.G= ones(cantidad_nodos-1,1)*G;

model.c= ones(cantidad_nodos-1,1)*c;
model.k= ones(cantidad_nodos-1,1)*k;
model.dx = dx;
model.dy = 1;
model.t = 0;
model.s = ones(cantidad_nodos-1,1)*s;  %ESPESOR DE  LA PLACA

%------------------- Solucion analitica y PHI----------------------------------------
x = nodos;
##---------Ejercicio 1a---------------
##T= -25*x.^2+65.*x + 10; 
##Terr= -25*xx.^2+65.*xx + 10; 

##---------Ejercicio 1b----------------
##T= (100*exp(-x).*(exp(2*x).+exp(4)))./(1.+ exp(4));
##Terr=(100*exp(-xx).*(exp(2*xx).+exp(4)))./(1.+ exp(4)); 

##---------Ejercicio 1c-----------------
##T= -25*x.^4+300*x.^3-1350*x.^2+1906*x+2345;
##T=T/3;

##Terr= -25*xx.^4+300*xx.^3-1350*xx.^2+1906*xx+2345; 
##Terr=Terr/3;

##---------Ejercicio 1d-----------------
##T= -36.6897*exp(-x)-3.3103*exp(x)+50;
##Terr= -36.6897*exp(-xx)-3.3103*exp(xx)+50;

##---------Ejercicio 1e-----------------
##T= -(x.^5)*(1/40)+(1225/3)*x-4600/3;
##Terr= -(xx.^5)*(1/40)+(1225/3)*xx-4600/3;

##Ejercicio 1 f
##T = (-5/4).*exp(-x-1).*(exp(x)-1).*(11*exp(x)+(11-30*exp(1)));
##Terr = (-5/4).*exp(-xx-1).*(exp(xx)-1).*(11*exp(xx)+(11-30*exp(1)));

##Ejercicio 1 g
T = 73.2433*sin(x)+ 50*cos(x);
Terr  = 73.2433*sin(xx)+ 50*cos(xx);

model.TT=Terr;
%---------------Parametros para elegir el modelo -------------------------------
## que modelo utilizar
## ts = -1 -> estacionario
## ts = 0 -> explícito
## ts = 1 -> implícito
## ts = 1/2 -> semi-implícito
model.ts = 1;

% Parametros para esquemas temporales
model.rho=rho;
model.cp=cp;
model.maxit = 2000;
model.tol = 0.00001;
model.t0=0;

model.dt = 0.01;

%---------------------------Metodo FVM -----------------------------------------


  [cells,K,F,Phi,model] = FVM_1D(nodos, centroides, model); 


%-------------Comparativa T - Phi en los centros de celdas----------------------
%-----------------------------Figuras-------------------------------------------
figure(1)
subplot(131), plot(x,T) , xlabel('x'), ylabel('T'), title('Analitica: Temperatura');
subplot(132), plot(xx,Phi) , xlabel('x'), ylabel('PHI'), title('Temperatura en centros de celda') ;
subplot(133), 
hold on 
plot(xx,Terr,'r-o-') , xlabel('x'), ylabel('PHI-T');
hold on
plot(xx,Phi,'b-x-') , xlabel('x');
legend('Analitica', 'PHI- FVM')

##subplot(133), plot(xx,Terror) , xlabel('x'), ylabel('E'), title('ERROR') ;
error = norm(Terr-Phi,2)/norm(Terr,2)

figure(2)
hold on
plot(xx,Terr,'r-o-')
hold on
plot(xx,Phi,'b--x-')
legend('Analitica', 'PHI- FVM')
