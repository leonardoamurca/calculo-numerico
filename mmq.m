% Definir os pontos de entrada da tabela
##xi = [-2.0, -0.5, 1.2, 2.1, 3.5, 5.4];
##yi = [4.4, 5.1, 3.2, 1.6, 0.1, -0.4];

##xi = [0, 0.15, 0.31, 0.5, 0.6, 0.75];
##yi = [1, 1.04, 1.031, 1.117, 1.223, 1.442];

xi = [0, 0.8, 1.8, 3.8, 5.8, 6.8, 7.8, 8.8, 9.8, 10.9, 12.8];
yi= [9.9, 14.3, 17.4, 30.6, 41.2, 51.9, 70.2, 93.1, 119.0, 146.2, 190.7];

% Definir o numero de pontos de entrada
n = size(xi)(2);

% Definir matriz A
Agrau1=[sum(xi.^2), sum(xi); sum(xi), n];

Agrau2=[
  sum(xi.^4), sum(xi.^3), sum(xi.^2);
  sum(xi.^3), sum(xi.^2), sum(xi);
  sum(xi.^2), sum(xi), n;
];

% Definir a matriz b
bgrau1 = [sum(xi.*yi), sum(yi)];
bgrau2 = [sum(xi.^2.*yi),sum(xi.*yi), sum(yi) ];

% Definir funcao para encontrar solucao de sistema linear
function resultado = gauss_seidel(A, b)
  % Definir numero maximo de iteracoes e tolerancia de erro
  maxit = 2000;
  tol = 0.000001;

  [m,n] = size(A);
  C = A;

  for i = 1:n
    C(i,i) = 0;
    w(i) = 0;
  end
  
  x = w';
  % Encontrando a matriz g
  for i = 1:n
    C(i,1:n) = C(i,1:n)/A(i,i);
    g(i,1) = b(i)/A(i,i);
  endfor

  % Definir contador de iteracoes (j)
  j = 0;
  
  % Método de Gauss-Seidel
  while (1)
    xo = x;
    for i=1:n
      x(i,1) = g(i,1)-C(i,:)*x;
      if x(i) ~= 0
        % Criterio de parada
        ea(i) = abs((x(i,1)-xo(i))/x(i,1));
      endif
    endfor
    if max(ea)<=tol || j >= maxit
      resultado = x;
      break
    endif
  endwhile
endfunction

resultadograu1 = gauss_seidel(Agrau1, bgrau1);
agrau1 = resultadograu1(1) ;
bgrau1 = resultadograu1(2);

resultadograu2 = gauss_seidel(Agrau2, bgrau2);
agrau2 = resultadograu2(1) ;
bgrau2 = resultadograu2(2);
cgrau2 = resultadograu2(3);

% Imprime equacao de 1 grau (y = ax + b)
disp(["Equacao 1º grau: " num2str(agrau1) "*x + " num2str(bgrau1)]);

% Imprime equacao de 2 grau (y = ax² + bx + c)
disp(["Equacao 2º grau: " num2str(agrau2) "*x^2 + " num2str(bgrau2) "*x + " num2str(cgrau2)]);

% Definir intervalo e equacoes das funcoes ajustadas
z=0:0.0001:13;
ygrau1 = agrau1.*z + bgrau1;
ygrau2=agrau2.*z.^2 + bgrau2.*z + cgrau2;

% Definir coeficiente de determinacao de cada polinomio
r2grau1 = 1 - (sum( (agrau1.*xi+bgrau1 - yi).^2 ))/ sum( (yi- (1/n)*sum(yi)).^2 );
r2grau2 = 1 - (sum( (agrau2.*xi.^2+bgrau2*xi + cgrau2 - yi).^2 ))/ sum( (yi- (1/n)*sum(yi)).^2 );

% Imprime coeficientes de determinacao
disp(["Coeficiente de Determinacao para equacao de 1º grau: " num2str(r2grau1)]);
disp(["Coeficiente de Determinacao para equacao de 2º grau: " num2str(r2grau2)]);

% Plotar reta ajustada e pontos iniciais
figure
plot(xi,yi,'*',z,ygrau1, z, ygrau2)
title('Ajuste de curvas para polinômio de 1o e 2º grau')
xlabel('x'); ylabel('y')
legend('Dados para ajuste','Reta ajustada','Polinômio de 2o grau ajustado','location','northwest')

