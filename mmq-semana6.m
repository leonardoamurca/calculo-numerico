% Utilizado no exercicio 1
##xi = [log(0.5), log(0.75), log(1), log(1.5), log(2), log(3)];
##yi= [-2.8, -0.6, 1, 3.2, 4.8, 7];

% Utilizado no exericio 2
xi = [-8, -6, -4, -2, 0, 2, 4];
yi = [log(30), log(10), log(9), log(6), log(5), log(4), log(4)];

% Definir o numero de pontos de entrada
n = size(xi)(2);

% Definir matriz A
Agrau1=[sum(xi.^2), sum(xi); sum(xi), n];

% Definir a matriz b
bgrau1 = [sum(xi.*yi), sum(yi)];

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

% Imprime equacao de 1 grau (y = a*ln(x) + b) (exercicio 1)
##disp(["Equacao 1º grau: " num2str(agrau1) "*ln(x) + " num2str(bgrau1)]);

% Imprime equacao de 1 grau (y = b*e^(a*x)) (exericio 2)
disp(["Equacao 1º grau: " num2str(bgrau1) "*e^" num2str(agrau1) "x" ]);

% Definir intervalo e equacoes das funcoes ajustadas (exercicio 1)
##z=0:0.0001:3.5;

% Definir intervalo e equacoes das funcoes ajustadas (exercicio 2)
z=-10:0.0001:5;

% Para exericio 1
##ygrau1 = agrau1.*log(z) + bgrau1;

% Para exericio 2
ygrau1 = bgrau1*e.^(agrau1.*z);

% Definir coeficiente de determinacao de cada polinomio (exercicio 1)
##r2grau1 = 1 - (sum( (agrau1.*xi+bgrau1 - yi).^2 ))/ sum( (yi- (1/n)*sum(yi)).^2 );

% Definir coeficiente de determinacao de cada polinomio (exercicio 2)
r2=((n*sum(xi.*yi)-sum(xi)*sum(yi))/sqrt(n*sum(xi.^2)-sum(xi.^2))/sqrt(n*sum(yi.^2)-sum(yi.^2)))^2

% Imprime coeficiente de determinacao
disp(["Coeficiente de Determinacao para equacao de 1º grau: " num2str(r2grau1)]);

% Definir pontos originais para plotar os pontos da tabela inicial (exercicio 1)
# xi = [0.5, 0.75, 1, 1.5, 2, 3];

% Definir pontos originais para plotar os pontos da tabela inicial (exercicio 2)
yi = [30, 10, 9, 6, 5, 4, 4];
% Plotar reta ajustada e pontos iniciais
figure
plot(xi,yi,'*',z,ygrau1)
title('Ajuste de curvas para polinômio de 1o (y = b*e^(a*x))')
xlabel('x'); ylabel('y')
legend('Dados para ajuste','Reta ajustada','location','northwest')

