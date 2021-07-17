
 Códigos desenvolvidos utilizando o Python na versão 3.9.1, a biblioteca NumPy na versão 1.20.1 e a biblioteca matplotlib na versão 3.4.2.                 
 Os arquivos podem ser executados diretamente no terminal de sua preferência.                                                                            
 Nenhum arquivo adicional é requerido.
 
 ========================================================================================
 
 Instruções para o arquivo A.py:
 
 1. Os argumentos devem ser fornecidos junto com o comando para executar o arquivo. Temos a seguinte forma para a linha de comando:
 
 		                                             python A.py  n  eps  ops
 
    Em que "n" é o parâmetro para o tamanho da matriz A, que será n x n. "eps" é o epsilon desejado (utilizar ponto ao invés de vírgula).
    "ops" indica se o algoritmo é com ou sem deslocamento espectral. Use "1" para indicar que deseja que o algoritmo seja executado com
    deslocamento e "0" para o algoritmo sem deslocamento.
 
 2. O script volta algumas informações, como:
    - Número de passos computados;
    - Os valores-gabarito para os autovalores e os valores de autovalores obtidos pelo algoritmo;
    - A matriz-gabarito dos autovetores e a matriz de autovalores obtida pela execução do algoritmo;
    - Array de erros relativos dos autovalores.
 
 ========================================================================================
 
 Instruções para o arquivo B.py:
 
 1. Os argumentos devem ser fornecidos junto com o comando para executar o arquivo. Temos a seguinte forma para a linha de comando:
 
 		                                             python B.py  ci  ops
 
    Em que "ci" indica quais são as condições iniciais a serem consideradas. Seus valores podem ser 1, 2 ou 3, de acordo com a lista:
 			(1) X(0) = (-2, -3, -1, -3, -1)
 			(2) X(0) = (1,  10, -4,  3, -2)
 			(3) X(0) correspondente ao modo de maior frequência
 
    "ops" indica se o algoritmo é com ou sem deslocamento espectral. Use "1" para indicar que deseja que o algoritmo seja executado com
    deslocamento e "0" para o algoritmo sem deslocamento.
 
 2. O script volta os gráficos de deslocamento (eixo y) versus tempo (eixo x) dos 5 corpos, um atrás do outro.
 
 
 ========================================================================================
 
 Instruções para o arquivo C.py:
 
 1. Os argumentos devem ser fornecidos junto com o comando para executar o arquivo. Temos a seguinte forma para a linha de comando:
 
 		                                             python C.py  ci  ops
 
    Em que "ci" indica quais são as condições iniciais a serem consideradas. Seus valores podem ser 1, 2 ou 3, de acordo com a lista:
 			(1) X(0) = (-2, -3, -1, -3, -1, -2, -3, -1, -3, -1)
 			(2) X(0) = (1,  10, -4,  3, -2, 1,  10, -4,  3, -2)
 			(3) X(0) correspondente ao modo de maior frequência
 
    "ops" indica se o algoritmo é com ou sem deslocamento espectral. Use "1" para indicar que deseja que o algoritmo seja executado com
    deslocamento e "0" para o algoritmo sem deslocamento.
 
 2. O script volta os gráficos de deslocamento (eixo y) versus tempo (eixo x) dos 10 corpos, um atrás do outro.
 
 
 ========================================================================================
