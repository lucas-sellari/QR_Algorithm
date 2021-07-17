import numpy as np
import math  as m
import sys

# Respostas analiticas para os autovalores e autovetores
def gabaritoA(n):
    vet = np.eye(n)
    val = np.zeros(n)

    for i in range(n):
        for j in range(n):
            vet[j][n - 1 - i] = m.sin((m.pi / (n + 1)) * (j + 1) * (i + 1))

    for l in range(n):
        val[n - 1 - l] = 2 * (1 - m.cos((m.pi / (n + 1)) * (l + 1)))

    return vet, val


# Funcao sgn
def sgn(d):
    if d >= 0:
        result =  1
    else:
        result = -1
    return result


# Funcao para a heuristica de wilkinson
def wilkinson(a_n1, a_n, b_n1):
    d_k = (a_n1 - a_n)/2
    if (d_k ** 2) >= (b_n1 ** 2):
        mu_k = a_n + d_k - sgn(d_k) * m.sqrt((d_k ** 2) - (b_n1 ** 2))
    else:
        mu_k = 0.0

    return mu_k


def QR_factorization(a, n, mu_k, v):
    c = np.zeros(n - 1)  # vetor de 0's de tamanho n-1
    s = np.zeros(n - 1)

    QkTo = np.eye(n)
    QkT = np.eye(n)

    # Vamos encontrar a matriz Qk e QkT:
    for k in range(0, n - 1, 1):
        Qk = np.eye(n)

        # Definindo os ck e sk
        ck = a[k][k] / m.sqrt(a[k][k] ** 2 + a[k + 1][k] ** 2)
        sk = -a[k + 1][k] / m.sqrt(a[k][k] ** 2 + a[k + 1][k] ** 2)

        # Definindo Qk
        for i in range(n):
            for j in range(n):
                if i == k and j == k:
                    Qk[i][j] = ck
                    Qk[i + 1][j + 1] = ck
                    Qk[i][j + 1] = -sk
                    Qk[i + 1][j] = sk

        QkT = Qk.transpose()
        QkT = np.dot(QkTo, QkT)
        QkTo = QkT

        a = np.dot(Qk, a)  # Ao final, esta sera a matriz R

    # Vamos determinar a matriz A^(k + 1):
    a = np.dot(a, QkT)

    # Ajuste da diagonal principal de A^(k + 1)
    for i in range(n):
        a[i][i] = a[i][i] + mu_k

    # Agora vamos calcular a matriz V^(k + 1):
    v = np.dot(v, QkT)

    return a, v


def QR_algorithm(a, n, eps):
    k = 0
    v = np.eye(n)     # matriz identidade de tamanho nxn

    #  "1" para QR com deslocamento espectral ou "0" para QR sem deslocamento espectral
    desloc = sys.argv[3]

    for l in range(n - 1, 0, -1):

        while m.fabs(a[l - 1][l - 2]) >= eps:

            if k > 0 and desloc == "1":
                mu_k = wilkinson(a[n - 2][n - 2], a[n - 1][n - 1], a[n - 1][n - 2])

            else:
                mu_k = 0.0

            # atualizando a diagonal principal, linha 6 do pseudocodigo
            for i in range(n):
                a[i][i] = a[i][i] - mu_k

            # fatoracao QR
            a, v = QR_factorization(a, n, mu_k, v)

            k = k + 1

    print("\nNúmero de passos computados: ", k)
    return a, v


if __name__ == '__main__':
    # n desejado
    n = int(sys.argv[1])

    # epsilon desejado (utilizar ponto no lugar da vírgula)
    eps = float(sys.argv[2])

    a = np.zeros((n, n))

    # Definindo a matriz A de acordo com o exercicio A do enunciado
    for i in range(1, n - 1, 1):
        a[i - 1][i] = -1.0
        a[i][i]     =  2.0
        a[i + 1][i] = -1.0

    a[0][0] = 2.0
    a[1][0] = -1.0
    a[n - 2][n - 1] = -1.0
    a[n - 1][n - 1] =  2.0

    a, v = QR_algorithm(a, n, eps) # chamando a funcao que faz o calculo numerico

    vet, val = gabaritoA(n)        # chamando a funcao gabarito

    print("\nGabarito autovalores:")
    for i in range(n):
        print("a[{}]: {}".format(i, val[i]))

    print("\nAutovalores obtidos:")
    for i in range(n):
        print("a[{}]: {}".format(i, np.round(a[i][i], 2)))

    print("\nGabarito autovetores: \n", vet)
    print("\nAutovetores obtidos(arredondado): \n", np.round(v, 2))

    # Calculo dos erros relativos dos autovalores, em percentual
    compval = np.zeros(n)

    for i in range(n):
        compval[i] = val[i] - a[i][i]
        compval[i] = (compval[i] / val[i]) * 100

    print("\nErro relativo dos autovalores (em %): \n", compval)