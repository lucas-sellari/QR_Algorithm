import numpy as np
import math  as m
import sys
import matplotlib.pyplot as plt

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

    desloc = sys.argv[2]

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

    return a, v


if __name__ == '__main__':
    n = 10
    mass = 2.0
    X  = np.zeros(n)
    K  = np.zeros(n + 1)
    A  = np.zeros((n, n))
    Y  = np.zeros((n, 1))
    Yo = np.zeros((n, 1))
    Xo = np.zeros((n, 1))
    desloc = np.zeros((n, 400)) # matriz que guarda os deslocamentos

    t = np.zeros(400)

    eps = 0.000001

    ci = int(sys.argv[1])

    # Montando o vetor das constantes elasticas
    for i in range(n + 1):
        K[i] = 40 + 2 * (-1)**(i + 1)

    # Montando a matriz A
    for i in range(n):
        A[i][i] = (1 / mass) * (K[i] + K[i + 1])

    for i in range(n - 1):
        A[i][i + 1] = -K[i + 1] * (1 / mass)
        A[i + 1][i] = -K[i + 1] * (1 / mass)

    # Calculando as amtrizes lambda e Q do enunciado
    w, Q = QR_algorithm(A, n, eps)

    # Montando o vetor de condicoes iniciais
    if ci == 1:
        Xo[0][0] = -2.0
        Xo[1][0] = -3.0
        Xo[2][0] = -1.0
        Xo[3][0] = -3.0
        Xo[4][0] = -1.0
        Xo[5][0] = -2.0
        Xo[6][0] = -3.0
        Xo[7][0] = -1.0
        Xo[8][0] = -3.0
        Xo[9][0] = -1.0

    if ci == 2:
        Xo[0][0] =  1.0
        Xo[1][0] = 10.0
        Xo[2][0] = -4.0
        Xo[3][0] =  3.0
        Xo[4][0] = -2.0
        Xo[5][0] = 1.0
        Xo[6][0] = 10.0
        Xo[7][0] = -4.0
        Xo[8][0] = 3.0
        Xo[9][0] = -2.0

    if ci == 3:
        l = max(w.min(), w.max(), key=abs)

        for i in range(n):
            if l == w[i][i]:
                ind = i

        for i in range(n):
            Xo[i][0] = Q[i][ind]

    # Matriz Yo de valores iniciais para Y:
    Qt = Q.transpose()
    Yo = np.dot(Qt, Xo)

    # Atualizando a primeira coluna da matriz de deslocamentos com as condicoes iniciais
    t[0] = 0
    for j in range(n):
        desloc[j][0] = Xo[j][0]
    
    # Calculando as solucoes do sistema Y''(t) + lam * Y(t) = 0, para o intervalo definido e
    # calculando as solucoes de X''(t) + A * X(t) = 0
    for i in range(1, 400, 1):
        for j in range(n):
            Y[j][0] = Yo[j][0] * m.cos(m.sqrt(w[j][j]) * (i + 1) * 0.025)

        X = np.dot(Q, Y)

        for j in range(n):
            desloc[j][i] = X[j][0]

        t[i] = 0.025 * i  # vetor de tempos

    for i in range(n):
        x_coord = t
        y_coord = desloc[i, :]
        plt.plot(x_coord, y_coord)
        plt.xlabel('tempo (s)')
        plt.ylabel('deslocamento (m)')
        plt.title('Evolução da solução do corpo {}'.format(i + 1))
        plt.show()