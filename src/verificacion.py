# Verificación independiente en precisión extendida. 

import math 
from fractions import Fraction 

from gauss import resolver_con_factorizacion, resolver_transpuesta 

def leer_texto(ruta_a, ruta_b): # leemos el archivo con la matriz A y el vector b
    with open(ruta_a, encoding="utf-8") as f:
        A = [linea.strip().split(",") for linea in f if linea.strip()]
    with open(ruta_b, encoding= "utf-8") as f:
        b = [linea.strip() for linea in f if linea.strip()]
    return A, b

def a_fraction(A_txt, b_txt): #convierte los valores a objetos para trabajar con aritmética
    return [[Fraction(v) for v in fila] for fila in A_txt], [Fraction(v) for v in b_txt]
def a_float(A_txt, b_txt): 
    return [[float(v) for v in fila] for fila in A_txt], [float(v) for v in b_txt]

def residuo_exacto(Aq, bq, x): #Calcula el residuo exacto r = b - A x 
    xq = [Fraction(v) for v in x]
    return [bq[i] - sum(Aq[i][j]* xq[j] for j in range(len(xq))) for i in range(len(bq))]


def residuo_doble(A, b, x): 
    return [b[i] - sum(A[i][j] * x[j] for j in range(len(x))) for i in range(len(b))]


def norma_inf(v): #define valor absoluto máximo de sus componentes
    return max(abs(t) for t in v)


def residuo_relativo(Aq, bq, x, r): # buscamos comparar con el residuo relativo en doble precisión
    nA = max(sum(abs(v) for v in fila) for fila in Aq)
    return float(norma_inf(r)) / (float(nA) * norma_inf(x) + float(norma_inf(bq)))

def refinar(F, Aq, bq, x0, pasadas=3): # El historial guarda en cada pasada ||r||_inf ,antes de corregirlo, enconteces usamos ||d||_inf/||x||_inf y la corrección componente a componente.
    x= list(x0)
    historial = []
    for k in range(1, pasadas + 1):
        r = residuo_exacto(Aq, bq, x)
        rf = [float(t) for t in r]               
        d = resolver_con_factorizacion(F, rf)
        historial.append({
            "pasada": k,
            "norma_r_antes": float(norma_inf(r)),
            "correccion_rel": norma_inf(d) / norma_inf(x),
            "d": d,
            "x_antes": list(x),
        })
        x = [x[i] + d[i] for i in range(len(x))]
    r = residuo_exacto(Aq, bq, x)
    historial.append({"pasada": "final", "norma_r_antes": float(norma_inf(r)),
                      "correccion_rel": None, "d": None, "x_antes": list(x)})
    return x, historial

def solucion_exacta(Aq, bq):  # devuelve la solución exacta x como una lista de objetos Fraction.
    n =len(Aq)
    a= [fila[:] + [bq[i]] for i, fila in enumerate(Aq)]
    for k in range(n):
        p = next(i for i in range(k, n) if a[i][k] != 0)
        a[k], a[p] = a[p], a[k]
        for i in range(k + 1, n):
            m = a[i][k] / a[k][k]
            if m:
                for j in range(k, n + 1):
                    a[i][j] -= m * a[k][j]
    x = [Fraction(0)] * n
    for i in range(n - 1, -1, -1):
        x[i] = (a[i][n] - sum(a[i][j] * x[j] for j in range(i + 1, n))) / a[i][i]
    return x


def determinante_exacto(Aq): # el determinante exacto de A en aritmética racional, usando eliminación gaussiana sin pivoteo
    n = len(Aq)
    a = [fila[:] for fila in Aq]
    det = Fraction(1)
    for k in range(n):
        p = next(i for i in range(k, n) if a[i][k] != 0)
        if p != k:
            a[k], a[p] = a[p], a[k]
            det = -det
        det *= a[k][k]
        for i in range(k + 1, n):
            m = a[i][k] / a[k][k]
            for j in range(k, n):
                a[i][j] -= m * a[k][j]
    return det


def digitos_correctos(x, xq): # Devuelve una lista con los dígitos correctos para cada componente.
    out = []
    for xi, ri in zip(x, xq):
        err = abs(Fraction(xi) - ri) / abs(ri)
        out.append(17.0 if err == 0 else -math.log10(float(err)))
    return out

def norma_1(A):
    n = len(A)
    return max(sum(abs(A[i][j]) for i in range(n)) for j in range(n)) #regresa la norma 1 de la matriz A que es el máximo de las sumas de los valores absolutos de cada columna


def estimar_inv_norma_1(F): #Devuelve una estimación
    n = F.n
    x = [1.0 / n] * n
    est = 0.0
    j_ant = -1
    for _ in range(5):
        y = resolver_con_factorizacion(F, x)
        est = sum(abs(v) for v in y)
        xi = [1.0 if v >= 0 else -1.0 for v in y]
        z = resolver_transpuesta(F, xi)
        j = max(range(n), key=lambda t: abs(z[t]))
        if abs(z[j]) <= sum(z[t] * x[t] for t in range(n)) or j == j_ant:
            break
        x = [0.0] * n
        x[j] = 1.0
        j_ant = j
    alt = [(-1) ** i * (1 + i / (n - 1)) for i in range(n)] # Agregamos pasos para protegernos de la dependencia de la estimación 
    y = resolver_con_factorizacion(F, alt)
    est_alt = 2 * sum(abs(v) for v in y) / (3 * n)
    return max(est, est_alt)

def casi_dependencia(A, i, j, k): # buscamos una combinación lineal de las filas j y k que se aproxime a la fila i
    fi=[Fraction(v) for v in A[i]]
    fj=[Fraction(v) for v in A[j]]
    fk=[Fraction(v) for v in A[k]]

    def dot(u, v):
        return sum(p * q for p, q in zip(u, v))

    # Calculamos los productos punto necesarios para resolver el sistema de ecuaciones normales que nos dará los coeficientes alfa y beta
    g11, g12, g22 = dot(fj, fj), dot(fj, fk), dot(fk, fk)
    h1, h2 = dot(fj, fi), dot(fk, fi)
    det = g11 * g22 - g12 * g12
    alfa = (h1 * g22 - h2 * g12) / det
    beta = (g11 * h2 - g12 * h1) / det
    resto = [fi[t] - alfa * fj[t] - beta * fk[t] for t in range(len(fi))]
    rel = math.sqrt(float(dot(resto, resto))) / math.sqrt(float(dot(fi, fi)))
    return float(alfa), float(beta), rel # Retornamos los valores de alfa, beta y la relación entre la norma del residuo y la norma de la fila i
