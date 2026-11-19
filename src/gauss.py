from conteo import Contador


class Factorizacion: 
    # Guarda todo lo que produce la eliminacion, para usarlo despues
    def __init__(self, n):
        self.n= n
        self.U=None
        self.c=None
        self.L=[[0.0] * n for _ in range(n)]
        self.perm = list(range(n))
        self.max_mult = 0.0
        self.max_mult_pos = None
        self.max_elem_inicial = 0.0
        self.max_elem_durante = 0.0
        self.pivotes = []

    @property
    def crecimiento(self):
        return self.max_elem_durante/ self.max_elem_inicial

def _max_abs(M,desde=0):
    # Función que determina el mayor valor absoluto
    return max(abs(v) for fila in M[desde:] for v in fila[desde:])


def eliminar(A, b, estrategia="escalado", contador=None, traza=None):
    #Deja en cero todo lo que está bajo la diagonal, eligiendo pivotes y contando cada operación

    if estrategia not in ("ninguno", "parcial", "escalado"):
        raise ValueError("estrategia desconocida: " + estrategia)
    n= len(A)
    cont= contador if contador is not None else Contador()
    a=[list(map(float, fila)) for fila in A]
    c=list(map(float, b))
    F=Factorizacion(n)
    F.max_elem_inicial = _max_abs(a)
    F.max_elem_durante = F.max_elem_inicial

    s = None
    if estrategia=="escalado":
        cont.iniciar("escalas")
        s= []
        for i in range(n):
            m=abs(a[i][0])
            for j in range(1, n):
                cont.comp()
                if abs(a[i][j]) > m:
                    m = abs(a[i][j])
            s.append(m)

    for k in range(n-1):
        cont.iniciar(f"paso {k+ 1}")

        p = k
        if estrategia== "parcial":
            mejor= abs(a[k][k])
            for i in range(k+1, n):
                cont.comp()
                if abs(a[i][k])> mejor:
                    mejor,p =abs(a[i][k]), i
        elif estrategia == "escalado":
            cont.div_esc()
            mejor = abs(a[k][k])/s[k]
            for i in range(k +1, n):
                cont.div_esc()
                r = abs(a[i][k])/s[i]
                cont.comp()
                if r > mejor:
                    mejor, p = r, i
        if p != k:
            cont.inter()
            a[k],a[p] = a[p], a[k]
            c[k], c[p] =c[p], c[k]
            F.perm[k], F.perm[p] = F.perm[p],F.perm[k]
            F.L[k], F.L[p] = F.L[p], F.L[k]
            if s is not None:
                s[k], s[p] = s[p], s[k]

        piv = a[k][k]
        if piv==0.0:
            raise ZeroDivisionError(f"pivote nulo en el paso {k+1}")
        F.pivotes.append((k + 1, F.perm[k] +1, piv))
        fila_piv = a[k]

        for i in range(k + 1, n):
            cont.div()
            m = a[i][k] / piv
            F.L[i][k] = m
            if abs(m) > F.max_mult:
                F.max_mult = abs(m)
                F.max_mult_pos = (k + 1, F.perm[i] + 1)
            fila= a[i]
            fila[k] = 0.0
            for j in range(k+ 1, n):
                cont.mul()
                t= m*fila_piv[j]
                cont.resta()
                fila[j] = fila[j] - t
            cont.mul()
            t = m * c[k]
            cont.resta()
            c[i] = c[i]-t

        F.max_elem_durante = max(F.max_elem_durante, _max_abs(a,k +1))
        if traza is not None:
            traza.append({
                "paso": k + 1,
                "pivote_fila_original": F.perm[k] + 1,
                "pivote": piv,
                "perm": [q +1 for q in F.perm],
                "aumentada": [a[i][:] +[c[i]] for i in range(n)],
            })

    F.U= a
    F.c= c
    return F

def sustitucion_regresiva(U, c, contador=None):
    #Despeja las incognitas del sistema triangular, desde la última hasta la primera
    n= len(U)
    cont = contador if contador is not None else Contador()
    cont.iniciar("sustitucion regresiva")
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        acc = c[i]
        for j in range(i + 1, n):
            cont.mul()
            t = U[i][j] * x[j]
            cont.resta()
            acc = acc - t
        cont.div()
        x[i] = acc / U[i][i]
    return x

def resolver_con_factorizacion(F, rhs, contador=None): 
    #Resuelve el sistema con otro lado derecho usando la L y la U ya calculadas
    n = F.n
    cont= contador if contador is not None else Contador()
    cont.iniciar("sustitucion progresiva (L)")
    y= [float(rhs[F.perm[i]]) for i in range(n)]
    for k in range(n - 1):
        for i in range(k + 1, n):
            cont.mul()
            t = F.L[i][k] * y[k]
            cont.resta()
            y[i] = y[i] - t
    return sustitucion_regresiva(F.U, y, cont)


def resolver_transpuesta(F, rhs):
    # Resuelve el sistema con la matriz transpuesta usando la misma L y U, sin calcular la inversa
    n =F.n
    U,L =F.U,F.L
    w= [0.0] * n
    for i in range(n):
        acc= float(rhs[i])
        for j in range(i):
            acc -= U[j][i] * w[j]
        w[i] = acc/U[i][i]
    v= w[:]
    for i in range(n-1,-1,-1):
        acc=v[i]
        for j in range(i+1, n):
            acc-= L[j][i] * v[j]
        v[i]=acc
    z= [0.0] * n
    for i in range(n):
        z[F.perm[i]] = v[i]
    return z


def resolver(A,b, estrategia="escalado", contador=None, traza=None):
    # Resuelve el sistema completo
    cont = contador if contador is not None else Contador()
    F = eliminar(A, b, estrategia, cont, traza)
    x = sustitucion_regresiva(F.U, F.c, cont)
    return x, F
