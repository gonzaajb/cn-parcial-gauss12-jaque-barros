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
