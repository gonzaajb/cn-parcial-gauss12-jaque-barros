# Instrumentación del costo operacional Convención de conteo: 1 flop = 1 operación aritmética elemental de punto flotante, es decir, una multiplicación, una división, una suma o una resta.

#Lo que no se cuenta como flop:
#    1 - asignar un 0.0 exacto en la subdiagonal
#    2 - el valor absoluto
#    3 - copiar o leer datos. 

#El sobrecosto de estabilidad se cuenta aparte en campos propios, para no mezclarlo con el costo de la eliminación:
#    - comparaciones para buscar el pivote,
#    - divisiones |a(i,k)| / s(i) del pivoteo escalado,
#    - comparaciones para calcular las escalas s(i),
#    - intercambios de filas. 

# Cada operación se registra con una llamada al contador justo donde se hace en gauss.py, ningún total sale de una fórmula escrita a mano.


CAMPOS_FLOP = ("mul", "div", "suma", "resta")
CAMPOS_PIVOTE = ("comp", "div_esc", "inter")


class Contador: #contador de operaciones
    

    def __init__(self): #función de inicialización
        self.fases={}          
        self.orden=[]          
        self.fase =None
        self.iniciar("preparacion")

    def iniciar(self, nombre): #función para iniciar una fase de conteo
        if nombre not in self.fases:
            self.fases[nombre] = {c: 0 for c in CAMPOS_FLOP + CAMPOS_PIVOTE}
            self.orden.append(nombre)
        self.fase = nombre

    # Aritmética de la eliminación
    def mul(self, k=1): 
        self.fases[self.fase]["mul"] += k
    def div(self, k=1):
        self.fases[self.fase]["div"] += k
    def suma(self, k=1):
        self.fases[self.fase]["suma"] += k
    def resta(self, k=1):
        self.fases[self.fase]["resta"] += k
    def comp(self, k=1):
        self.fases[self.fase]["comp"] += k
    def div_esc(self, k=1):
        self.fases[self.fase]["div_esc"] += k
    def inter(self, k=1):
        self.fases[self.fase]["inter"] += k
    def flops(self, fase): 
        f = self.fases[fase]
        return sum(f[c] for c in CAMPOS_FLOP)
    def total(self, campo): #total de un campo 
        return sum(f[campo] for f in self.fases.values())
    def total_flops(self): # total de flops 
        return sum(self.total(c) for c in CAMPOS_FLOP)
    def total_sobrecosto(self): # total de sobrecosto de pivoteo 
        return sum(self.total(c) for c in CAMPOS_PIVOTE)


def tabla_markdown(cont, titulo): #función para generar una tabla en formato markdown con los resultados del conteo de operaciones
    lin = [f"### {titulo}", ""]
    lin.append("| Fase | Mult. | Div. | Mult.+Div. | Sumas | Restas | Sumas+Restas | Flops | Comparaciones | Div. escala | Intercambios |")
    lin.append("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    for nombre in cont.orden:
        f = cont.fases[nombre]
        if sum(f.values()) == 0:
            continue
        md = f["mul"] + f["div"]
        sr = f["suma"] + f["resta"]
        lin.append(
            f"| {nombre} | {f['mul']} | {f['div']} | {md} | {f['suma']} | {f['resta']} | {sr} | "
            f"{md + sr} | {f['comp']} | {f['div_esc']} | {f['inter']} |"
        )
    T = {c: cont.total(c) for c in CAMPOS_FLOP + CAMPOS_PIVOTE}
    lin.append(
        f"| **Total** | **{T['mul']}** | **{T['div']}** | **{T['mul'] + T['div']}** | **{T['suma']}** | "
        f"**{T['resta']}** | **{T['suma'] + T['resta']}** | **{cont.total_flops()}** | **{T['comp']}** | "
        f"**{T['div_esc']}** | **{T['inter']}** |"
    )
    lin.append("")
    return "\n".join(lin)
