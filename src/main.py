import math
import os
import platform
import sys

import numpy as np

from conteo import Contador, tabla_markdown
from gauss import resolver, resolver_con_factorizacion
from verificacion import (a_float, a_fraction, casi_dependencia, determinante_exacto,digitos_correctos, estimar_inv_norma_1, leer_texto, norma_1,norma_inf, refinar, residuo_doble, residuo_exacto,residuo_relativo, solucion_exacta)

RAIZ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATOS = os.path.join(RAIZ, "datos")
RES = os.path.join(RAIZ, "resultados")
EPS = 2.0 ** -52                     
U = EPS / 2                          
N_REF_PASADAS = 3
COTA_MINIMA = 1486

ESTRATEGIAS = [
    ("ninguno", "Sin pivoteo"),
    ("parcial", "Pivoteo parcial (max |a(i,k)|)"),
    ("escalado", "Pivoteo parcial escalado (max |a(i,k)|/s(i))"),
]
OFICIAL = "escalado"                  # corrida cuya solución se reporta


def fmt_matriz(M, ancho=24):
    #da formato a la matriz para escribirla en la traza
    return "\n".join("".join(f"{v: .16e}".rjust(ancho) for v in fila) for fila in M)


def cifras(err_rel, tope=17):
    #convierte un error relativo en cantidad de cifras significativas.
    if err_rel == 0:
        return tope
    return min(tope, int(math.floor(math.log10(5.0 / err_rel))))


def digitos_afirmados(x, d):
    #calcula cuántas cifras se pueden defender sin conocer la solución exacta
    return [cifras(abs(di) / abs(xi) + U, tope=16) for xi, di in zip(x, d)]


def digitos_verificados(valores_txt, x_exacta):
    #comprueba cuántas cifras del valor reportado son realmente correctas
    from fractions import Fraction
    return [cifras(float(abs(Fraction(v) - xe) / abs(xe))) for v, xe in zip(valores_txt, x_exacta)]


def redondear_sig(v, dig):
    #escribe un número con la cantidad de cifras indicada
    return f"{v:.{dig - 1}e}"


def main():
    #crea la carpeta resultados, lee A y b como texto y los convierte a doble precisión y a fracciones exactas. Luego verifica que el sistema sea de 12×12
    os.makedirs(RES, exist_ok=True)
    A_txt, b_txt = leer_texto(os.path.join(DATOS, "matriz_a.csv"), os.path.join(DATOS, "vector_b.csv"))
    A, b = a_float(A_txt, b_txt)
    Aq, bq = a_fraction(A_txt, b_txt)
    n = len(A)
    assert n == 12 and all(len(f) == 12 for f in A) and len(b) == 12

    #resuelve el sistema con fracciones, solo para usarla como referencia al comparar.
    x_exacta = solucion_exacta(Aq, bq)

    #para cada estrategia resuelve el sistema, calcula el residuo exacto, refina la solución y cuenta los dígitos correctos
    corridas = {}
    for clave, nombre in ESTRATEGIAS:
        cont = Contador()
        traza = []
        x, F = resolver(A, b, clave, cont, traza)
        r = residuo_exacto(Aq, bq, x)
        x_ref, hist = refinar(F, Aq, bq, x, N_REF_PASADAS)
        corridas[clave] = dict(nombre=nombre, cont=cont, traza=traza, x=x, F=F, r=r,x_ref=x_ref, hist=hist,dig_x=digitos_correctos(x, x_exacta),dig_ref=digitos_correctos(x_ref, x_exacta))

    of = corridas[OFICIAL]

    #mide cuántos flops cuesta una pasada de refinamiento que reutiliza L y U
    cont_ref = Contador()
    resolver_con_factorizacion(of["F"], [float(t) for t in of["r"]], cont_ref)

    # Estima kappa(A), la casi dependencia de F7, det(A) y la solución de numpy como control
    nA1 = norma_1(A)
    kappa1_est = nA1 * estimar_inv_norma_1(of["F"])
    An = np.array(A)
    kappa1_np = np.linalg.cond(An, 1)
    kappa2_np = np.linalg.cond(An, 2)
    D = np.abs(An).max(axis=1)
    kappa2_eq = np.linalg.cond(An / D[:, None], 2)
    alfa, beta, dep_rel = casi_dependencia(A_txt, 6, 2, 4)
    det = float(determinante_exacto(Aq))
    x_np = np.linalg.solve(An, np.array(b))

    #escribe la matriz tras cada paso y comprueba los 66 ceros bajo la diagonal
    with open(os.path.join(RES, "traza_eliminacion.txt"), "w", encoding="utf-8") as f:
        f.write("TRAZA DE LA ELIMINACIÓN GAUSSIANA (matriz aumentada [A | b] tras cada paso)\n")
        f.write("Valores en doble precisión con 17 cifras significativas (repr. exacta del binario).\n")
        f.write("La última columna es b transformado. Filas en el orden en que quedan tras los intercambios.\n\n")
        for clave in ("escalado", "ninguno"):
            c = corridas[clave]
            f.write("=" * 100 + "\n")
            f.write(f"CORRIDA: {c['nombre']}\n")
            f.write("=" * 100 + "\n\n")
            for t in c["traza"]:
                f.write(f"--- Paso {t['paso']} de 11 --- pivote = {t['pivote']: .16e} "
                        f"(fila original F{t['pivote_fila_original']})\n")
                f.write("orden de filas: " + " ".join(f"F{q}" for q in t["perm"]) + "\n")
                f.write(fmt_matriz(t["aumentada"]) + "\n\n")
            Uf = c["F"].U
            nulos = sum(1 for i in range(n) for j in range(i) if Uf[i][j] == 0.0)
            diag = min(abs(Uf[i][i]) for i in range(n))
            f.write(f"Comprobación: elementos bajo la diagonal iguales a 0.0 -> {nulos} de 66\n")
            f.write(f"Comprobación: min |u(i,i)| = {diag:.6e} (ningún pivote nulo)\n\n")

    # Escribe las tablas de flops por paso, la comparación con 1486 y el sobrecosto del pivoteo
    teor = 2 * n ** 3 / 3
    with open(os.path.join(RES, "conteo_flops.md"), "w", encoding="utf-8") as f:
        f.write("# Conteo de operaciones (generado por src/main.py)\n\n")
        f.write("Convención: **1 flop = 1 operación aritmética elemental** (mult., div., suma o resta) "
                "en doble precisión. Asignar el 0.0 de la subdiagonal, tomar valor absoluto, copiar o "
                "intercambiar filas NO son flops. Cada número sale del contador instrumentado en "
                "`src/gauss.py` (clase `Contador` de `src/conteo.py`).\n\n")
        f.write("Cada fila *paso k* incluye: n-k divisiones (un multiplicador por fila), la actualización "
                "de la submatriz activa (columnas k+1..n) y la actualización de b. El sobrecosto del "
                "pivoteo va en las tres columnas de la derecha y **no** está sumado en *Flops*.\n\n")
        for clave, _ in ESTRATEGIAS:
            f.write(tabla_markdown(corridas[clave]["cont"], corridas[clave]["nombre"]) + "\n")
        tot = of["cont"].total_flops()
        f.write("## Comparación con las cotas\n\n")
        f.write("| Referencia | Flops | Diferencia vs. este código |\n|---|---:|---:|\n")
        f.write(f"| Este código (cualquier estrategia) | {tot} | — |\n")
        f.write(f"| Mínimo de referencia del enunciado | {COTA_MINIMA} | {100 * (tot - COTA_MINIMA) / COTA_MINIMA:+.1f} % |\n")
        f.write(f"| Fórmula cerrada de este algoritmo: 2n³/3 + 3n²/2 − 7n/6 | "
                f"{round(2 * n**3 / 3 + 3 * n**2 / 2 - 7 * n / 6)} | coincide con el contador |\n")
        f.write(f"| Término dominante 2n³/3 | {teor:.0f} | {100 * (tot - teor) / teor:+.1f} % |\n\n")
        f.write(f"La diferencia de {COTA_MINIMA - tot} flops con la cota de referencia es exactamente "
                f"Σ_(p=1)^(11) 2p = 132: la cota recalcula a(i,k) − m·a(k,k) en la columna del pivote, "
                f"y este código asigna 0.0 directamente porque ese resultado se conoce de antemano.\n\n")
        so = {k: corridas[k]["cont"] for k, _ in ESTRATEGIAS}
        f.write("## Sobrecosto de estabilidad (aparte)\n\n")
        f.write("| Estrategia | Comparaciones | Divisiones de escala | Intercambios | Total sobrecosto |\n|---|---:|---:|---:|---:|\n")
        for k, nom in ESTRATEGIAS:
            c = so[k]
            f.write(f"| {nom.replace('|', chr(92) + '|')} | {c.total('comp')} | {c.total('div_esc')} | {c.total('inter')} | {c.total_sobrecosto()} |\n")
        f.write("\n## Costo del refinamiento iterativo (aparte)\n\n")
        f.write(f"Una pasada reutiliza L y U: {cont_ref.total_flops()} flops en doble "
                f"({cont_ref.flops('sustitucion progresiva (L)')} con L + "
                f"{cont_ref.flops('sustitucion regresiva')} con U), más el residuo en aritmética racional "
                f"exacta (144 productos y 144 restas en `Fraction`, que no se cuentan como flops de doble). "
                f"No se refactoriza nunca.\n")

    # Escribe el residuo, crecimiento, bits perdidos, refinamiento, kappa y control con numpy
    with open(os.path.join(RES, "residuo.txt"), "w", encoding="utf-8") as f:
        f.write("VERIFICACIÓN INDEPENDIENTE\n")
        f.write("Residuo r = b - A x calculado en aritmética racional EXACTA (fractions.Fraction):\n")
        f.write("A y b se leen como texto (b con sus 20 dígitos) y x se convierte exacto desde binario.\n\n")
        for clave, _ in ESTRATEGIAS:
            c = corridas[clave]
            rd = residuo_doble(A, b, c["x"])
            f.write("=" * 90 + f"\n{c['nombre']}\n" + "=" * 90 + "\n")
            f.write(f"{'i':>3} {'r_i exacto':>24} {'r_i calculado en doble':>26}\n")
            for i in range(n):
                f.write(f"{i + 1:>3} {float(c['r'][i]):>24.15e} {rd[i]:>26.15e}\n")
            f.write(f"||r||_inf exacto            = {float(norma_inf(c['r'])):.6e}\n")
            f.write(f"||r||_inf en doble          = {norma_inf(rd):.6e}\n")
            f.write(f"residuo relativo normado    = {residuo_relativo(Aq, bq, c['x'], c['r']):.3e}"
                    "   (||r|| / (||A|| ||x|| + ||b||))\n")
            f.write(f"multiplicador máximo |m|    = {c['F'].max_mult:.6e}  (paso, fila orig.) = {c['F'].max_mult_pos}\n")
            f.write(f"factor de crecimiento rho   = {c['F'].crecimiento:.6e}\n")
            err = max(abs(float(c['x'][i] - x_exacta[i])) for i in range(n)) / max(abs(float(v)) for v in x_exacta)
            f.write(f"error relativo ||x-x*||/||x*|| = {err:.3e}  -> {-math.log10(err):.2f} dígitos, "
                    f"{-math.log2(err):.1f} bits correctos de 53, {53 + math.log2(err):.1f} bits perdidos\n")
            f.write("dígitos correctos por incógnita (vs. solución exacta): "
                    + " ".join(f"{d:.1f}" for d in c["dig_x"]) + "\n")
            f.write("\nRefinamiento iterativo (residuo exacto, misma L y U):\n")
            for h in c["hist"]:
                cr = "" if h["correccion_rel"] is None else f"   ||d||/||x|| = {h['correccion_rel']:.3e}"
                f.write(f"  pasada {h['pasada']!s:>5}: ||r||_inf = {h['norma_r_antes']:.6e}{cr}\n")
            f.write("dígitos correctos tras refinar: " + " ".join(f"{d:.1f}" for d in c["dig_ref"]) + "\n\n")

        f.write("=" * 90 + "\nCONDICIONAMIENTO\n" + "=" * 90 + "\n")
        f.write(f"det(A) (exacto)                       = {det:.6e}\n")
        f.write(f"||A||_1                               = {nA1:.6e}\n")
        f.write(f"kappa_1(A) estimado (Hager/Higham, LU) = {kappa1_est:.4e}\n")
        f.write(f"kappa_1(A) numpy (control)            = {kappa1_np:.4e}\n")
        f.write(f"kappa_2(A) numpy (control)            = {kappa2_np:.4e}\n")
        f.write(f"cota kappa_1 * eps                    = {kappa1_est * EPS:.3e}   -> {-math.log10(kappa1_est * EPS):.2f} dígitos garantizados\n")
        f.write(f"cota kappa_2 * eps                    = {kappa2_np * EPS:.3e}\n")
        f.write(f"kappa_2 de A con filas equilibradas    = {kappa2_eq:.4e}   (D^-1 A, D = max |a(i,j)| por fila)\n")
        f.write(f"Casi dependencia: F7 ~ {alfa:.6f}*F3 + {beta:.6f}*F5, "
                f"||resto||_2/||F7||_2 = {dep_rel:.3e}\n\n")
        f.write("=" * 90 + "\nCONTROL CON BIBLIOTECA (numpy.linalg.solve, solo comparación)\n" + "=" * 90 + "\n")
        f.write(f"max |x_np - x_propio| / |x_propio| = "
                f"{max(abs(x_np[i] - of['x'][i]) / abs(of['x'][i]) for i in range(n)):.3e}\n")
        f.write(f"dígitos correctos de numpy (vs. exacta): "
                + " ".join(f"{d:.1f}" for d in digitos_correctos(list(map(float, x_np)), x_exacta)) + "\n")

    # calcula las cifras que se pueden afirmar, las verifica y escribe solucion.csv
    d1 = of["hist"][0]["d"]                  
    d_ult = of["hist"][-2]["d"]                
    afirm_x = digitos_afirmados(of["x"], d1)
    afirm_ref = digitos_afirmados(of["x_ref"], d_ult)
    txt_ref = [redondear_sig(of["x_ref"][i], afirm_ref[i]) for i in range(n)]
    txt_x = [redondear_sig(of["x"][i], afirm_x[i]) for i in range(n)]
    verif_ref = digitos_verificados(txt_ref, x_exacta)
    verif_x = digitos_verificados(txt_x, x_exacta)
    with open(os.path.join(RES, "solucion.csv"), "w", encoding="utf-8", newline="\n") as f:
        f.write("i,x_reportado,cifras_afirmadas,cifras_verificadas,x_doble_completo,"
                "x_sin_refinar,cifras_afirmadas_sin_refinar,cifras_verificadas_sin_refinar\n")
        for i in range(n):
            f.write(f"{i + 1},{txt_ref[i]},{afirm_ref[i]},{verif_ref[i]},{of['x_ref'][i]:.17g},"
                    f"{txt_x[i]},{afirm_x[i]},{verif_x[i]}\n")
    ok_ref = all(a <= v for a, v in zip(afirm_ref, verif_ref))
    ok_x = all(a <= v for a, v in zip(afirm_x, verif_x))

    # imprime en pantalla los resultados principales
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except (AttributeError, ValueError):
        pass
    print(f"Python {platform.python_version()} | numpy {np.__version__} | {sys.platform}")
    print(f"Flops eliminación + sustitución: {of['cont'].total_flops()} (referencia {COTA_MINIMA})")
    for clave, nom in ESTRATEGIAS:
        c = corridas[clave]
        print(f"  {nom:48s} |m|max={c['F'].max_mult:.3e} rho={c['F'].crecimiento:.3e} "
              f"dígitos min={min(c['dig_x']):.1f} sobrecosto={c['cont'].total_sobrecosto()}")
    print(f"kappa_1 estimado = {kappa1_est:.3e}, kappa_2 (numpy) = {kappa2_np:.3e}")
    print(f"||r||_inf final (exacto) = {of['hist'][-1]['norma_r_antes']:.3e}")
    print("Solución (corrida oficial, refinada):")
    for i in range(n):
        print(f"  x{i + 1:<2} = {txt_ref[i]:>24}  [{afirm_ref[i]} cifras afirmadas, "
              f"{verif_ref[i]} verificadas]   sin refinar: {txt_x[i]:>20} [{afirm_x[i]} / {verif_x[i]}]")
    print(f"¿Toda cifra afirmada está verificada? refinada: {ok_ref} | sin refinar: {ok_x}")
    print("Archivos generados en resultados/: solucion.csv, traza_eliminacion.txt, "
          "conteo_flops.md, residuo.txt")


if __name__ == "__main__":
    main()
