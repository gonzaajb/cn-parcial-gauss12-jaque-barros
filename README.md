# cn-parcial-gauss12-jaque-barros

Resolvemos un sistema lineal denso de 12x12 con eliminación gaussiana hecha por nosotros. Las filas van de 10⁻⁸ a 10⁵ así que se tuvo mucho cuidado con el pivoteo. Contamos cada flop y comparamos tres estrategias: sin pivoteo, con pivoteo parcial y con pivoteo parcial escalado

## Autoría y asignatura

- Estudiantes: Gonzalo Jaque y Benjamín Barros. 
- Carrera: Ingeniería Civil Informática, UCM
- Asignatura: Computación Numérica INF-321. Prueba Parcial "El pivote"
- Docentes: Juan Manuel Cabezas Contardo y Pedro Pesante
- Talca, noviembre de 2026

## Requisitos

Necesitas Python 3.10 o superior. Lo probamos con Python 3.12.10 en Windows 11 y funcionó sin problemas.

Las dependencias están en `requirements.txt`. En realidad es solo `numpy` y lo usamos nomás como control, la resolución del sistema es código propio. Para la precisión extendida usamos `fractions` que ya viene con Python así que no hay que instalar nada extra para eso|

```bash
python -m venv .venv
.venv\Scripts\activate        # en Linux/Mac se debe usar source .venv/bin/activate
pip install -r requirements.txt
```

## Ejecución

Desde la raíz del repo corre:

```bash
python src/main.py
```

Se demora unos segundos nada más y genera todo lo de `resultados/`

| Archivo | Qué tiene |
|---|---|
| `solucion.csv` | las 12 incógnitas con las cifras que afirmamos y las que verificamos contra la solución exacta |
| `traza_eliminacion.txt` | la matriz aumentada después de cada uno de los 11 pasos, con pivoteo escalado y sin pivoteo |
| `conteo_flops.md` | la tabla de operaciones paso a paso para las tres estrategias, el pivoteo se cuenta aparte |
| `residuo.txt` | el residuo r = b − Ax en aritmética exacta, la norma infinito, el refinamiento, kappa(A) y el control con numpy |

## Estructura de carpetas

```
datos/        A y b con todas sus cifras tal como vienen en el enunciado
src/          gauss.py hace la eliminación, conteo.py cuenta los flops, verificacion.py calcula el
              residuo exacto y el refinamiento, main.py ejecuta todo
resultados/   acá cae todo lo que genera main.py
```

En la rama `pruebas` también está `pruebas/probar_gauss.py`. Ese script revisa el contador y la LU con matrices aleatorias de tamaño 2 a 15.

## Resultado resumido

Esta es la solución que reportamos usando pivoteo parcial escalado más refinamiento iterativo, redondeada a 10 cifras:

| i | xᵢ | i | xᵢ |
|---|---:|---|---:|
| 1 | 0.7384915627 | 7 | 1.473820592 |
| 2 | −1.205937417 | 8 | −0.03849271640 |
| 3 | 0.004387195200 | 9 | 0.9152648037 |
| 4 | −0.09172645380 | 10 | −2.607139485 |
| 5 | 2.846109372 | 11 | 0.05187362940 |
| 6 | −0.006281540900 | 12 | −1.084937162 |

- Flops totales contando eliminación, actualización de b y sustitución regresiva: 1354. La cota de referencia es como 1486 así que quedamos un 8.9% por debajo. Aparte el pivoteo escalado suma 198 comparaciones, 77 divisiones de escala y 8 intercambios.
- La norma infinito del residuo calculado exacto da 2.603e-10 para la solución refinada. Sin refinar queda en 2.291e-10.

Las 16 cifras completas de cada incógnita están en `resultados/solucion.csv`.

## Informe

Todo el análisis de costo, mantisa, condicionamiento y verificación está en [informe.md](informe.md).

## Licencia / uso

Este es un trabajo académico que hicimos en pareja para INF-321 en la Universidad Católica del Maule. Puedes leerlo y correrlo como referencia pero no lo entregues como si fuera tuyo
