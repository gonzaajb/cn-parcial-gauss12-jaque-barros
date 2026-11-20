# Conteo de operaciones (generado por src/main.py)

Convención: **1 flop = 1 operación aritmética elemental** (mult., div., suma o resta) en doble precisión. Asignar el 0.0 de la subdiagonal, tomar valor absoluto, copiar o intercambiar filas NO son flops. Cada número sale del contador instrumentado en `src/gauss.py` (clase `Contador` de `src/conteo.py`).

Cada fila *paso k* incluye: n-k divisiones (un multiplicador por fila), la actualización de la submatriz activa (columnas k+1..n) y la actualización de b. El sobrecosto del pivoteo va en las tres columnas de la derecha y **no** está sumado en *Flops*.

### Sin pivoteo

| Fase | Mult. | Div. | Mult.+Div. | Sumas | Restas | Sumas+Restas | Flops | Comparaciones | Div. escala | Intercambios |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| paso 1 | 132 | 11 | 143 | 0 | 132 | 132 | 275 | 0 | 0 | 0 |
| paso 2 | 110 | 10 | 120 | 0 | 110 | 110 | 230 | 0 | 0 | 0 |
| paso 3 | 90 | 9 | 99 | 0 | 90 | 90 | 189 | 0 | 0 | 0 |
| paso 4 | 72 | 8 | 80 | 0 | 72 | 72 | 152 | 0 | 0 | 0 |
| paso 5 | 56 | 7 | 63 | 0 | 56 | 56 | 119 | 0 | 0 | 0 |
| paso 6 | 42 | 6 | 48 | 0 | 42 | 42 | 90 | 0 | 0 | 0 |
| paso 7 | 30 | 5 | 35 | 0 | 30 | 30 | 65 | 0 | 0 | 0 |
| paso 8 | 20 | 4 | 24 | 0 | 20 | 20 | 44 | 0 | 0 | 0 |
| paso 9 | 12 | 3 | 15 | 0 | 12 | 12 | 27 | 0 | 0 | 0 |
| paso 10 | 6 | 2 | 8 | 0 | 6 | 6 | 14 | 0 | 0 | 0 |
| paso 11 | 2 | 1 | 3 | 0 | 2 | 2 | 5 | 0 | 0 | 0 |
| sustitucion regresiva | 66 | 12 | 78 | 0 | 66 | 66 | 144 | 0 | 0 | 0 |
| **Total** | **638** | **78** | **716** | **0** | **638** | **638** | **1354** | **0** | **0** | **0** |

### Pivoteo parcial (max |a(i,k)|)

| Fase | Mult. | Div. | Mult.+Div. | Sumas | Restas | Sumas+Restas | Flops | Comparaciones | Div. escala | Intercambios |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| paso 1 | 132 | 11 | 143 | 0 | 132 | 132 | 275 | 11 | 0 | 1 |
| paso 2 | 110 | 10 | 120 | 0 | 110 | 110 | 230 | 10 | 0 | 1 |
| paso 3 | 90 | 9 | 99 | 0 | 90 | 90 | 189 | 9 | 0 | 1 |
| paso 4 | 72 | 8 | 80 | 0 | 72 | 72 | 152 | 8 | 0 | 1 |
| paso 5 | 56 | 7 | 63 | 0 | 56 | 56 | 119 | 7 | 0 | 1 |
| paso 6 | 42 | 6 | 48 | 0 | 42 | 42 | 90 | 6 | 0 | 1 |
| paso 7 | 30 | 5 | 35 | 0 | 30 | 30 | 65 | 5 | 0 | 0 |
| paso 8 | 20 | 4 | 24 | 0 | 20 | 20 | 44 | 4 | 0 | 1 |
| paso 9 | 12 | 3 | 15 | 0 | 12 | 12 | 27 | 3 | 0 | 1 |
| paso 10 | 6 | 2 | 8 | 0 | 6 | 6 | 14 | 2 | 0 | 0 |
| paso 11 | 2 | 1 | 3 | 0 | 2 | 2 | 5 | 1 | 0 | 1 |
| sustitucion regresiva | 66 | 12 | 78 | 0 | 66 | 66 | 144 | 0 | 0 | 0 |
| **Total** | **638** | **78** | **716** | **0** | **638** | **638** | **1354** | **66** | **0** | **9** |

### Pivoteo parcial escalado (max |a(i,k)|/s(i))

| Fase | Mult. | Div. | Mult.+Div. | Sumas | Restas | Sumas+Restas | Flops | Comparaciones | Div. escala | Intercambios |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| escalas | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 132 | 0 | 0 |
| paso 1 | 132 | 11 | 143 | 0 | 132 | 132 | 275 | 11 | 12 | 1 |
| paso 2 | 110 | 10 | 120 | 0 | 110 | 110 | 230 | 10 | 11 | 0 |
| paso 3 | 90 | 9 | 99 | 0 | 90 | 90 | 189 | 9 | 10 | 1 |
| paso 4 | 72 | 8 | 80 | 0 | 72 | 72 | 152 | 8 | 9 | 1 |
| paso 5 | 56 | 7 | 63 | 0 | 56 | 56 | 119 | 7 | 8 | 1 |
| paso 6 | 42 | 6 | 48 | 0 | 42 | 42 | 90 | 6 | 7 | 1 |
| paso 7 | 30 | 5 | 35 | 0 | 30 | 30 | 65 | 5 | 6 | 1 |
| paso 8 | 20 | 4 | 24 | 0 | 20 | 20 | 44 | 4 | 5 | 1 |
| paso 9 | 12 | 3 | 15 | 0 | 12 | 12 | 27 | 3 | 4 | 1 |
| paso 10 | 6 | 2 | 8 | 0 | 6 | 6 | 14 | 2 | 3 | 0 |
| paso 11 | 2 | 1 | 3 | 0 | 2 | 2 | 5 | 1 | 2 | 0 |
| sustitucion regresiva | 66 | 12 | 78 | 0 | 66 | 66 | 144 | 0 | 0 | 0 |
| **Total** | **638** | **78** | **716** | **0** | **638** | **638** | **1354** | **198** | **77** | **8** |

## Comparación con las cotas

| Referencia | Flops | Diferencia vs. este código |
|---|---:|---:|
| Este código (cualquier estrategia) | 1354 | — |
| Mínimo de referencia del enunciado | 1486 | -8.9 % |
| Fórmula cerrada de este algoritmo: 2n³/3 + 3n²/2 − 7n/6 | 1354 | coincide con el contador |
| Término dominante 2n³/3 | 1152 | +17.5 % |

La diferencia de 132 flops con la cota de referencia es exactamente Σ_(p=1)^(11) 2p = 132: la cota recalcula a(i,k) − m·a(k,k) en la columna del pivote, y este código asigna 0.0 directamente porque ese resultado se conoce de antemano.

## Sobrecosto de estabilidad (aparte)

| Estrategia | Comparaciones | Divisiones de escala | Intercambios | Total sobrecosto |
|---|---:|---:|---:|---:|
| Sin pivoteo | 0 | 0 | 0 | 0 |
| Pivoteo parcial (max \|a(i,k)\|) | 66 | 0 | 9 | 75 |
| Pivoteo parcial escalado (max \|a(i,k)\|/s(i)) | 198 | 77 | 8 | 283 |

## Costo del refinamiento iterativo (aparte)

Una pasada reutiliza L y U: 276 flops en doble (132 con L + 144 con U), más el residuo en aritmética racional exacta (144 productos y 144 restas en `Fraction`, que no se cuentan como flops de doble). No se refactoriza nunca.
