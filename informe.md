# Informe técnico — El pivote: eliminación gaussiana de un sistema denso de 12×12

Todos los números de este informe salen de correr `python src/main.py`. Los archivos que genera quedan en la carpeta `resultados/`.

## 1. Identificación

- Estudiantes: Gonzalo Jaque · Benjamín Barros
- Sección: S2
- Asignatura: Computación Numérica (INF-321) — Ingeniería Civil Informática, UCM
- Docentes: Juan Manuel Cabezas Contardo · Pedro Pesante
- Fecha: Talca, noviembre de 2026
- Repositorio: https://github.com/gonzaajb/cn-parcial-gauss12-jaque-barros

---

## 2. Planteamiento

Teníamos que resolver el sistema **Ax = b** con 12 ecuaciones y 12 incógnitas. La matriz A es **densa**, o sea, ninguno de sus 144 números es cero, y el vector b viene con 20 cifras significativas.

Como no hay ceros, no hay atajos: la eliminación gaussiana va a costar sí o sí del orden de 2n³/3 operaciones. Lo que sí podemos hacer es no gastar operaciones de más, es decir, no repetir cálculos.

El verdadero problema del sistema es que sus filas tienen tamaños muy distintos:

| Filas | Tamaño aproximado de sus números |
|---|---|
| F1 | 10⁻⁵ (y a(1,1) = 1.3407×10⁻⁸) |
| F5, F9 | 10⁻⁵ |
| F3, F7, F11 | 10⁻³ – 10⁻² |
| F4, F10 | 10³ |
| F2, F8, F12 | 10⁴ – 10⁵ |
| F6 | 10⁵ – 10⁶ |

Entre el número más chico y el más grande hay unos **14 órdenes de magnitud**. ¿Por qué eso es un problema? El computador guarda cada número con 53 bits de mantisa (más o menos 16 cifras decimales). Si a un número chico le restamos uno gigante, para hacer la resta el computador tiene que "alinear" los dos números, y las últimas cifras del número chico simplemente no caben y se pierden. Esa información no se puede recuperar después.

Además hay dos complicaciones extra:

- El primer pivote, a(1,1), es pequeño.
- La fila F7 es casi una combinación de F3 y F5 (casi, pero no del todo).

Todo esto hace que el número de condición sea enorme: κ₂(A) ≈ 3.08×10¹³.

---

## 3. Estrategia de eliminación

La función `eliminar` (en `src/gauss.py`) hace lo siguiente en cada paso k = 1, …, 11:

1. **Elige el pivote** según la estrategia que se esté probando:
   - sin pivoteo (se usa el que esté en la diagonal);
   - pivoteo parcial (se elige el mayor |a(i,k)| de la columna);
   - pivoteo parcial escalado (se elige el mayor |a(i,k)|/s(i), donde s(i) es el número más grande de la fila i). Las escalas s(i) se calculan **una sola vez** al principio, con la matriz original.
2. Calcula un solo multiplicador por fila, m(i,k) = a(i,k)/a(k,k), y lo usa para toda la fila y también para b. No se vuelve a calcular en cada columna.
3. Solo recorre las columnas que todavía importan (de k+1 a n). Ya sabemos que a(i,k) tiene que quedar en cero, así que le ponemos 0.0 directamente en vez de calcular a(i,k) − m·a(k,k). Con esto nos ahorramos 2(n−k) operaciones por paso y además el cero queda exacto. Si lo calculáramos, por el redondeo podría quedar algo como 10⁻²⁰ en vez de 0.
4. Actualiza b en el mismo ciclo, con el mismo multiplicador.
5. Guarda los multiplicadores en L (junto con los intercambios de filas). Así, si después hay que resolver con otro lado derecho, como en el refinamiento iterativo, usamos `resolver_con_factorizacion` y reutilizamos L y U sin factorizar de nuevo.

Al final se hace la sustitución regresiva con U. **Nunca calculamos la inversa A⁻¹.** Incluso para estimar κ(A) usamos la misma factorización LU, con el estimador de Hager–Higham (`verificacion.estimar_inv_norma_1`).

**¿Por qué decimos que no gastamos operaciones de más?** Porque cada operación que hace el código sirve para una de estas cuatro cosas: calcular un multiplicador que se necesita, actualizar un número de la parte de la matriz que todavía falta, actualizar b, o despejar una incógnita. No operamos sobre los ceros que ya hicimos, no repetimos multiplicadores y no factorizamos dos veces.

---

## 4. Conteo de operaciones

**Qué contamos como "flop":** cualquier multiplicación, división, suma o resta en doble precisión. No contamos poner un 0.0, sacar valor absoluto, copiar o intercambiar filas.

Los números no los sacamos a mano. Los cuenta el objeto `Contador` (`src/conteo.py`), que suma 1 justo en la línea de `gauss.py` donde se hace cada operación. En la rama `pruebas` revisamos con matrices aleatorias de tamaño 2 a 15 que el contador siempre da lo mismo que la fórmula.

Este es el conteo paso a paso. Es igual para las tres estrategias, porque el pivoteo lo contamos aparte:

| Paso | Mult. | Div. | **Mult.+Div.** | Restas | **Sumas+Restas** | **Flops** |
|---|---:|---:|---:|---:|---:|---:|
| 1 | 132 | 11 | 143 | 132 | 132 | 275 |
| 2 | 110 | 10 | 120 | 110 | 110 | 230 |
| 3 | 90 | 9 | 99 | 90 | 90 | 189 |
| 4 | 72 | 8 | 80 | 72 | 72 | 152 |
| 5 | 56 | 7 | 63 | 56 | 56 | 119 |
| 6 | 42 | 6 | 48 | 42 | 42 | 90 |
| 7 | 30 | 5 | 35 | 30 | 30 | 65 |
| 8 | 20 | 4 | 24 | 20 | 20 | 44 |
| 9 | 12 | 3 | 15 | 12 | 12 | 27 |
| 10 | 6 | 2 | 8 | 6 | 6 | 14 |
| 11 | 2 | 1 | 3 | 2 | 2 | 5 |
| Sustitución regresiva | 66 | 12 | 78 | 66 | 66 | 144 |
| **Total** | **638** | **78** | **716** | **638** | **638** | **1354** |

**De dónde salen estos números:** en el paso k quedan p = 12 − k filas bajo el pivote. Para cada una hay que hacer 1 división (el multiplicador), p multiplicaciones en la matriz y 1 en b, y la misma cantidad de restas. En total: p divisiones, p² + p multiplicaciones y p² + p restas. Sumando todos los pasos y la sustitución regresiva queda:

  **2n³/3 + 3n²/2 − 7n/6 = 1152 + 216 − 14 = 1354 flops.**

**Comparación con el enunciado.** El enunciado da como referencia ≈ 1486 flops. Nosotros llegamos a 1354, o sea, un **8.9 % menos**, que está dentro del ±10 % permitido. La diferencia es justo 132 = n(n−1) operaciones, y se explica fácil: la referencia calcula a(i,k) − m·a(k,k) para cada elemento que se hace cero (1 multiplicación + 1 resta), y nosotros ponemos el cero directamente.

**Costo extra del pivoteo** (va aparte, no está sumado en los 1354):

| Estrategia | Comparaciones | Divisiones de escala | Intercambios | Total |
|---|---:|---:|---:|---:|
| Sin pivoteo | 0 | 0 | 0 | 0 |
| Parcial | 66 | 0 | 9 | 75 |
| Parcial escalado | 198 (132 para las escalas + 66 para buscar) | 77 | 8 | 283 |

El parcial coincide con lo que dice el enunciado ("1486 + 66 comparaciones"). El escalado necesita 132 comparaciones más (para encontrar el máximo de cada fila, una sola vez) y 77 divisiones |a(i,k)|/s(i).

**Refinamiento iterativo:** cada vuelta cuesta 276 flops en doble precisión (132 con L y 144 con U), más el cálculo del residuo en aritmética exacta. No se vuelve a factorizar.

---

## 5. Análisis de mantisa y pivoteo

Corrimos el mismo código con las tres estrategias y comparamos cada resultado con la **solución exacta** del sistema. Esa solución exacta la calculamos con fracciones (`fractions.Fraction`) y la usamos solo para revisar, no para resolver.

| Corrida | Multiplicador más grande \|m\| | Crecimiento ρ | Error relativo ‖x−x*‖∞/‖x*‖∞ | Cifras correctas (peor – mejor incógnita) | Bits buenos (de 53) | Bits perdidos |
|---|---:|---:|---:|---|---:|---:|
| Sin pivoteo | 2.75×10¹³ (paso 1, F6) | 2.64×10³ | 2.08×10⁻¹⁰ | 7.1 – 10.2 | 32.2 | **20.8** |
| Parcial | 0.620 | 1.00 | 5.52×10⁻¹³ | 9.7 – 12.8 | 40.7 | 12.3 |
| Parcial escalado | 2.47×10⁹ (paso 3, F6) | 1.78 | 2.36×10⁻¹³ | 10.1 – 13.1 | 41.9 | **11.1** |

- **ρ (factor de crecimiento):** cuánto crecen los números durante la eliminación. Es el número más grande que aparece en el proceso dividido por el más grande de A. Si ρ es grande, es mala señal.
- **Bits perdidos:** 53 + log₂(error relativo). Mientras más bits se pierden, menos cifras buenas quedan.

### 5.1 ¿Por qué no sirve usar a(1,1) como pivote?

En la primera columna está a(1,1) = 1.3407×10⁻⁸ y también a(6,1) = 3.69×10⁵. Si usamos a(1,1) como pivote, el multiplicador para la fila 6 es m(6,1) ≈ 2.75×10¹³, un número gigante.

Entonces, en el paso 1, a la fila 6 se le resta m·(fila 1). Eso es más o menos 2.75×10¹³ · 5×10⁻⁵ ≈ 10⁹, mientras que los números originales de la fila 6 son del orden de 10⁵ a 10⁶. Es como sumarle unos pocos pesos a mil millones: los números originales de la fila 6 quedan entre 10 y 13 bits "por debajo" y se pierden al hacer la resta. En la práctica, la fila 6 se transforma casi en una copia de la fila 1 y pierde su propia información. Lo mismo les pasa a F2, F4, F8, F10 y F12.

Esto se nota en los resultados: la diagonal de U sin pivoteo va desde 1.3×10⁻⁸ hasta 3.7×10⁷, y el crecimiento es ρ = 2.6×10³ ≈ 2¹¹·⁴, o sea, unos 11 bits. Esos 11 bits son casi exactamente la diferencia que medimos entre no pivotear y pivotear escalado (20.8 − 11.1 ≈ 9.7 bits).

El sistema igual "se resuelve" sin pivoteo, porque ningún pivote llega a ser cero, pero se botan unas 3 cifras decimales sin ninguna necesidad.

### 5.2 Pivoteo parcial vs. pivoteo escalado

En el paso 1, el pivoteo parcial elige F6 porque tiene el número más grande de la columna. El problema es que F6 tiene números grandes en *todas* sus columnas. Si lo comparamos con el resto de su propia fila, su número de la columna 1 no es nada especial (0.38 de su máximo).

El escalado, en cambio, elige F8, porque comparada con su propia fila es la mejor opción (|a(8,1)|/s(8) = 1).

Hay algo que al principio parece raro: con escalado aparece un multiplicador **enorme** (2.47×10⁹, en el paso 3, cuando se usa F1 para eliminar F6), y aun así casi no hay crecimiento (ρ = 1.78). La explicación es que lo que importa no es solo el tamaño de m, sino m comparado con las escalas de las filas (m·s(k)/s(i)). F1 tiene escala 10⁻⁴ y F6 tiene escala 10⁶, así que al multiplicar, m·a(1,j) queda del mismo tamaño que a(6,j) y no se pierde información.

En resumen: el problema nunca fue la fila 1 en sí, sino que a(1,1) es muy chico **comparado con el resto de su propia fila**. Con escalado, F1 se puede usar como pivote en el paso 3 sin problemas.

Resultado: el escalado gana 1.2 bits respecto al parcial, y la incógnita peor resuelta pasa de 9.7 a 10.1 cifras.

---

## 6. Condicionamiento

| Cantidad | Valor | De dónde sale |
|---|---:|---|
| det(A) | 1.3646×10⁹ (exacto, la matriz no es singular) | `verificacion.determinante_exacto` |
| κ₁(A) estimado (Hager–Higham, con la LU) | 4.4675×10¹³ | nuestro código |
| κ₁(A) | 4.4675×10¹³ | numpy, para comparar (coincide) |
| κ₂(A) | 3.0847×10¹³ | numpy, para comparar |
| κ₂(D⁻¹A), con las filas equilibradas | **8.13×10³** | numpy, para comparar |

**Lo que dice la teoría.** La cota clásica es ‖δx‖/‖x‖ ≲ κ(A)·ε = 3.08×10¹³ · 2.22×10⁻¹⁶ ≈ **6.8×10⁻³** (con κ₁ da 9.9×10⁻³). Según esto, solo podríamos asegurar unas **2 cifras** correctas. Pero en la práctica obtuvimos entre 10 y 13. ¿Por qué tanta diferencia?

**La fila 7 es casi una combinación de otras.** Haciendo un ajuste por mínimos cuadrados (en aritmética exacta) encontramos que

  F7 ≈ 4.17725 · F3 + 583.247 · F5,  con ‖resto‖₂/‖F7‖₂ = 2.28×10⁻³.

El resto, en valor absoluto, es ≈ 1.9×10⁻⁴, que es justo la perturbación "de orden 10⁻⁴" que menciona el enunciado. Esto también aparece en la eliminación escalada: F7 queda como la última fila de U, con u(12,12) = 1.57×10⁻⁴, unas 400 veces más chico que los números normales de F7 (6.3×10⁻²). Ese pivote tan pequeño es la "marca" que deja la casi dependencia.

**Por qué la cota exagera tanto.** De los 13 órdenes de magnitud de κ₂(A), casi 10 vienen solamente de que las filas tienen tamaños distintos. Si dividimos cada fila por su número más grande (equilibrar las filas), κ₂ baja de 3.1×10¹³ a 8.1×10³. Lo que queda (≈ 10³–10⁴) sí es mal condicionamiento de verdad, y lo causa la casi dependencia de F7.

Si el algoritmo respeta los tamaños de cada fila (como el pivoteo escalado), el error depende de ese condicionamiento equilibrado y no del κ(A) original. Con κ₂(D⁻¹A)·ε ≈ 1.8×10⁻¹² esperaríamos **unas 11 o 12 cifras**, y medimos entre 10.1 y 13.1, así que calza bien.

La cota κ·ε se pone en el peor caso posible: supone que el error va en la peor dirección y que es del tamaño de la fila más grande (10⁶). Pero en realidad el error de redondeo de cada fila es proporcional a su propio tamaño. Como dice el enunciado, la cota es para el peor caso y lo que realmente manda es el residuo.

Además, en la corrida escalada el error hacia atrás es ‖r‖/(‖A‖‖x‖+‖b‖) = 1.1×10⁻¹⁷, que es menor que ε. Eso significa que el algoritmo es estable hacia atrás.

---

## 7. Resultados

La solución que entregamos se obtuvo con **pivoteo parcial escalado + 2 vueltas de refinamiento iterativo**. Con la primera vuelta ya se llega al máximo que permite la doble precisión; la segunda sirve para confirmar que ya no cambia nada (su corrección es menor que u). Una tercera vuelta deja x exactamente igual.

Cómo contamos las cifras.

- **Cifras afirmadas:** las calculamos **sin mirar la solución exacta**. La corrección dᵢ que entrega el refinamiento nos dice más o menos cuánto error tiene xᵢ, y a eso le sumamos el error de redondeo u = 1.1×10⁻¹⁶.
- **Cifras verificadas:** comparamos el valor que imprimimos con la solución exacta en fracciones.

Usamos la definición de Burden & Faires: xᵢ tiene t cifras significativas si |x − x*|/|x*| ≤ 5×10⁻ᵗ. Nunca decimos que hay más de 16 cifras, porque la doble precisión no alcanza para 17.

| i | xᵢ reportado | Cifras afirmadas | Cifras verificadas |
|---|---:|:---:|:---:|
| 1 | 7.384915627000000×10⁻¹ | 16 | ≥ 16 |
| 2 | −1.205937416800000 | 16 | ≥ 16 |
| 3 | 4.387195200000000×10⁻³ | 16 | ≥ 16 |
| 4 | −9.172645380000000×10⁻² | 16 | ≥ 16 |
| 5 | 2.846109372500000 | 16 | ≥ 16 |
| 6 | −6.281540900000000×10⁻³ | 16 | ≥ 16 |
| 7 | 1.473820591600000 | 16 | ≥ 16 |
| 8 | −3.849271640000000×10⁻² | 16 | ≥ 16 |
| 9 | 9.152648037000000×10⁻¹ | 16 | ≥ 16 |
| 10 | −2.607139485200000 | 16 | ≥ 16 |
| 11 | 5.187362940000000×10⁻² | 16 | ≥ 16 |
| 12 | −1.084937162500000 | 16 | ≥ 16 |

**Sin refinamiento** (solo con la eliminación escalada), las cifras afirmadas son 13, 12, 10, 12, 13, 10, 13, 12, 13, 13, 12 y 13 para x₁ hasta x₁₂, y todas se cumplen al compararlas con la solución exacta (ver `resultados/solucion.csv`). Las que peor salen son x₃ y x₆ (unas 10 cifras), que justamente son las incógnitas más pequeñas (del orden de 10⁻³).

**Algo curioso que notamos:** la solución exacta coincide, hasta ~10⁻¹⁷, con números de 10 decimales (0.7384915627, −1.2059374168, …). Parece que la solución se inventó primero con 10 cifras y después se calculó b con 20 dígitos para que el sistema la conservara. Si b viniera con solo 8 cifras, ese error de entrada (10⁻⁸) multiplicado por κ arruinaría completamente el resultado.

---

## 8. Verificación

### 8.1 Residuo calculado de forma exacta

Calculamos el residuo r = b − Ax con `fractions.Fraction`, o sea, con fracciones exactas. A y b se leen como **texto** (así b mantiene sus 20 cifras y nunca pasa por `float`), y x se convierte a fracción exacta desde su valor binario. De esta forma, el residuo no tiene ningún redondeo.

| Corrida | ‖r‖∞ exacto | ‖r‖∞ calculado en doble | Error hacia atrás |
|---|---:|---:|---:|
| Sin pivoteo | 1.772×10⁻⁷ | 1.774×10⁻⁷ | 8.3×10⁻¹⁵ |
| Parcial | 5.373×10⁻¹⁰ | 3.492×10⁻¹⁰ | 2.5×10⁻¹⁷ |
| Parcial escalado | 2.291×10⁻¹⁰ | 5.821×10⁻¹⁰ | 1.1×10⁻¹⁷ |

Fíjense que, en las corridas con pivoteo, el residuo calculado en doble precisión sale mal: se equivoca hasta 2.5 veces. Esto pasa porque al restar b − Ax en doble también se pierden cifras (cancelación). Por eso es obligatorio calcularlo con más precisión.

### 8.2 Refinamiento iterativo

Usamos el residuo exacto y siempre la misma L y U:

| Vuelta | ‖r‖∞ antes de corregir | ‖d‖∞/‖x‖∞ |
|---|---:|---:|
| 1 | 2.291×10⁻¹⁰ | 2.36×10⁻¹³ |
| 2 | 2.603×10⁻¹⁰ | 6.0×10⁻¹⁷ |
| 3 | 2.603×10⁻¹⁰ | 6.0×10⁻¹⁷ (x ya no cambia) |
| final | 2.603×10⁻¹⁰ | — |

Con una sola vuelta, las cifras correctas pasan de 10–13 a más de 16 en las tres corridas. Incluso la corrida sin pivoteo, que partía con 7 cifras, llega a más de 16.

Puede parecer raro que el residuo **no baje** (se queda en 2.6×10⁻¹⁰), pero es normal. Como x está guardado en doble precisión, hay un "piso" de residuo que no se puede bajar: ‖A‖∞‖x‖∞·u ≈ 7.2×10⁶ · 2.85 · 1.1×10⁻¹⁶ ≈ 2.3×10⁻⁹. Nuestro residuo ya está por debajo de eso. Lo que mejora con el refinamiento es el **error**, no el residuo. Cuando la corrección relativa llega a 6×10⁻¹⁷ (menor que u), significa que ya no hay nada más que mejorar.

### 8.3 Comparación con numpy

`numpy.linalg.solve` (que usa LAPACK con pivoteo parcial) da un resultado que difiere del nuestro (escalado) en 7.3×10⁻¹² relativo, y logra entre 10.0 y 13.1 cifras. Es prácticamente lo mismo que nuestro código sin refinar. Solo lo usamos para comparar.

### 8.4 ¿Quedó bien triangular la matriz?

En `traza_eliminacion.txt` está la matriz aumentada después de cada uno de los 11 pasos, tanto para la corrida escalada como para la sin pivoteo. En las dos, al final se comprueba que "66 de 66 elementos bajo la diagonal = 0.0" y que no hubo ningún pivote igual a cero.

---

## 9. Conclusiones

- **Lo que elegimos: pivoteo parcial escalado + refinamiento iterativo con residuo exacto.** La eliminación cuesta 1354 flops con cualquier estrategia. El escalado agrega 283 operaciones (≈ 21 %, que parece mucho solo porque n es chico), pero a cambio entrega la mejor solución con una sola factorización y el menor error hacia atrás. Una vuelta de refinamiento cuesta 276 flops (≈ 20 %) y lleva la solución de 10–13 cifras a más de 16, que es lo máximo posible en doble precisión. Aun sumando todo, sale bastante más barato que resolver usando la inversa.
- **Lo más barato no conviene.** No pivotear no cuesta nada extra, pero se pierden unos 21 bits de mantisa. Pivotear cuesta O(n²) y protege la parte O(n³) del cálculo, así que no vale la pena ahorrárselo.
- **κ(A) no sirvió para predecir el error.** κ(A) ≈ 3×10¹³ es enorme casi solo porque las filas tienen tamaños muy distintos. Lo que de verdad explica lo que pasó es el condicionamiento con las filas equilibradas (≈ 8×10³) y el residuo exacto.
- **¿Qué haríamos distinto si n = 1000?**
  - (a) La eliminación costaría unos 2n³/3 ≈ 6.7×10⁸ flops. Usaríamos una LU **por bloques** (BLAS-3/LAPACK) para aprovechar mejor la memoria caché. El número de operaciones es el mismo, pero con matrices tan grandes lo que más demora es mover los datos en memoria.
  - (b) En vez de hacer pivoteo escalado con divisiones en cada paso, **equilibraríamos las filas una sola vez** al principio (D⁻¹A, cuesta O(n²)) y después usaríamos pivoteo parcial normal. Con n grande, el costo del pivoteo (O(n²)) es casi nada comparado con la eliminación (O(n³)).
  - (c) Calcular el residuo con `Fraction` sería carísimo, porque las fracciones se vuelven enormes. Lo calcularíamos en doble-doble o en `float128` y usaríamos **refinamiento en precisión mixta**.
  - (d) Seguiríamos estimando κ con Hager–Higham, que solo necesita unas pocas resoluciones O(n²), y nunca usando la inversa.

---

## 10. Reproducibilidad

| Ítem | Valor |
|---|---|
| Python | 3.12.10 (funciona con 3.10 o superior) |
| Librerías | numpy 1.26.4 (solo para comparar); `fractions`, `math` y `csv`, que vienen con Python |
| Sistema operativo | Windows 11 |
| Cómo correrlo | `pip install -r requirements.txt` y después `python src/main.py` desde la carpeta principal del repositorio |
| Archivos que genera | `resultados/solucion.csv`, `traza_eliminacion.txt`, `conteo_flops.md`, `residuo.txt` |
| Pruebas | en la rama `pruebas`: `python pruebas/probar_gauss.py` |

Comprobamos que `python src/main.py` funciona en un entorno virtual nuevo, instalando solo lo que está en `requirements.txt`.

---

### Uso de asistentes de IA

 Este trabajo se desarrolló con apoyo de **Claude Code (Anthropic)**, que se usó para: verificar los cálculos y la organización de este informe.

### Fuentes

- Burden, R. L., Faires, J. D. y Burden, A. M. *Análisis Numérico*, 10.ª ed. Cengage, 2017. Cap. 1 (cifras significativas) y cap. 6 (eliminación gaussiana, pivoteo parcial y parcial escalado).
- Higham, N. J. *Accuracy and Stability of Numerical Algorithms*, 2.ª ed. SIAM, 2002. Cap. 9 (factor de crecimiento), cap. 12 (refinamiento iterativo) y cap. 15 (estimación de condición).
- Hager, W. W. "Condition estimates". *SIAM J. Sci. Stat. Comput.* 5(2), 1984, pp. 311–316.
- Higham, N. J. "FORTRAN codes for estimating the one-norm of a real or complex matrix, with applications to condition estimation". *ACM Trans. Math. Softw.* 14(4), 1988, pp. 381–396.
- Golub, G. H. y Van Loan, C. F. *Matrix Computations*, 4.ª ed. Johns Hopkins, 2013. Cap. 3.
- IEEE Std 754-2019, *Standard for Floating-Point Arithmetic*.
