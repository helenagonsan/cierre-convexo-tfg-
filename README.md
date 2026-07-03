# El cierre convexo: Teoría y aplicaciones — Anexos digitales

Código que acompaña al Trabajo de Fin de Grado **«El cierre convexo: Teoría y
aplicaciones»**.

- **Autora:** Helena González Sánchez
- **Tutor:** Ángel María Martín del Rey
- **Grado en Matemáticas — Universidad de Salamanca**
- **Curso académico 2025–2026**

## Descripción

Este repositorio reúne los anexos digitales de la memoria: la implementación de
los cinco algoritmos clásicos para el cálculo de la envolvente convexa en el
plano y los scripts que generan las figuras de visualización paso a paso.

Todos los algoritmos se construyen sobre un único predicado geométrico
elemental —el predicado de orientación, evaluado como el signo de un
determinante y exacto en aritmética entera—, de modo que el resultado no depende
de errores de redondeo en coma flotante.

## Estructura del repositorio

```
.
├── implementacion_computacional/
│   └── envolvente_convexa.py        Los cinco algoritmos (biblioteca estándar)
├── visualizacion_figuras/
│   ├── README.txt                   Instrucciones específicas de las figuras
│   ├── viz_gift_wrapping.py
│   ├── viz_quickhull.py
│   ├── viz_graham.py
│   ├── viz_incremental_directo.py
│   ├── viz_incremental_preorden.py
│   └── viz_divide_venceras.py
├── requirements.txt
└── README.md
```

## Requisitos

- **Anexo de implementación** (`implementacion_computacional/`): Python 3, sin
  dependencias externas. Se ejecuta con la biblioteca estándar.
- **Anexo de visualización** (`visualizacion_figuras/`): Python 3, `numpy` y
  `matplotlib`.

Instalación de las dependencias de la parte de visualización:

```bash
pip install -r requirements.txt
```

## Uso

### Algoritmos

Cada función recibe una lista de puntos `P`, donde cada punto es una tupla
`(x, y)`, y devuelve los índices de los vértices de la envolvente convexa en
orden antihorario. El fichero se puede ejecutar directamente:

```bash
python3 implementacion_computacional/envolvente_convexa.py
```

Esto imprime la envolvente de un conjunto de ejemplo con los cinco algoritmos y,
a continuación, ejecuta una verificación de concordancia sobre 2000 conjuntos
aleatorios en posición general.

### Figuras

Cada script de visualización es autónomo, muestra la figura en pantalla y la
guarda como PDF en la carpeta desde la que se ejecuta. Por ejemplo:

```bash
cd visualizacion_figuras
python3 viz_quickhull.py
```

Véase `visualizacion_figuras/README.txt` para el detalle de cada figura.

## Algoritmos incluidos

| Algoritmo           | Complejidad         |
| ------------------- | ------------------- |
| Gift Wrapping (Jarvis) | O(nh)            |
| QuickHull           | O(n log n) esperado |
| Graham Scan         | O(n log n)          |
| Incremental         | O(n²) / O(n log n)  |
| Divide y Vencerás   | O(n log n)          |

(donde *n* es el número de puntos y *h* el número de vértices de la envolvente).

## Validación

El bloque de prueba de `envolvente_convexa.py` comprueba que los cinco
algoritmos producen el mismo conjunto de vértices sobre 2000 conjuntos de puntos
generados aleatoriamente en posición general, tomando Gift Wrapping como
referencia.

## Licencia

Véase el fichero `LICENSE`.
