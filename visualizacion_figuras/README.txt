========================================================================
Material de visualización — Anexo de figuras
Trabajo de Fin de Grado: "El cierre convexo: Teoría y aplicaciones"
Autora: Helena González Sánchez — Universidad de Salamanca
========================================================================

CONTENIDO

Scripts de visualización paso a paso de los algoritmos de envolvente convexa
en el plano. Cada script genera la figura del anexo correspondiente y, al
ejecutarlo, muestra la figura en pantalla y la guarda como PDF.

  viz_gift_wrapping.py        Gift Wrapping (marcha de Jarvis)
  viz_quickhull.py            QuickHull
  viz_graham.py               Graham Scan
  viz_incremental_directo.py  Incremental, variante directa (orden de entrada)
  viz_incremental_preorden.py Incremental, variante con preordenamiento por x
  viz_divide_venceras.py      Divide y Vencerás (fase de fusión: tangentes comunes)

PDFs que genera cada script (en la carpeta desde la que se ejecuta):
  jarvis_pasos.pdf, quickhull_pasos.pdf, graham_pasos.pdf,
  incremental_directo.pdf, incremental_preorden.pdf, dc_fusion.pdf


REQUISITOS

  - Python 3
  - numpy
  - matplotlib

  Instalación de las dependencias (Windows):
      py -m pip install numpy matplotlib
  En macOS o Linux:
      pip3 install numpy matplotlib


EJECUCIÓN

  Cada script es autónomo (incluye las primitivas geométricas que necesita) y
  se ejecuta de forma independiente, por ejemplo:

      python3 viz_quickhull.py

  Al ejecutarlo se abre una ventana con la figura y se guarda el PDF
  correspondiente en la misma carpeta.


NOTAS DE COHERENCIA CON LA MEMORIA

  - Las primitivas geométricas (orientacion, distancia_sq) de estos scripts son
    idénticas a las del anexo de implementación de los algoritmos.

  - Las dos figuras del algoritmo Incremental ilustran las dos variantes
    descritas en el análisis de complejidad de la memoria: la directa, de coste
    O(n^2), y la optimizada con preordenamiento por x, de coste O(n log n). El
    anexo de implementación contiene una única versión del algoritmo (la robusta
    sin preordenar); ambas producen la misma envolvente.

  - Por simplicidad de la instrumentación, los scripts del Incremental usan la
    inicialización directa con los tres primeros puntos, válida en posición
    general. La inicialización robusta completa (búsqueda del primer triple no
    colineal) es la del anexo de implementación; ambas coinciden cuando esos
    tres puntos no están alineados, como en los ejemplos.

  - El script de Divide y Vencerás ilustra la fase de fusión (cálculo de las
    tangentes comunes inferior y superior), que es el núcleo del algoritmo.
