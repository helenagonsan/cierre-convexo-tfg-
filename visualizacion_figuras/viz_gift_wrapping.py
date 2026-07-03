# -*- coding: utf-8 -*-
"""
Visualización paso a paso del algoritmo Gift Wrapping (marcha de Jarvis).

Trabajo de Fin de Grado: "El cierre convexo: Teoria y aplicaciones"
Autora: Helena González Sánchez — Universidad de Salamanca

Script autónomo: genera la figura del anexo de visualización con una subfigura
por iteración (recta dirigida de referencia, segmentos candidatos, arista
elegida y envolvente parcial).

Ejecución:  python3 viz_gift_wrapping.py
Requisitos: numpy, matplotlib
"""

import numpy as np
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------
# Primitivas geométricas (idénticas a las del anexo de implementación)
# ---------------------------------------------------------------------

def orientacion(a, b, c):
    """
    Doble del área con signo del triángulo (a, b, c).
      > 0 : c queda a la izquierda de la recta dirigida a -> b
      < 0 : c queda a la derecha
      = 0 : a, b, c colineales
    """
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])


def distancia_sq(p, q):
    """Distancia euclídea al cuadrado (evita la raíz; suficiente para comparar)."""
    return (p[0] - q[0])**2 + (p[1] - q[1])**2


# ---------------------------------------------------------------------
# Gift Wrapping instrumentado: registra el estado de cada iteración
# ---------------------------------------------------------------------

def jarvis_pasos(P):
    n = len(P)
    inicio = min(range(n), key=lambda i: (P[i][1], -P[i][0]))   # mínima y, desempate máxima x
    hull, pasos = [], []
    actual = inicio
    while True:
        hull.append(actual)
        siguiente = (actual + 1) % n
        for i in range(n):
            if i == actual:
                continue
            o = orientacion(P[actual], P[siguiente], P[i])
            if o < 0 or (o == 0 and distancia_sq(P[actual], P[i]) > distancia_sq(P[actual], P[siguiente])):
                siguiente = i

        # Dirección de referencia d_i (solo para la figura):
        #   i = 0  -> horizontal hacia la derecha
        #   i >= 1 -> dirección de la arista de entrada  p_{i-1} -> p_i
        k = len(hull) - 1
        if k == 0:
            d = np.array([1.0, 0.0])
        else:
            v = np.asarray(P[actual], float) - np.asarray(P[hull[k - 1]], float)
            d = v / np.linalg.norm(v)

        pasos.append({
            "i": k, "actual": actual, "siguiente": siguiente,
            "hull_parcial": hull.copy(), "inicio": inicio, "direccion": d,
        })
        actual = siguiente
        if actual == inicio:
            break
    return hull, pasos


# ---------------------------------------------------------------------
# Dibujo de la rejilla de iteraciones
# ---------------------------------------------------------------------

def dibujar_pasos(P, pasos, ncols=3, archivo="jarvis_pasos.pdf"):
    P = np.asarray(P, float)

    # ----- estilo -----
    COL_PUNTO   = "#1f2937"
    COL_CAND    = "#cbd5e1"
    COL_ELEGIDA = "#dc2626"
    COL_DIR     = "#0f766e"
    LW_HULL, LW_ELEGIDA, LW_CAND, LW_DIR = 1.8, 2.4, 0.7, 1.6
    S_PUNTO = 24

    ancho = P[:, 0].max() - P[:, 0].min()
    alto  = P[:, 1].max() - P[:, 1].min()
    L = 0.20 * np.hypot(ancho, alto)
    # El margen debe contener la flecha completa, no solo los puntos.
    margen = max(0.18 * max(ancho, alto), 1.2 * L)

    npasos = len(pasos)
    nrows = (npasos + ncols - 1) // ncols
    fig, axs = plt.subplots(nrows, ncols, figsize=(3.6 * ncols, 3.6 * nrows))
    axs = np.atleast_1d(axs).ravel()

    for k, paso in enumerate(pasos):
        ax = axs[k]
        i = paso["i"]
        actual, siguiente = paso["actual"], paso["siguiente"]
        hull_parcial, inicio, d = paso["hull_parcial"], paso["inicio"], paso["direccion"]
        pa = P[actual]

        for j in range(len(P)):
            if j != actual and j not in hull_parcial:
                ax.plot([pa[0], P[j, 0]], [pa[1], P[j, 1]],
                        ls="--", color=COL_CAND, lw=LW_CAND, zorder=1)

        for t in range(len(hull_parcial) - 1):
            a, b = hull_parcial[t], hull_parcial[t + 1]
            ax.plot([P[a, 0], P[b, 0]], [P[a, 1], P[b, 1]],
                    color=COL_PUNTO, lw=LW_HULL, zorder=2)

        # dirección d_i (flecha + etiqueta). annotation_clip=False evita recortes.
        tip = pa + L * d
        ax.annotate("", xy=tip, xytext=pa,
                    arrowprops=dict(arrowstyle="-|>", color=COL_DIR, lw=LW_DIR),
                    zorder=5, annotation_clip=False)
        ax.annotate(rf"$\vec{{d}}_{{{i}}}$", xy=tip, textcoords="offset points",
                    xytext=(4, 4), fontsize=8, color=COL_DIR, zorder=6,
                    annotation_clip=False)

        ps = P[siguiente]
        ax.plot([pa[0], ps[0]], [pa[1], ps[1]],
                color=COL_ELEGIDA, lw=LW_ELEGIDA, zorder=4)

        ax.scatter(P[:, 0], P[:, 1], color=COL_PUNTO, s=S_PUNTO, zorder=3)

        ax.annotate(r"$p_0$", P[inicio], textcoords="offset points",
                    xytext=(-18, -6), fontsize=11, color=COL_PUNTO, zorder=6)
        if actual != inicio:
            ax.annotate(rf"$p_{{{i}}}$", pa, textcoords="offset points",
                        xytext=(-18, -6), fontsize=11, color=COL_PUNTO, zorder=6)

        ax.set_title(f"Iteración {k + 1}", fontsize=10)
        ax.set_aspect("equal")
        ax.set_xlim(P[:, 0].min() - margen, P[:, 0].max() + margen)
        ax.set_ylim(P[:, 1].min() - margen, P[:, 1].max() + margen)
        ax.axis("off")

    for k in range(npasos, len(axs)):
        axs[k].axis("off")

    fig.suptitle("Ejecución de Gift Wrapping", fontsize=14)
    plt.tight_layout()
    plt.savefig(archivo, bbox_inches="tight")
    plt.show()


if __name__ == "__main__":
    P = np.array([
        [4, 0],   # vértice más bajo -> punto de partida
        [8, 2], [8, 6], [4, 8], [0, 6], [1, 2],          # resto de vértices de la envolvente
        [3, 3], [5, 5], [5, 3], [3, 5], [4, 4], [6, 4],  # puntos interiores
    ], dtype=float)
    hull, pasos = jarvis_pasos(P)
    print("Vértices de CH(P):", hull)
    dibujar_pasos(P, pasos, archivo="jarvis_pasos.pdf")
    print("Figura guardada en jarvis_pasos.pdf")
