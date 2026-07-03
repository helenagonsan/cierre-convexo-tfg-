# -*- coding: utf-8 -*-
"""
Visualización paso a paso del algoritmo Graham Scan.

Trabajo de Fin de Grado: "El cierre convexo: Teoria y aplicaciones"
Autora: Helena González Sánchez — Universidad de Salamanca

Script autónomo: genera la figura del anexo de visualización con un panel por
paso (orden angular desde el ancla, barrido con pila mostrando el triplete
evaluado, los puntos extraídos y la envolvente final).

Ejecución:  python3 viz_graham.py
Requisitos: numpy, matplotlib
"""

import functools
import numpy as np
import matplotlib.pyplot as plt

plt.ioff()


# ---------------------------------------------------------------------
# Primitiva geométrica (idéntica a la del anexo de implementación)
# ---------------------------------------------------------------------

def orientacion(a, b, c):
    """
    Doble del área con signo del triángulo (a, b, c).
      > 0 : c queda a la izquierda de la recta dirigida a -> b
      < 0 : c queda a la derecha
      = 0 : a, b, c colineales
    """
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])


# ---------------------------------------------------------------------
# Graham Scan instrumentado: registra orden angular y barrido
# ---------------------------------------------------------------------

def graham_scan_pasos(P):
    P = np.asarray(P, float)
    n = len(P)
    p0 = min(range(n), key=lambda i: (P[i][1], -P[i][0]))

    def comparar(i, j):
        o = orientacion(P[p0], P[i], P[j])
        if o > 0:
            return -1
        if o < 0:
            return 1
        di = (P[i][0]-P[p0][0])**2 + (P[i][1]-P[p0][1])**2
        dj = (P[j][0]-P[p0][0])**2 + (P[j][1]-P[p0][1])**2
        return -1 if di < dj else 1

    resto = sorted([i for i in range(n) if i != p0],
                   key=functools.cmp_to_key(comparar))

    pasos = [{"tipo": "orden", "p0": p0, "orden": list(resto)}]
    pila = [p0, resto[0]]
    pasos.append({"tipo": "init", "p0": p0, "pila": list(pila),
                  "procesado": resto[0], "extraidos": [],
                  "tope": pila[-1], "segundo": pila[-2] if len(pila) >= 2 else None})

    for i in resto[1:]:
        extraidos = []
        while len(pila) >= 2 and orientacion(P[pila[-2]], P[pila[-1]], P[i]) <= 0:
            extraidos.append(pila.pop())
        tope = pila[-1]
        segundo = pila[-2] if len(pila) >= 2 else None
        pila.append(i)
        pasos.append({"tipo": "barrido", "p0": p0, "pila": list(pila),
                      "procesado": i, "extraidos": extraidos,
                      "tope": tope, "segundo": segundo})
    return pila, pasos


# ---------------------------------------------------------------------
# Dibujo paso a paso
# ---------------------------------------------------------------------

def dibujar_graham_pasos(P, pasos, hull, ncols=3, archivo="graham_pasos.pdf"):
    plt.close('all')
    P = np.asarray(P, float)
    COL_PROC, COL_PILA = "#dc2626", "#1f2937"
    COL_RAYO, COL_DESC, COL_P0 = "#0f766e", "#9ca3af", "#2563eb"
    COL_REJ = "#f0a0a0"          # triplete rechazado (rojo tenue)
    S_PT = 26
    orden_vertices = {idx: k for k, idx in enumerate(hull)}

    xmin, xmax = P[:, 0].min(), P[:, 0].max()
    ymin, ymax = P[:, 1].min(), P[:, 1].max()
    mx, my = 0.15 * (xmax - xmin), 0.15 * (ymax - ymin)
    xlim, ylim = (xmin - mx, xmax + mx), (ymin - my, ymax + my)

    npan = len(pasos) + 1
    nrows = (npan + ncols - 1) // ncols
    fig, axs = plt.subplots(nrows, ncols, figsize=(3.8 * ncols, 3.8 * nrows))
    axs = np.atleast_1d(axs).ravel()

    def edge(ax, i, j, color, lw, z, ls="-"):
        ax.plot([P[i, 0], P[j, 0]], [P[i, 1], P[j, 1]], color=color, lw=lw, ls=ls, zorder=z)

    def etq(ax, i, txt, color, dx=6, dy=-10, fs=11, bold=False):
        ax.annotate(txt, P[i], textcoords="offset points", xytext=(dx, dy),
                    fontsize=fs, color=color, annotation_clip=False, zorder=9,
                    fontweight=("bold" if bold else "normal"))

    for k, paso in enumerate(pasos):
        ax = axs[k]
        p0 = paso["p0"]
        ax.scatter(P[:, 0], P[:, 1], color=COL_PILA, s=S_PT-6, zorder=1)
        ax.scatter([P[p0, 0]], [P[p0, 1]], color=COL_P0, s=S_PT+24, zorder=6)
        etq(ax, p0, "$p_0$", COL_P0, dx=-20, dy=-12, bold=True)

        if paso["tipo"] == "orden":
            for rango, idx in enumerate(paso["orden"], start=1):
                edge(ax, p0, idx, COL_RAYO, 0.7, 2, ls="--")
                etq(ax, idx, str(rango), COL_PILA, dx=5, dy=4, fs=8)
            ax.set_title(f"Paso {k+1}: orden angular desde $p_0$", fontsize=9)

        else:
            pila = paso["pila"]; proc = paso["procesado"]; extr = paso["extraidos"]
            tope, segundo = paso["tope"], paso["segundo"]

            # Cadena actual de la pila.
            for t in range(len(pila) - 1):
                edge(ax, pila[t], pila[t+1], COL_PILA, 2.2, 5)

            # Rayo angular a p0 (indica el orden, NO el predicado).
            edge(ax, p0, proc, COL_RAYO, 1.2, 3, ls="--")

            # Triplete RECHAZADO: p_t -> (último extraído) -> p_i  (giro a la derecha/colineal).
            hay_leyenda = False
            if extr:
                e_last = extr[-1]
                ax.plot([P[tope, 0], P[e_last, 0]], [P[tope, 1], P[e_last, 1]],
                        color=COL_REJ, lw=1.8, ls=(0, (4, 2)), zorder=4)
                ax.plot([P[e_last, 0], P[proc, 0]], [P[e_last, 1], P[proc, 1]],
                        color=COL_REJ, lw=1.8, ls=(0, (4, 2)), zorder=4,
                        label="triplete rechazado")
                ax.scatter(P[extr, 0], P[extr, 1], color=COL_DESC, marker="x",
                           s=S_PT+10, zorder=5, label="extraídos")
                hay_leyenda = True

            # Triplete ACEPTADO: p_{t-1} -> p_t -> p_i  (giro a la izquierda).
            if segundo is not None:
                edge(ax, segundo, tope, COL_PROC, 2.6, 6)
                edge(ax, tope, proc, COL_PROC, 2.6, 6)

            ax.scatter(P[pila, 0], P[pila, 1], color=COL_PILA, s=S_PT+6, zorder=7)
            ax.scatter([P[proc, 0]], [P[proc, 1]], color=COL_PROC, s=S_PT+30, zorder=8)

            # Etiquetas de los tres puntos del predicado.
            etq(ax, proc, "$p_i$", COL_PROC, dx=6, dy=6, bold=True)
            etq(ax, tope, "$p_t$", COL_PILA, dx=6, dy=-12, bold=True)
            if segundo is not None and segundo != p0:
                etq(ax, segundo, "$p_{t-1}$", COL_PILA, dx=-30, dy=-6, bold=True)

            if hay_leyenda:
                ax.legend(fontsize=7, loc="upper left", framealpha=0.9)

            tipo = "inicial" if paso["tipo"] == "init" else "barrido"
            ax.set_title(rf"Paso {k+1}: {tipo}, evalúa $D(p_{{t-1}}, p_t, p_i)$", fontsize=9)

        ax.set_xlim(*xlim); ax.set_ylim(*ylim)
        ax.set_aspect("equal"); ax.axis("off")

    # Panel final: CH cerrada.
    ax = axs[len(pasos)]
    ax.scatter(P[:, 0], P[:, 1], color=COL_PILA, s=S_PT-6, zorder=1)
    ciclo = hull + [hull[0]]
    ax.fill(P[ciclo, 0], P[ciclo, 1], color=COL_PROC, alpha=0.07, zorder=1)
    for t in range(len(hull)):
        edge(ax, hull[t], hull[(t+1) % len(hull)], COL_PILA, 2.2, 5)
    ax.scatter(P[hull, 0], P[hull, 1], color=COL_PROC, s=S_PT+18, zorder=7)
    for idx in hull:
        etq(ax, idx, f"$p_{{{orden_vertices[idx]}}}$", COL_PILA, dx=8, dy=-8, fs=10, bold=True)
    ax.set_title(f"Paso {len(pasos)+1}: envolvente final", fontsize=9)
    ax.set_xlim(*xlim); ax.set_ylim(*ylim)
    ax.set_aspect("equal"); ax.axis("off")

    for k in range(npan, len(axs)):
        axs[k].axis("off")

    fig.suptitle("Ejecución de Graham Scan", fontsize=14)
    plt.tight_layout(rect=[0, 0, 1, 0.97])
    plt.savefig(archivo, bbox_inches="tight")
    plt.show()


if __name__ == "__main__":
    P = np.array([
        [4, 0],                    # p0 (y mínima)
        [7, 1], [9, 4], [8, 7],    # cadena derecha
        [6, 9], [3, 8],            # cadena superior
        [1, 5], [2, 2],            # cadena izquierda
        [5, 3], [4, 5], [6, 5],    # interiores (provocan extracciones)
    ], dtype=float)

    hull, pasos = graham_scan_pasos(P)
    print("Vértices de la envolvente (antihorario):", hull)
    dibujar_graham_pasos(P, pasos, hull, archivo="graham_pasos.pdf")
    print("Figura guardada en graham_pasos.pdf")
