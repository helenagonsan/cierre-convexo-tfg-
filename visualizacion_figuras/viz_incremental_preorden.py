# -*- coding: utf-8 -*-
"""
Visualización paso a paso del algoritmo Incremental, variante con
preordenamiento por la coordenada x.

Trabajo de Fin de Grado: "El cierre convexo: Teoria y aplicaciones"
Autora: Helena González Sánchez — Universidad de Salamanca

Script autónomo. Ilustra la variante optimizada del algoritmo incremental
descrita en el análisis de complejidad: los puntos se preordenan por x antes de
incorporarse uno a uno, lo que admite resolución mediante búsqueda binaria y
coste O(n log n). La figura complementa a la de la variante directa (orden de
entrada), que ilustra la versión de coste O(n^2).

Nota sobre la inicialización: con el fin de no recargar la instrumentación, este
script de figuras emplea la inicialización directa con los tres primeros puntos
del orden por x, válida bajo la hipótesis de posición general. La versión robusta
de la inicialización es la del anexo de implementación; ambas coinciden cuando
esos tres puntos no están alineados, como en el ejemplo.

Ejecución:  python3 viz_incremental_preorden.py
Requisitos: numpy, matplotlib
"""

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
# Incremental instrumentado, con preordenamiento por x
# ---------------------------------------------------------------------

def incremental_pasos(P):
    P = np.asarray(P, float)
    n = len(P)
    orden = sorted(range(n), key=lambda i: (P[i][0], P[i][1]))   # preorden por x (desempate y)

    a, b, c = orden[0], orden[1], orden[2]
    H = [a, b, c] if orientacion(P[a], P[b], P[c]) > 0 else [a, c, b]

    pasos = [{"tipo": "init", "H": list(H)}]
    for k in orden[3:]:
        m = len(H)
        visibles = [i for i in range(m)
                    if orientacion(P[H[i]], P[H[(i+1) % m]], P[k]) < 0]
        interior = (len(visibles) == 0)

        H_antes = list(H)
        tangencias = []
        if not interior:
            for i in range(m):
                ant = orientacion(P[H[(i-1) % m]], P[H[i]], P[k])
                sig = orientacion(P[H[i]], P[H[(i+1) % m]], P[k])
                if (ant < 0) != (sig < 0):
                    tangencias.append(H[i])
            nuevo = []
            for i in range(m):
                ant = orientacion(P[H[(i-1) % m]], P[H[i]], P[k])
                sig = orientacion(P[H[i]], P[H[(i+1) % m]], P[k])
                if ant >= 0:
                    nuevo.append(H[i])
                if ant < 0 and sig >= 0:
                    nuevo.append(k)
                    nuevo.append(H[i])
            H = nuevo

        removidos = [v for v in H_antes if v not in H]   # vértices eliminados en este paso
        pasos.append({"tipo": "iter", "k": k, "interior": interior,
                      "H_antes": H_antes, "H_despues": list(H),
                      "visibles": [(H_antes[i], H_antes[(i+1) % m]) for i in visibles],
                      "tangencias": tangencias, "removidos": removidos})
    return H, pasos


# ---------------------------------------------------------------------
# Dibujo paso a paso
# ---------------------------------------------------------------------

def dibujar_incremental_pasos(P, pasos, hull, ncols=3, archivo="incremental_preorden.pdf"):
    plt.close('all')
    P = np.asarray(P, float)
    COL_PROC, COL_HULL = "#dc2626", "#1f2937"
    COL_VIS, COL_TAN, COL_NUEVA = "#dc2626", "#2563eb", "#16a34a"
    COL_ELIM = "#d97706"
    S_PT = 26
    orden_vertices = {idx: k for k, idx in enumerate(hull)}

    xmin, xmax = P[:, 0].min(), P[:, 0].max()
    ymin, ymax = P[:, 1].min(), P[:, 1].max()
    mx, my = 0.15 * (xmax - xmin), 0.15 * (ymax - ymin)
    xlim, ylim = (xmin - mx, xmax + mx), (ymin - my, ymax + my)

    npan = len(pasos) + 1
    nrows = (npan + ncols - 1) // ncols
    fig, axs = plt.subplots(nrows, ncols, figsize=(3.6 * ncols, 3.6 * nrows))
    axs = np.atleast_1d(axs).ravel()

    def cadena(ax, H, color, lw, z, ls="-"):
        m = len(H)
        for i in range(m):
            a, b = H[i], H[(i+1) % m]
            ax.plot([P[a, 0], P[b, 0]], [P[a, 1], P[b, 1]],
                    color=color, lw=lw, ls=ls, zorder=z)

    def etq(ax, i, txt, color, dx=6, dy=-10, fs=11, bold=False):
        ax.annotate(txt, P[i], textcoords="offset points", xytext=(dx, dy),
                    fontsize=fs, color=color, annotation_clip=False, zorder=9,
                    fontweight=("bold" if bold else "normal"))

    for k_idx, paso in enumerate(pasos):
        ax = axs[k_idx]
        ax.scatter(P[:, 0], P[:, 1], color=COL_HULL, s=S_PT-6, zorder=1)

        if paso["tipo"] == "init":
            H = paso["H"]
            ax.fill(P[H + [H[0]], 0], P[H + [H[0]], 1], color=COL_PROC, alpha=0.06, zorder=1)
            cadena(ax, H, COL_HULL, 2.2, 5)
            ax.scatter(P[H, 0], P[H, 1], color=COL_HULL, s=S_PT+6, zorder=7)
            ax.set_title(f"Paso {k_idx+1}: triángulo inicial $H_2$", fontsize=9)
        else:
            k = paso["k"]; H_antes = paso["H_antes"]
            cadena(ax, H_antes, COL_HULL, 2.0, 4)
            if paso["interior"]:
                ax.scatter([P[k, 0]], [P[k, 1]], color="#9ca3af", marker="x",
                           s=S_PT+30, zorder=8, label="descartado (interior)")
                ax.legend(fontsize=7, loc="upper left", framealpha=0.9)
                ax.set_title(f"Paso {k_idx+1}: $p_k$ interior, se descarta", fontsize=9)
            else:
                for (a, b) in paso["visibles"]:
                    ax.plot([P[a, 0], P[b, 0]], [P[a, 1], P[b, 1]],
                            color=COL_VIS, lw=2.6, zorder=6)
                for j, tg in enumerate(paso["tangencias"]):
                    ax.plot([P[tg, 0], P[k, 0]], [P[tg, 1], P[k, 1]],
                            color=COL_NUEVA, lw=2.2, zorder=6)
                    ax.scatter([P[tg, 0]], [P[tg, 1]], color=COL_TAN, s=S_PT+18, zorder=7,
                               label="punto de tangencia" if j == 0 else None)
                if paso["removidos"]:
                    ax.scatter(P[paso["removidos"], 0], P[paso["removidos"], 1],
                               color=COL_ELIM, marker="x", s=S_PT+24, zorder=8,
                               label="eliminado de H")
                ax.scatter([P[k, 0]], [P[k, 1]], color=COL_PROC, s=S_PT+30, zorder=8)
                etq(ax, k, "$p_k$", COL_PROC, dx=6, dy=6, bold=True)
                if paso["tangencias"] or paso["removidos"]:
                    ax.legend(fontsize=7, loc="upper left", framealpha=0.9)
                titulo_extra = " (elimina vértices)" if paso["removidos"] else ""
                ax.set_title(f"Paso {k_idx+1}: $p_k$ exterior{titulo_extra}", fontsize=9)

        ax.set_xlim(*xlim); ax.set_ylim(*ylim)
        ax.set_aspect("equal"); ax.axis("off")

    ax = axs[len(pasos)]
    ax.scatter(P[:, 0], P[:, 1], color=COL_HULL, s=S_PT-6, zorder=1)
    ax.fill(P[hull + [hull[0]], 0], P[hull + [hull[0]], 1], color=COL_PROC, alpha=0.07, zorder=1)
    cadena(ax, hull, COL_HULL, 2.2, 5)
    ax.scatter(P[hull, 0], P[hull, 1], color=COL_PROC, s=S_PT+18, zorder=7)
    for idx in hull:
        etq(ax, idx, f"$p_{{{orden_vertices[idx]}}}$", COL_HULL, dx=8, dy=-8, fs=10, bold=True)
    ax.set_title(f"Paso {len(pasos)+1}: envolvente final", fontsize=9)
    ax.set_xlim(*xlim); ax.set_ylim(*ylim)
    ax.set_aspect("equal"); ax.axis("off")

    for k_idx in range(npan, len(axs)):
        axs[k_idx].axis("off")

    fig.suptitle("Ejecución del algoritmo Incremental (preorden por $x$)", fontsize=14)
    plt.tight_layout(rect=[0, 0, 1, 0.97])
    plt.savefig(archivo, bbox_inches="tight")
    plt.show()


if __name__ == "__main__":
    P = np.array([
        [0, 3], [1, 6], [2, 1],     # tres de menor x -> triángulo inicial
        [3, 7], [4, 0], [5, 2],
        [5.5, 4], [6, 8], [7, -1], [8, 5],
    ], dtype=float)

    hull, pasos = incremental_pasos(P)
    print("Vértices de la envolvente (antihorario):", hull)
    dibujar_incremental_pasos(P, pasos, hull, ncols=3, archivo="incremental_preorden.pdf")
    print("Figura guardada en incremental_preorden.pdf")
