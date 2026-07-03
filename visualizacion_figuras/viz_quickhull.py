# -*- coding: utf-8 -*-
"""
Visualización paso a paso del algoritmo QuickHull.

Trabajo de Fin de Grado: "El cierre convexo: Teoria y aplicaciones"
Autora: Helena González Sánchez — Universidad de Salamanca

Script autónomo: genera la figura del anexo de visualización con un panel por
paso (división inicial por la recta ab, triángulos de descarte con sus
subconjuntos A y B, y envolvente final).

Ejecución:  python3 viz_quickhull.py
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
# QuickHull instrumentado: registra los eventos de la recursión
# ---------------------------------------------------------------------

def quickhull_pasos(P):
    P = np.asarray(P, float)
    n = len(P)
    eventos = []
    a = max(range(n), key=lambda i: (P[i][0], P[i][1]))
    b = min(range(n), key=lambda i: (P[i][0], -P[i][1]))
    S_inf = [i for i in range(n) if orientacion(P[a], P[b], P[i]) > 0]
    S_sup = [i for i in range(n) if orientacion(P[a], P[b], P[i]) < 0]
    eventos.append({"tipo": "init", "a": a, "b": b,
                    "S_inf": list(S_inf), "S_sup": list(S_sup)})
    cad_inf = _qh_pasos(P, a, b, S_inf, eventos)
    cad_sup = _qh_pasos(P, b, a, S_sup, eventos)
    hull = [a] + cad_inf + [b] + cad_sup
    hull = [hull[0]] + hull[1:][::-1]   # la cadena se construye horaria; se invierte
    return hull, eventos


def _qh_pasos(P, p, q, T, eventos):
    if not T:
        eventos.append({"tipo": "arista", "p": p, "q": q})
        return []
    c = max(T, key=lambda i: orientacion(P[p], P[q], P[i]))
    A = [i for i in T if orientacion(P[p], P[c], P[i]) > 0]
    B = [i for i in T if orientacion(P[c], P[q], P[i]) > 0]
    descartados = [i for i in T if i != c and i not in A and i not in B]
    eventos.append({"tipo": "triangulo", "p": p, "q": q, "c": c,
                    "T": list(T), "A": list(A), "B": list(B),
                    "descartados": descartados})
    return _qh_pasos(P, p, c, A, eventos) + [c] + _qh_pasos(P, c, q, B, eventos)


# ---------------------------------------------------------------------
# Dibujo paso a paso
# ---------------------------------------------------------------------

def dibujar_quickhull_pasos(P, eventos, hull, ncols=3, archivo="quickhull_pasos.pdf"):
    plt.close('all')
    P = np.asarray(P, float)
    COL_BASE, COL_C = "#0f766e", "#dc2626"
    COL_A, COL_B, COL_DESC = "#2563eb", "#d97706", "#9ca3af"
    COL_HULL = "#1f2937"
    S_PT = 26
    orden_vertices = {idx: i for i, idx in enumerate(hull)}

    ev_init = next(ev for ev in eventos if ev["tipo"] == "init")
    global_a, global_b = ev_init["a"], ev_init["b"]

    previas, paneles = set(), []
    for ev in eventos:
        if ev["tipo"] == "arista":
            previas.add(tuple(sorted((ev["p"], ev["q"]))))
        elif ev["tipo"] == "init":
            paneles.append({"tipo": "init", **ev, "previas": sorted(previas)})
        elif ev["tipo"] == "triangulo":
            if not ev["A"]:
                previas.add(tuple(sorted((ev["p"], ev["c"]))))
            if not ev["B"]:
                previas.add(tuple(sorted((ev["c"], ev["q"]))))
            paneles.append({"tipo": "triangulo", **ev, "previas": sorted(previas)})
    paneles.append({"tipo": "final", "previas": sorted(previas)})
    print(f"Mostrando los {len(paneles)} paneles (todos los pasos).")

    xmin, xmax = P[:, 0].min(), P[:, 0].max()
    ymin, ymax = P[:, 1].min(), P[:, 1].max()
    mx, my = 0.15 * (xmax - xmin), 0.15 * (ymax - ymin)
    xlim, ylim = (xmin - mx, xmax + mx), (ymin - my, ymax + my)

    npan = len(paneles)
    nrows = (npan + ncols - 1) // ncols
    fig, axs = plt.subplots(nrows, ncols, figsize=(3.8 * ncols, 3.8 * nrows))
    axs = np.atleast_1d(axs).ravel()

    def edge(ax, i, j, color, lw, z, ls="-"):
        ax.plot([P[i, 0], P[j, 0]], [P[i, 1], P[j, 1]], color=color, lw=lw, ls=ls, zorder=z)
    def etq_hull(ax, i, color, dx=8, dy=-8):
        if i in orden_vertices:
            ax.annotate(f"$p_{{{orden_vertices[i]}}}$", P[i], textcoords="offset points",
                        xytext=(dx, dy), fontsize=10, color=color, annotation_clip=False,
                        zorder=9, fontweight='bold')
    def etq_generica(ax, i, txt, color, dx=-16, dy=-6):
        ax.annotate(txt, P[i], textcoords="offset points", xytext=(dx, dy),
                    fontsize=11, color=color, annotation_clip=False, zorder=9)
    def hull_previo(ax, previas):
        for (i, j) in previas:
            edge(ax, i, j, COL_HULL, 2.2, 5)

    for k, pan in enumerate(paneles):
        ax = axs[k]
        ax.scatter(P[:, 0], P[:, 1], color=COL_HULL, s=S_PT-6, zorder=1)

        if pan["tipo"] == "init":
            a, b = pan["a"], pan["b"]
            edge(ax, a, b, COL_HULL, 1.6, 2)
            if pan["S_inf"]:
                ax.scatter(P[pan["S_inf"], 0], P[pan["S_inf"], 1], color=COL_A,
                           s=S_PT, zorder=3, label=r"$S_{\mathrm{inf}}$")
            if pan["S_sup"]:
                ax.scatter(P[pan["S_sup"], 0], P[pan["S_sup"], 1], color=COL_B,
                           s=S_PT, zorder=3, label=r"$S_{\mathrm{sup}}$")
            ax.scatter(P[[a, b], 0], P[[a, b], 1], color=COL_C, s=S_PT+18, zorder=4)
            etq_generica(ax, a, "$a$", COL_C, dx=6); etq_generica(ax, b, "$b$", COL_C)
            ax.set_title(f"Paso {k+1}: división por $ab$", fontsize=9)
            if pan["S_inf"] or pan["S_sup"]:
                ax.legend(fontsize=7, loc="upper left", framealpha=0.9)

        elif pan["tipo"] == "triangulo":
            hull_previo(ax, pan["previas"])
            nivel_es_raiz = (pan["p"], pan["q"]) in [(global_a, global_b), (global_b, global_a)]
            ext_a, ext_b, c = pan["p"], pan["q"], pan["c"]
            tri = np.array([P[ext_a], P[ext_b], P[c]])
            ax.fill(tri[:, 0], tri[:, 1], color=COL_C, alpha=0.08, zorder=1)
            if tuple(sorted((ext_a, c))) not in pan["previas"]:
                edge(ax, ext_a, c, COL_C, 1.2, 2, ls=":")
            if tuple(sorted((c, ext_b))) not in pan["previas"]:
                edge(ax, c, ext_b, COL_C, 1.2, 2, ls=":")
            ax.annotate("", xy=P[ext_b], xytext=P[ext_a],
                        arrowprops=dict(arrowstyle="-|>", color=COL_BASE, lw=1.6), zorder=3)
            if pan["A"]:
                ax.scatter(P[pan["A"], 0], P[pan["A"], 1], color=COL_A, s=S_PT, zorder=4, label="$A$")
            if pan["B"]:
                ax.scatter(P[pan["B"], 0], P[pan["B"], 1], color=COL_B, s=S_PT, zorder=4, label="$B$")
            if pan["descartados"]:
                ax.scatter(P[pan["descartados"], 0], P[pan["descartados"], 1], color=COL_DESC,
                           marker="x", s=S_PT+6, zorder=4, label="descartados")
            ax.scatter(P[[ext_a, ext_b], 0], P[[ext_a, ext_b], 1], color=COL_HULL, s=S_PT+6, zorder=7)
            ax.scatter([P[c, 0]], [P[c, 1]], color=COL_C, s=S_PT+30, zorder=7)
            etq_generica(ax, ext_a, "$a$", COL_HULL, dx=6, dy=-10)
            etq_generica(ax, ext_b, "$b$", COL_HULL, dx=-20, dy=-10)
            etq_generica(ax, c, "$c$", COL_C, dx=6, dy=6)
            if nivel_es_raiz:
                arg = r"S_{\mathrm{inf}}" if (ext_a, ext_b) == (global_a, global_b) else r"S_{\mathrm{sup}}"
                a1 = "a" if ext_a == global_a else "b"
                a2 = "b" if ext_b == global_b else "a"
                ax.set_title(rf"Paso {k+1}: $\mathrm{{QuickHull}}({a1}, {a2}, {arg})$", fontsize=8.5)
            else:
                ax.set_title(rf"Paso {k+1}: $\mathrm{{QuickHull}}(a, b, T)$, $|T|={len(pan['T'])}$",
                             fontsize=8.5)
            if pan["A"] or pan["B"] or pan["descartados"]:
                ax.legend(fontsize=7, loc="upper left", framealpha=0.9)

        else:  # final
            ciclo = hull + [hull[0]]
            ax.fill(P[ciclo, 0], P[ciclo, 1], color=COL_C, alpha=0.07, zorder=1)
            hull_previo(ax, pan["previas"])
            ax.scatter(P[hull, 0], P[hull, 1], color=COL_C, s=S_PT+18, zorder=7)
            for idx in hull:
                etq_hull(ax, idx, COL_HULL)
            ax.set_title(f"Paso {k+1}: envolvente final", fontsize=9)

        ax.set_xlim(*xlim); ax.set_ylim(*ylim)
        ax.set_aspect("equal"); ax.axis("off")

    for k in range(npan, len(axs)):
        axs[k].axis("off")
    fig.suptitle("Ejecución de QuickHull", fontsize=14)
    plt.tight_layout(rect=[0, 0, 1, 0.97])
    plt.savefig(archivo, bbox_inches="tight")
    plt.show()


if __name__ == "__main__":
    P = np.array([
        [10, 5], [0, 5],          # a (x máx), b (x mín)
        [3, 0], [7, 1],           # vértices inferiores
        [2, 9], [7, 10],          # vértices superiores
        [4, 4], [5, 3],           # interiores
        [6, 6], [5, 7],           # interiores
        [4, 3], [6.5, 2], [3, 8], # interiores varios
        [9, 0.5],                 # exterior inferior derecho
        [8, 9.5],                 # exterior superior derecho
        [2, -1],                  # exterior a la derecha de bc (inferior)
        [4, -2]                   # exterior a la derecha de bc (inferior)
    ], dtype=float)

    hull, eventos = quickhull_pasos(P)
    print("Vértices de la envolvente (antihorario):", hull)
    dibujar_quickhull_pasos(P, eventos, hull, archivo="quickhull_pasos.pdf")
    print("Figura guardada en quickhull_pasos.pdf")
