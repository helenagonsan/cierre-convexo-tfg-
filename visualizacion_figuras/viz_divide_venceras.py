# -*- coding: utf-8 -*-
"""
Visualización de la fase de fusión del algoritmo Divide y Vencerás.

Trabajo de Fin de Grado: "El cierre convexo: Teoría y aplicaciones"
Autora: Helena González Sánchez — Universidad de Salamanca

Script autónomo. Ilustra el cálculo de las tangentes comunes inferior y superior
que fusionan dos envolventes convexas disjuntas A (izquierda) y B (derecha) en
CH(A ∪ B). Cada panel muestra un estado del avance de los extremos de la tangente
sobre las cadenas de ambos polígonos, hasta alcanzar la doble tangencia.

Ejecución:  python3 viz_divide_venceras.py
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


def graham_scan(P):
    """Graham Scan auxiliar para construir las dos envolventes A y B de partida."""
    P = np.asarray(P, float); n = len(P)
    if n < 3:
        return list(range(n))
    p0 = min(range(n), key=lambda i: (P[i][1], -P[i][0]))
    def cmp(i, j):
        o = orientacion(P[p0], P[i], P[j])
        if o > 0: return -1
        if o < 0: return 1
        di = (P[i][0]-P[p0][0])**2 + (P[i][1]-P[p0][1])**2
        dj = (P[j][0]-P[p0][0])**2 + (P[j][1]-P[p0][1])**2
        return -1 if di < dj else 1
    resto = sorted([i for i in range(n) if i != p0], key=functools.cmp_to_key(cmp))
    pila = [p0, resto[0]]
    for i in resto[1:]:
        while len(pila) >= 2 and orientacion(P[pila[-2]], P[pila[-1]], P[i]) <= 0:
            pila.pop()
        pila.append(i)
    return pila


# ---------------------------------------------------------------------
# Tangentes instrumentadas: registran cada estado (i, j) del avance
#   inferior: a avanza horario en A (i-1), b antihorario en B (j+1)
#   superior: a avanza antihorario en A (i+1), b horario en B (j-1)
# ---------------------------------------------------------------------

def tangente_inferior_pasos(P, A, B):
    nA, nB = len(A), len(B)
    i = max(range(nA), key=lambda k: (P[A[k]][0], P[A[k]][1]))
    j = min(range(nB), key=lambda k: (P[B[k]][0], P[B[k]][1]))
    rec = [(i, j)]
    while True:
        mov = False
        while orientacion(P[B[j]], P[A[i]], P[A[(i-1) % nA]]) > 0:
            i = (i-1) % nA; rec.append((i, j)); mov = True
        while orientacion(P[A[i]], P[B[j]], P[B[(j+1) % nB]]) < 0:
            j = (j+1) % nB; rec.append((i, j)); mov = True
        if not mov:
            return i, j, rec


def tangente_superior_pasos(P, A, B):
    nA, nB = len(A), len(B)
    i = max(range(nA), key=lambda k: (P[A[k]][0], P[A[k]][1]))
    j = min(range(nB), key=lambda k: (P[B[k]][0], P[B[k]][1]))
    rec = [(i, j)]
    while True:
        mov = False
        while orientacion(P[B[j]], P[A[i]], P[A[(i+1) % nA]]) < 0:
            i = (i+1) % nA; rec.append((i, j)); mov = True
        while orientacion(P[A[i]], P[B[j]], P[B[(j-1) % nB]]) > 0:
            j = (j-1) % nB; rec.append((i, j)); mov = True
        if not mov:
            return i, j, rec


def fusionar(P, A, B):
    nA, nB = len(A), len(B)
    inf_a, inf_b, _ = tangente_inferior_pasos(P, A, B)
    sup_a, sup_b, _ = tangente_superior_pasos(P, A, B)
    # CCW: desde a_sup por A (cadena exterior izquierda) hasta a_inf,
    # luego desde b_inf por B (cadena exterior derecha) hasta b_sup.
    res = [A[sup_a]]; k = sup_a
    while k != inf_a:
        k = (k+1) % nA; res.append(A[k])
    k = inf_b; res.append(B[k])
    while k != sup_b:
        k = (k+1) % nB; res.append(B[k])
    return res


# ---------------------------------------------------------------------
# Dibujo de la fusión paso a paso
# ---------------------------------------------------------------------

def dibujar_fusion(P, A, B, ncols=3, archivo="dc_fusion.pdf"):
    plt.close('all')
    P = np.asarray(P, float)
    COL_A, COL_B = "#2563eb", "#d97706"
    COL_SEG, COL_TAN, COL_HULL = "#dc2626", "#16a34a", "#1f2937"
    S_PT = 28

    inf_a, inf_b, rec_inf = tangente_inferior_pasos(P, A, B)
    sup_a, sup_b, rec_sup = tangente_superior_pasos(P, A, B)
    hull = fusionar(P, A, B)

    paneles = []
    for t, (i, j) in enumerate(rec_inf):
        paneles.append(("inf", i, j, t == len(rec_inf)-1))
    for t, (i, j) in enumerate(rec_sup):
        paneles.append(("sup", i, j, t == len(rec_sup)-1))
    paneles.append(("final", None, None, True))

    xmin, xmax = P[:, 0].min(), P[:, 0].max()
    ymin, ymax = P[:, 1].min(), P[:, 1].max()
    mx, my = 0.12*(xmax-xmin), 0.12*(ymax-ymin)
    xlim, ylim = (xmin-mx, xmax+mx), (ymin-my, ymax+my)
    xsep = (P[A[max(range(len(A)), key=lambda k: P[A[k]][0])]][0]
            + P[B[min(range(len(B)), key=lambda k: P[B[k]][0])]][0]) / 2

    npan = len(paneles)
    nrows = (npan + ncols - 1) // ncols
    fig, axs = plt.subplots(nrows, ncols, figsize=(3.8*ncols, 3.8*nrows))
    axs = np.atleast_1d(axs).ravel()

    def poligono(ax, H, color, lw=2.0, z=3):
        m = len(H)
        for t in range(m):
            a, b = H[t], H[(t+1) % m]
            ax.plot([P[a, 0], P[b, 0]], [P[a, 1], P[b, 1]], color=color, lw=lw, zorder=z)

    for k, (clase, i, j, ultimo) in enumerate(paneles):
        ax = axs[k]
        ax.axvline(xsep, color="0.85", ls="--", lw=1, zorder=0)
        ax.scatter(P[:, 0], P[:, 1], color=COL_HULL, s=S_PT-8, zorder=4)

        if clase != "final":
            poligono(ax, A, COL_A)
            poligono(ax, B, COL_B)
            a_idx, b_idx = A[i], B[j]
            color_seg = COL_TAN if ultimo else COL_SEG
            ls = "-" if ultimo else "--"
            ax.plot([P[a_idx, 0], P[b_idx, 0]], [P[a_idx, 1], P[b_idx, 1]],
                    color=color_seg, lw=2.6, ls=ls, zorder=6)
            ax.scatter([P[a_idx, 0], P[b_idx, 0]], [P[a_idx, 1], P[b_idx, 1]],
                       color=color_seg, s=S_PT+30, zorder=7)
            ax.annotate("$a$", P[a_idx], textcoords="offset points", xytext=(-16, -4),
                        fontsize=11, color=color_seg, zorder=9, annotation_clip=False)
            ax.annotate("$b$", P[b_idx], textcoords="offset points", xytext=(8, -4),
                        fontsize=11, color=color_seg, zorder=9, annotation_clip=False)
            nombre = "inferior" if clase == "inf" else "superior"
            primero = (clase == "inf" and k == 0) or (clase == "sup" and paneles[k-1][0] == "inf")
            if ultimo:
                titulo = f"Paso {k+1}: tangente {nombre} hallada"
            elif primero:
                titulo = f"Paso {k+1}: tangente {nombre}, candidato inicial"
            else:
                titulo = f"Paso {k+1}: tangente {nombre}, avance"
            ax.set_title(titulo, fontsize=9)
        else:
            ax.fill(P[hull + [hull[0]], 0], P[hull + [hull[0]], 1],
                    color=COL_TAN, alpha=0.06, zorder=1)
            poligono(ax, hull, COL_HULL, lw=2.4, z=5)
            ax.plot([P[A[inf_a], 0], P[B[inf_b], 0]], [P[A[inf_a], 1], P[B[inf_b], 1]],
                    color=COL_TAN, lw=2.4, zorder=6)
            ax.plot([P[A[sup_a], 0], P[B[sup_b], 0]], [P[A[sup_a], 1], P[B[sup_b], 1]],
                    color=COL_TAN, lw=2.4, zorder=6)
            ax.scatter(P[hull, 0], P[hull, 1], color=COL_SEG, s=S_PT+18, zorder=7)
            ax.set_title(f"Paso {k+1}: fusión $CH(A \\cup B)$", fontsize=9)

        ax.set_xlim(*xlim); ax.set_ylim(*ylim)
        ax.set_aspect("equal"); ax.axis("off")

    for k in range(npan, len(axs)):
        axs[k].axis("off")

    fig.suptitle("Fusión de Divide y Vencerás: tangentes comunes", fontsize=14)
    plt.tight_layout(rect=[0, 0, 1, 0.97])
    plt.savefig(archivo, bbox_inches="tight")
    plt.show()


if __name__ == "__main__":
    # Polígonos diseñados A (izquierda) y B (derecha); cada uno con un punto interior.
    ptsA = np.array([[0,4],[1,1],[3,0],[4,3],[3,6],[1,6],[2,3]], dtype=float)   # (2,3) interior
    ptsB = np.array([[7,3],[8,0],[10,1],[11,4],[9,7],[7,6],[9,3]], dtype=float) # (9,3) interior
    P = np.vstack([ptsA, ptsB])

    A = graham_scan(ptsA)                              # CCW, índices en P (parte A)
    B = [len(ptsA) + i for i in graham_scan(ptsB)]     # CCW, índices en P (parte B)

    print("A (antihorario):", A)
    print("B (antihorario):", B)
    print("Fusión CH(A ∪ B):", fusionar(P, A, B))
    dibujar_fusion(P, A, B, ncols=3, archivo="dc_fusion.pdf")
    print("Figura guardada en dc_fusion.pdf")
