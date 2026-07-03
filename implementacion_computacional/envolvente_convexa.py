# -*- coding: utf-8 -*-
"""
Implementación en Python de los algoritmos para el cálculo de la envolvente
convexa en el plano.

Trabajo de Fin de Grado: "El cierre convexo: Teoría y aplicaciones"
Autora: Helena González Sánchez
Tutor:  Ángel María Martín del Rey
Grado en Matemáticas — Universidad de Salamanca
Curso académico 2025-2026

Este fichero corresponde al Anexo de implementación de la memoria. Contiene los
cinco algoritmos estudiados (Gift Wrapping, QuickHull, Graham Scan, Incremental
y Divide y Vencerás), construidos sobre un único predicado geométrico elemental:
el predicado de orientación, evaluado como el signo de un determinante y exacto
en aritmética entera.

Todas las funciones reciben una lista de puntos P, donde cada punto es una tupla
(x, y), y devuelven los índices de los vértices de CH(P) en orden antihorario.

Requisitos: Python 3.x (sólo biblioteca estándar). El bloque de prueba final no
requiere dependencias externas.
"""

import functools


# =====================================================================
# Primitivas geométricas comunes a los cinco algoritmos
# =====================================================================

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


# =====================================================================
# Gift Wrapping (marcha de Jarvis)
# =====================================================================

def gift_wrapping(P):
    """
    Gift Wrapping (marcha de Jarvis).
    Devuelve los índices de los vértices de CH(P) en orden antihorario.
    """
    n = len(P)
    if n < 3:
        return list(range(n))

    inicio = min(range(n), key=lambda i: (P[i][1], -P[i][0]))

    hull = []
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
        actual = siguiente
        if actual == inicio:
            break
    return hull


# =====================================================================
# QuickHull
# =====================================================================

def _qh(P, a, b, T):
    """Cadena de vértices de CH entre a y b; T = puntos a la izquierda de ->ab."""
    if not T:
        return []
    c = max(T, key=lambda i: orientacion(P[a], P[b], P[i]))
    A = [i for i in T if orientacion(P[a], P[c], P[i]) > 0]
    B = [i for i in T if orientacion(P[c], P[b], P[i]) > 0]
    return _qh(P, a, c, A) + [c] + _qh(P, c, b, B)


def quickhull(P):
    """
    QuickHull en R^2 (posición general).
    Devuelve los índices de los vértices de CH(P) en orden antihorario.
    """
    n = len(P)
    if n < 3:
        return list(range(n))

    a = max(range(n), key=lambda i: (P[i][0], P[i][1]))    # x máx (desempate y máx)
    b = min(range(n), key=lambda i: (P[i][0], -P[i][1]))   # x mín (desempate y máx)

    S_inf = [i for i in range(n) if orientacion(P[a], P[b], P[i]) > 0]
    S_sup = [i for i in range(n) if orientacion(P[a], P[b], P[i]) < 0]

    hull = [a] + _qh(P, a, b, S_inf) + [b] + _qh(P, b, a, S_sup)
    hull = hull[::-1]   # la cadena se construye horaria; se invierte a antihoraria
    return hull


# =====================================================================
# Graham Scan
# =====================================================================

def graham_scan(P):
    """
    Graham Scan en R^2 (posición general).
    Devuelve los índices de los vértices de CH(P) en orden antihorario.
    """
    n = len(P)
    if n < 3:
        return list(range(n))

    # Ancla p0: y mínima, desempate x máxima.
    p0 = min(range(n), key=lambda i: (P[i][1], -P[i][0]))

    # Orden angular antihorario respecto de p0; empate -> más cercano primero.
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

    # Barrido con pila. Giro estrictamente a la izquierda (convención estricta).
    pila = [p0, resto[0]]
    for i in resto[1:]:
        while len(pila) >= 2 and orientacion(P[pila[-2]], P[pila[-1]], P[i]) <= 0:
            pila.pop()
        pila.append(i)
    return pila


# =====================================================================
# Algoritmo Incremental
# =====================================================================

def incremental(P):
    """
    Algoritmo incremental en R^2 (posición general).
    Procesa los puntos en el orden dado (sin reordenar). Si los primeros puntos
    son colineales, localiza el primer punto que rompe la colinealidad y forma
    el triángulo inicial con los dos extremos del tramo alineado y dicho punto.
    Devuelve los índices de los vértices de CH(P) en orden antihorario.
    """
    n = len(P)
    if n < 3:
        return list(range(n))

    # Inicialización: primer punto no colineal con P[0], P[1].
    t = 2
    while t < n and orientacion(P[0], P[1], P[t]) == 0:
        t += 1
    if t == n:
        # Todos los puntos son colineales: la envolvente es el segmento extremo.
        a = min(range(n), key=lambda i: (P[i][0], P[i][1]))
        b = max(range(n), key=lambda i: (P[i][0], P[i][1]))
        return [a, b]

    # Extremos del tramo colineal P[0..t-1]; triángulo inicial antihorario.
    a = min(range(t), key=lambda i: (P[i][0], P[i][1]))
    b = max(range(t), key=lambda i: (P[i][0], P[i][1]))
    if orientacion(P[a], P[b], P[t]) > 0:
        H = [a, b, t]
    else:
        H = [a, t, b]

    procesados = {a, b, t}
    for k in range(n):
        if k in procesados:
            continue
        m = len(H)
        # Pertenencia: p_k es interior si D >= 0 en todas las aristas.
        interior = all(orientacion(P[H[i]], P[H[(i+1) % m]], P[k]) >= 0
                       for i in range(m))
        if interior:
            continue

        # Punto exterior: conservar los vértices cuya arista entrante no es
        # visible, e insertar p_k entre los dos puntos de tangencia.
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
    return H


# =====================================================================
# Divide y Vencerás
# =====================================================================

def divide_venceras(P):
    """
    Divide y Vencerás en R^2 (posición general, x distintas).
    Devuelve los índices de los vértices de CH(P) en orden antihorario.
    """
    n = len(P)
    if n < 2:
        return list(range(n))
    orden = sorted(range(n), key=lambda i: (P[i][0], P[i][1]))   # preorden por x
    return _dc(P, orden)


def _dc(P, idxs):
    if len(idxs) <= 3:                       # caso base: segmento o triángulo antihorario
        if len(idxs) <= 2:
            return list(idxs)
        a, b, c = idxs
        return [a, b, c] if orientacion(P[a], P[b], P[c]) > 0 else [a, c, b]
    m = (len(idxs) + 1) // 2                  # ceil(|S| / 2)
    A = _dc(P, idxs[:m])                       # mitad izquierda (antihoraria)
    B = _dc(P, idxs[m:])                       # mitad derecha (antihoraria)
    return _fusionar(P, A, B)


def _avanza(P, pivote, actual, siguiente, signo):
    """
    Decide si avanzar de 'actual' a 'siguiente' sobre una cadena durante la
    búsqueda de tangentes. Avanza si 'siguiente' queda del lado indicado por
    'signo' respecto de la recta dirigida pivote -> actual (signo = +1 a la
    izquierda, -1 a la derecha). En empate colineal (orientacion = 0), avanza
    hacia el vértice más lejano del pivote, descartando los puntos intermedios
    de un tramo alineado.
    """
    o = orientacion(P[pivote], P[actual], P[siguiente])
    if signo * o > 0:
        return True
    if o == 0:
        return distancia_sq(P[pivote], P[siguiente]) > distancia_sq(P[pivote], P[actual])
    return False


def _tangente_inferior(P, A, B):
    """Tangente común inferior: a avanza horario en A, b antihorario en B."""
    nA, nB = len(A), len(B)
    i = max(range(nA), key=lambda k: (P[A[k]][0], P[A[k]][1]))   # x máx de A
    j = min(range(nB), key=lambda k: (P[B[k]][0], P[B[k]][1]))   # x mín de B
    for _ in range(nA + nB + 1):                                 # tope de seguridad
        movido = False
        while _avanza(P, B[j], A[i], A[(i - 1) % nA], +1):
            i = (i - 1) % nA
            movido = True
        while _avanza(P, A[i], B[j], B[(j + 1) % nB], -1):
            j = (j + 1) % nB
            movido = True
        if not movido:
            break
    return i, j


def _tangente_superior(P, A, B):
    """Tangente común superior: a avanza antihorario en A, b horario en B."""
    nA, nB = len(A), len(B)
    i = max(range(nA), key=lambda k: (P[A[k]][0], P[A[k]][1]))
    j = min(range(nB), key=lambda k: (P[B[k]][0], P[B[k]][1]))
    for _ in range(nA + nB + 1):                                 # tope de seguridad
        movido = False
        while _avanza(P, B[j], A[i], A[(i + 1) % nA], -1):
            i = (i + 1) % nA
            movido = True
        while _avanza(P, A[i], B[j], B[(j - 1) % nB], +1):
            j = (j - 1) % nB
            movido = True
        if not movido:
            break
    return i, j


def _fusionar(P, A, B):
    nA, nB = len(A), len(B)
    inf_a, inf_b = _tangente_inferior(P, A, B)
    sup_a, sup_b = _tangente_superior(P, A, B)
    # Concatenación antihoraria (Proposición de fusión): desde a_sup por A hasta
    # a_inf, y desde b_inf por B hasta b_sup.
    res = [A[sup_a]]
    k = sup_a
    while k != inf_a:
        k = (k + 1) % nA
        res.append(A[k])
    k = inf_b
    res.append(B[k])
    while k != sup_b:
        k = (k + 1) % nB
        res.append(B[k])
    return res


# =====================================================================
# Bloque de prueba
# =====================================================================

if __name__ == "__main__":
    import random

    ALGORITMOS = {
        "Gift Wrapping":    gift_wrapping,
        "QuickHull":        quickhull,
        "Graham Scan":      graham_scan,
        "Incremental":      incremental,
        "Divide y Vencerás": divide_venceras,
    }

    def conjunto_vertices(hull):
        return frozenset(hull)

    # ---- Ejemplo ilustrativo ----
    P = [(0, 0), (1, 4), (5, 5), (6, 1), (3, 2), (2, 3), (4, 1), (3, 4)]
    print("Conjunto de puntos de ejemplo:")
    for i, p in enumerate(P):
        print(f"  P[{i}] = {p}")
    print()
    print("Vértices de la envolvente convexa (índices, orden antihorario):")
    for nombre, fn in ALGORITMOS.items():
        print(f"  {nombre:18s}: {fn(P)}")
    print()

    # ---- Verificación de concordancia sobre casos aleatorios ----
    # En posición general (sin tres puntos colineales), los cinco algoritmos
    # deben coincidir en el conjunto de vértices de CH(P).
    def en_posicion_general(P):
        n = len(P)
        for i in range(n):
            for j in range(i + 1, n):
                for k in range(j + 1, n):
                    if orientacion(P[i], P[j], P[k]) == 0:
                        return False
        return True

    random.seed(0)
    casos, discrepancias = 0, 0
    while casos < 2000:
        n = random.randint(3, 30)
        Q, vistos = [], set()
        while len(Q) < n:
            pt = (random.randint(-100, 100), random.randint(-100, 100))
            if pt not in vistos:
                vistos.add(pt)
                Q.append(pt)
        if not en_posicion_general(Q):
            continue
        casos += 1
        referencia = conjunto_vertices(gift_wrapping(Q))
        for fn in ALGORITMOS.values():
            if conjunto_vertices(fn(Q)) != referencia:
                discrepancias += 1
                break

    print(f"Verificación sobre {casos} conjuntos aleatorios en posición general:")
    if discrepancias == 0:
        print("  Los cinco algoritmos coinciden en todos los casos.")
    else:
        print(f"  ATENCIÓN: {discrepancias} discrepancias detectadas.")
