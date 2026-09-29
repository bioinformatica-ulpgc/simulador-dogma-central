"""Visualización de cada etapa en imágenes (PNG) con matplotlib.

Genera cuatro figuras en la carpeta ``resultados/``:

    1_replicacion.png    Horquilla de replicación (instantánea) + moléculas hijas
    2_transcripcion.png  Burbuja de transcripción + ARNm resultante
    3_traduccion.png     Ribosoma con sus sitios E/P/A, ARNt y proteína final
    4_flujo_completo.png ADN -> ARNm -> proteína alineados

Los colores de las bases son los mismos que en la consola.
"""

import os

import matplotlib

matplotlib.use("Agg")  # genera los PNG sin abrir ventanas

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, Ellipse, FancyBboxPatch, Patch, RegularPolygon

from .traduccion import AMINOACIDOS, CODIGO_GENETICO, anticodon

# --------------------------------------------------------------------------
# Estilo
# --------------------------------------------------------------------------
COLOR_BASE = {
    "A": "#2E9E44",  # verde
    "T": "#D43C3C",  # rojo
    "U": "#B03FB0",  # magenta
    "G": "#E0A000",  # ámbar
    "C": "#1A8FB5",  # cian azulado
}
COLOR_PARENTAL = "#5B6270"
COLOR_NUEVA = "#2F6FDB"
COLOR_CEBADOR = "#8E2C8E"
COLOR_TEXTO = "#222222"
COLOR_SUAVE = "#8A8F98"

ENZIMA_COLOR = {
    "helicasa": "#F28E2B",
    "topo": "#9C6ADE",
    "ssb": "#B8BCC4",
    "pol3": "#3AA0A0",
    "primasa": "#C0508C",
    "pol1": "#6C9A2B",
    "ligasa": "#C9A227",
    "arnpol": "#F28E2B",
    "ribosoma": "#7FA7D9",
}

# Clasificación de aminoácidos por las propiedades de su cadena lateral
CLASE_AA = {
    **dict.fromkeys("AVLIMFWP", ("Apolar (hidrófobo)", "#E8A33D")),
    **dict.fromkeys("GSTCYNQ", ("Polar sin carga", "#5DAE8B")),
    **dict.fromkeys("KRH", ("Básico (+)", "#4C78C9")),
    **dict.fromkeys("DE", ("Ácido (-)", "#D1495B")),
}

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "axes.titlesize": 13,
    "axes.titleweight": "bold",
    "axes.titlecolor": COLOR_TEXTO,
})


def _tamano_letra(n):
    return max(5, min(10, 420 / max(n, 1)))


def _ancho_figura(n):
    return min(40, max(13, 0.3 * n + 4))


def _preparar_eje(ax, n, ymin, ymax, izq=-3.2):
    ax.set_xlim(izq, n + 2.2)
    ax.set_ylim(ymin, ymax)
    ax.set_aspect("auto")
    ax.axis("off")


def _base(ax, x, y, letra, fs, arn=False, color=None, alpha=1.0, texto=True):
    """Dibuja un nucleótido como un cuadrado de color con su letra."""
    letra_mayus = letra.upper()
    if arn and letra_mayus == "T":
        letra_mayus = "U"
    fc = color or COLOR_BASE.get(letra_mayus, "#999999")
    caja = FancyBboxPatch(
        (x - 0.42, y - 0.34), 0.84, 0.68,
        boxstyle="round,pad=0.02,rounding_size=0.12",
        fc=fc, ec=COLOR_CEBADOR if arn else "white",
        lw=1.6 if arn else 0.6, hatch="////" if arn else None, alpha=alpha,
        zorder=3,
    )
    ax.add_patch(caja)
    if texto:
        ax.text(x, y, letra_mayus, ha="center", va="center", fontsize=fs,
                color="white", fontweight="bold", zorder=4, alpha=alpha)


def _extremos(ax, x_izq, x_der, y, izq, der, fs=10):
    ax.text(x_izq, y, izq, ha="right", va="center", fontsize=fs,
            color=COLOR_TEXTO, fontweight="bold")
    ax.text(x_der, y, der, ha="left", va="center", fontsize=fs,
            color=COLOR_TEXTO, fontweight="bold")


def _minuscula(texto):
    return texto[:1].lower() + texto[1:]


def _etiqueta_fila(ax, y, texto, fs=9):
    ax.text(-3.0, y, texto, ha="left", va="center", fontsize=fs, color=COLOR_SUAVE)


def _enzima(ax, x, y, nombre, clave, ancho=2.4, alto=0.9, fs=8, forma="elipse",
            alpha=0.92):
    color = ENZIMA_COLOR[clave]
    if forma == "hexagono":
        p = RegularPolygon((x, y), numVertices=6, radius=alto * 0.75, fc=color,
                           ec="white", lw=1.5, alpha=alpha, zorder=6)
    else:
        p = Ellipse((x, y), ancho, alto, fc=color, ec="white", lw=1.5,
                    alpha=alpha, zorder=6)
    ax.add_patch(p)
    ax.text(x, y, nombre, ha="center", va="center", fontsize=fs, color="white",
            fontweight="bold", zorder=7)


def _guardar(fig, carpeta, nombre):
    os.makedirs(carpeta, exist_ok=True)
    ruta = os.path.join(carpeta, nombre)
    fig.savefig(ruta, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return ruta


def _leyenda_bases(ax, incluir_u=True, extra=()):
    bases = "ATGCU" if incluir_u else "ATGC"
    nombres = {"A": "Adenina", "T": "Timina", "G": "Guanina", "C": "Citosina",
               "U": "Uracilo"}
    handles = [Patch(fc=COLOR_BASE[b], label=f"{b} {nombres[b]}") for b in bases]
    handles += list(extra)
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.0),
              ncol=min(len(handles), 8), frameon=False, fontsize=9)


# --------------------------------------------------------------------------
# 1. REPLICACIÓN
# --------------------------------------------------------------------------
def _elegir_instantanea_replicacion(rep):
    """Momento con al menos dos fragmentos de Okazaki aún sin unir."""
    for e in rep.eventos:
        if e.enzima == "Primasa" and "fragmento de Okazaki 2" in e.accion:
            return e
    # Molécula corta con un solo fragmento: último evento de elongación
    elong = [e for e in rep.eventos if e.fase == "Elongación"]
    return elong[-1]


def figura_replicacion(rep, carpeta="resultados"):
    n = len(rep.hebra_superior)
    fs = _tamano_letra(n)
    e = _elegir_instantanea_replicacion(rep)
    abierto = e.abierto_hasta
    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(_ancho_figura(n), 11.5),
        gridspec_kw={"height_ratios": [1.45, 1]})

    # ---------------- Panel A: horquilla ----------------
    Y_SUP, Y_REZ, Y_LID, Y_INF = 2.6, 1.55, -1.55, -2.6
    Y_SUP_C, Y_INF_C = 0.45, -0.45

    def y_sup(i):
        return Y_SUP if i < abierto else Y_SUP_C

    def y_inf(i):
        return Y_INF if i < abierto else Y_INF_C

    # esqueletos azúcar-fosfato
    xs = list(range(n))
    ax1.plot(xs, [y_sup(i) for i in xs], color=COLOR_PARENTAL, lw=2.2, zorder=1)
    ax1.plot(xs, [y_inf(i) for i in xs], color=COLOR_PARENTAL, lw=2.2, zorder=1)
    # puentes de hidrógeno en la doble hélice aún cerrada
    for i in range(abierto, n):
        ax1.plot([i, i], [Y_INF_C + 0.34, Y_SUP_C - 0.34], color="#C9CDD3",
                 lw=1, zorder=2)

    for i in range(n):
        _base(ax1, i, y_sup(i), rep.hebra_superior[i], fs)
        _base(ax1, i, y_inf(i), rep.hebra_inferior[i], fs)
        if e.rezagada[i] != "·":
            _base(ax1, i, Y_REZ, e.rezagada[i], fs, arn=e.rezagada[i].islower())
            ax1.plot([i, i], [Y_REZ + 0.34, Y_SUP - 0.34], color="#C9CDD3", lw=1)
        if e.lider[i] != "·":
            _base(ax1, i, Y_LID, e.lider[i], fs, arn=e.lider[i].islower())
            ax1.plot([i, i], [Y_LID - 0.34, Y_INF + 0.34], color="#C9CDD3", lw=1)

    _extremos(ax1, -0.7, n - 0.3, Y_SUP_C if abierto < n else Y_SUP, "", "3'")
    _extremos(ax1, -0.7, n - 0.3, Y_INF_C if abierto < n else Y_INF, "", "5'")
    _extremos(ax1, -0.7, n, Y_SUP, "5'", "")
    _extremos(ax1, -0.7, n, Y_INF, "3'", "")
    _extremos(ax1, -0.7, n, Y_REZ, "3'", "")
    _extremos(ax1, -0.7, n, Y_LID, "5'", "")
    _etiqueta_fila(ax1, Y_SUP + 0.75, "Hebra parental (molde de la rezagada)")
    _etiqueta_fila(ax1, Y_INF - 0.75, "Hebra parental (molde de la líder)")

    # Cadena líder: flecha continua
    fin_lider = max(i for i in range(n) if e.lider[i] != "·")
    ax1.annotate("", xy=(fin_lider + 0.4, -0.8), xytext=(-0.4, -0.8),
                 arrowprops=dict(arrowstyle="-|>", color=COLOR_NUEVA, lw=2.2))
    ax1.text(fin_lider / 2, -0.45, "Cadena líder: síntesis CONTINUA 5'→3'",
             ha="center", va="center", fontsize=9.5, color=COLOR_NUEVA,
             fontweight="bold")

    # Cadena rezagada: una flecha por fragmento de Okazaki (sentido 5'->3' = hacia la izquierda)
    fragmentos_visibles = [f for f in rep.fragmentos
                           if any(e.rezagada[i] != "·" for i in range(f.inicio, f.fin))]
    for f in fragmentos_visibles:
        ax1.annotate("", xy=(f.inicio - 0.4, 0.85), xytext=(f.fin - 0.6, 0.85),
                     arrowprops=dict(arrowstyle="-|>", color=COLOR_CEBADOR, lw=2))
        ax1.text((f.inicio + f.fin - 1) / 2, 0.5, f"Fragmento de Okazaki {f.numero}",
                 ha="center", va="center", fontsize=8.5, color=COLOR_CEBADOR,
                 fontweight="bold")

    # Enzimas
    if abierto < n:
        _enzima(ax1, abierto - 0.6, 0, "Helicasa", "helicasa", alto=1.4,
                fs=7.5, forma="hexagono", alpha=0.88)
        x_topo = min(n - 2, abierto + max(4, n // 10))
        _enzima(ax1, x_topo, 1.5, "Topoisomerasa", "topo", ancho=3.6, alto=0.9, fs=7.5)
        ax1.annotate("", xy=(x_topo, 0.85), xytext=(x_topo, 1.05),
                     arrowprops=dict(arrowstyle="-|>", color=ENZIMA_COLOR["topo"]))
    for i in range(abierto):
        if e.rezagada[i] == "·":
            ax1.add_patch(Circle((i, Y_SUP - 0.62), 0.22, fc=ENZIMA_COLOR["ssb"],
                                 ec="white", zorder=5))
        if e.lider[i] == "·":
            ax1.add_patch(Circle((i, Y_INF + 0.62), 0.22, fc=ENZIMA_COLOR["ssb"],
                                 ec="white", zorder=5))
    _enzima(ax1, fin_lider + 1.6, Y_LID, "ADN pol III", "pol3", ancho=3.0, fs=7.5)
    if fragmentos_visibles:
        ultimo = fragmentos_visibles[-1]
        x_ceb = (ultimo.inicio_cebador + ultimo.fin - 1) / 2
        _enzima(ax1, ultimo.inicio_cebador - 1.8, Y_REZ, "ADN pol III ⟵", "pol3",
                ancho=3.4, fs=7.5, alpha=0.8)
        ax1.annotate("Primasa: acaba de poner el cebador de ARN\n"
                     f"del fragmento {ultimo.numero}", xy=(x_ceb, Y_REZ + 0.35),
                     xytext=(x_ceb - 6, Y_SUP + 1.3), fontsize=8.5,
                     color=ENZIMA_COLOR["primasa"], fontweight="bold",
                     arrowprops=dict(arrowstyle="->", color=ENZIMA_COLOR["primasa"]))
        if len(fragmentos_visibles) >= 2:
            previo = fragmentos_visibles[-2]
            x_prev = (previo.inicio_cebador + previo.fin - 1) / 2
            ax1.annotate(f"Cebador del fragmento {previo.numero}: después la ADN pol I\n"
                         "lo cambia por ADN y la ADN ligasa une los fragmentos",
                         xy=(x_prev, Y_REZ + 0.35), xytext=(0, Y_SUP + 1.3),
                         fontsize=8.5, color=ENZIMA_COLOR["pol1"], fontweight="bold",
                         arrowprops=dict(arrowstyle="->", color=ENZIMA_COLOR["pol1"]))
    ax1.text(n + 1.5, Y_SUP + 1.3, "Avance de la horquilla  ➜", ha="right",
             va="center", fontsize=10, color=ENZIMA_COLOR["helicasa"], fontweight="bold")
    ax1.text(-3.0, Y_SUP + 1.3, "ori", ha="left", va="center", fontsize=10,
             color=COLOR_TEXTO, fontweight="bold")

    ax1.set_title(f"Replicación del ADN — horquilla de replicación (paso {e.numero} "
                  f"de {len(rep.eventos)}: {_minuscula(e.accion)})", loc="left")
    _preparar_eje(ax1, n, Y_INF - 1.3, Y_SUP + 1.9)
    extra = [
        Patch(fc="white", ec=COLOR_CEBADOR, hatch="////", label="Cebador de ARN"),
        Line2D([], [], marker="o", ls="", mfc=ENZIMA_COLOR["ssb"], mec="white",
               ms=9, label="Proteínas SSB"),
    ]
    _leyenda_bases(ax1, extra=extra)

    # ---------------- Panel B: moléculas hijas ----------------
    filas = [
        (2.6, rep.hebra_superior, COLOR_PARENTAL, "5'", "3'", "Parental"),
        (1.8, rep.cadena_rezagada, COLOR_NUEVA, "3'", "5'", "Nueva (rezagada)"),
        (-0.4, rep.cadena_lider, COLOR_NUEVA, "5'", "3'", "Nueva (líder)"),
        (-1.2, rep.hebra_inferior, COLOR_PARENTAL, "3'", "5'", "Parental"),
    ]
    for y, seq, col, izq, der, etiqueta in filas:
        ax2.plot([0, n - 1], [y, y], color=col, lw=2, zorder=1)
        for i, b in enumerate(seq):
            _base(ax2, i, y, b, fs, color=col)
        _extremos(ax2, -0.7, n - 0.3, y, izq, der)
        ax2.text(-1.8, y, etiqueta, ha="right", va="center", fontsize=8.5, color=col,
                 fontweight="bold")
    for i in range(n):
        ax2.plot([i, i], [1.8 + 0.34, 2.6 - 0.34], color="#C9CDD3", lw=1)
        ax2.plot([i, i], [-1.2 + 0.34, -0.4 - 0.34], color="#C9CDD3", lw=1)
    ax2.text(n / 2, 3.35, "Molécula hija 1", ha="center", fontsize=10.5,
             fontweight="bold", color=COLOR_TEXTO)
    ax2.text(n / 2, 0.35, "Molécula hija 2", ha="center", fontsize=10.5,
             fontweight="bold", color=COLOR_TEXTO)
    ax2.set_title(
        f"Resultado: replicación semiconservativa — {len(rep.fragmentos)} fragmentos "
        f"de Okazaki, {1 + len(rep.fragmentos)} cebadores, {rep.ligaciones} uniones "
        f"de la ligasa, {rep.puentes_rotos} puentes de H rotos", loc="left")
    _preparar_eje(ax2, n, -2.0, 3.8, izq=-9)
    ax2.legend(handles=[Patch(fc=COLOR_PARENTAL, label="Hebra parental (se conserva)"),
                        Patch(fc=COLOR_NUEVA, label="Hebra nueva (sintetizada)")],
               loc="upper center", bbox_to_anchor=(0.5, 0.02), ncol=2,
               frameon=False, fontsize=9)
    fig.tight_layout(h_pad=3)
    return _guardar(fig, carpeta, "1_replicacion.png")


# --------------------------------------------------------------------------
# 2. TRANSCRIPCIÓN
# --------------------------------------------------------------------------
def figura_transcripcion(tra, carpeta="resultados"):
    n = len(tra.hebra_codificante)
    fs = _tamano_letra(n)
    elong = [e for e in tra.eventos if e.fase == "Elongación"]
    e = elong[len(elong) // 2] if elong else tra.eventos[-1]
    ini, fin = e.burbuja
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(_ancho_figura(n), 9.5),
                                   gridspec_kw={"height_ratios": [1.3, 1]})

    # ---------------- Panel A: burbuja ----------------
    Y_COD_C, Y_COD_A, Y_MOL, Y_HIB, Y_SAL = 0.5, 1.9, -0.5, -1.5, -2.7
    HIBRIDO = 8  # nt del híbrido ARN-ADN dentro de la burbuja

    def y_cod(i):
        return Y_COD_A if ini <= i < fin else Y_COD_C

    ax1.plot(range(n), [y_cod(i) for i in range(n)], color=COLOR_PARENTAL, lw=2.2)
    ax1.plot([0, n - 1], [Y_MOL, Y_MOL], color=COLOR_PARENTAL, lw=2.2)
    sintetizadas = [i for i in range(n) if e.arn[i] != "·"]

    ultimo_arn = (max(sintetizadas) + 1) if sintetizadas else 0

    def y_arn(i):
        return Y_HIB if i >= ultimo_arn - HIBRIDO else Y_SAL

    if sintetizadas:
        ax1.plot(sintetizadas, [y_arn(i) for i in sintetizadas],
                 color=COLOR_BASE["U"], lw=2)
    for i in range(n):
        _base(ax1, i, y_cod(i), tra.hebra_codificante[i], fs)
        _base(ax1, i, Y_MOL, tra.hebra_molde[i], fs)
        if not (ini <= i < fin):
            ax1.plot([i, i], [Y_MOL + 0.34, Y_COD_C - 0.34], color="#C9CDD3", lw=1)
        if e.arn[i] != "·":
            _base(ax1, i, y_arn(i), e.arn[i], fs)
            if y_arn(i) == Y_HIB:
                ax1.plot([i, i], [Y_HIB + 0.34, Y_MOL - 0.34], color="#C9CDD3", lw=1)
    if fin > ini:
        ax1.add_patch(Ellipse(((ini + fin - 1) / 2, -0.1), fin - ini + 3, 4.6,
                              fc=ENZIMA_COLOR["arnpol"], ec=ENZIMA_COLOR["arnpol"],
                              alpha=0.16, lw=2, zorder=0))
        ax1.text((ini + fin - 1) / 2, 3.0, "ARN polimerasa  ➜", ha="center",
                 fontsize=10, fontweight="bold", color=ENZIMA_COLOR["arnpol"])
    _extremos(ax1, -0.7, n - 0.3, Y_COD_C, "5'", "3'")
    _extremos(ax1, -0.7, n - 0.3, Y_MOL, "3'", "5'")
    if sintetizadas:
        _extremos(ax1, sintetizadas[0] - 0.7, sintetizadas[-1] + 0.7,
                  y_arn(sintetizadas[0]), "5'", "")
        ax1.text(sintetizadas[-1] + 0.7, Y_HIB, "3'", ha="left", va="center",
                 fontsize=10, fontweight="bold")
    _etiqueta_fila(ax1, Y_COD_C + 0.9, "Hebra codificante")
    _etiqueta_fila(ax1, Y_MOL - 0.75, "Hebra molde")
    _etiqueta_fila(ax1, Y_SAL - 0.7, "ARNm saliendo de la polimerasa")
    ax1.text(n + 1.8, -3.5,
             "Complementariedad ADN molde → ARN:   A→U   T→A   G→C   C→G",
             ha="right", fontsize=9.5, color=COLOR_TEXTO,
             bbox=dict(boxstyle="round", fc="#F4F5F7", ec="#D5D8DD"))
    ax1.set_title(f"Transcripción — burbuja de transcripción (paso {e.numero} de "
                  f"{len(tra.eventos)}: {_minuscula(e.accion)})", loc="left")
    _preparar_eje(ax1, n, -4.1, 3.5)
    _leyenda_bases(ax1)

    # ---------------- Panel B: resultado ----------------
    filas = [(1.6, tra.hebra_codificante, "5'", "3'", "Hebra codificante"),
             (0.8, tra.hebra_molde, "3'", "5'", "Hebra molde"),
             (-0.6, tra.arnm, "5'", "3'", "ARNm")]
    for y, seq, izq, der, etiqueta in filas:
        for i, b in enumerate(seq):
            _base(ax2, i, y, b, fs)
        _extremos(ax2, -0.7, n - 0.3, y, izq, der)
        ax2.text(-1.8, y, etiqueta, ha="right", va="center", fontsize=8.5,
                 color=COLOR_TEXTO, fontweight="bold")
    for i in range(n):
        ax2.plot([i, i], [0.8 + 0.34, 1.6 - 0.34], color="#C9CDD3", lw=1)
        if tra.hebra_codificante[i] == "T":
            ax2.add_patch(FancyBboxPatch((i - 0.5, -1.05), 1.0, 2.9,
                                         boxstyle="round,pad=0,rounding_size=0.2",
                                         fc="none", ec=COLOR_BASE["U"], lw=1,
                                         ls="--", zorder=5))
    comp = tra.composicion
    ax2.set_title(f"Resultado: ARNm de {len(tra.arnm)} nt (A={comp['A']}, U={comp['U']}, "
                  f"G={comp['G']}, C={comp['C']}) — igual a la hebra codificante "
                  "cambiando T por U (recuadros)", loc="left")
    _preparar_eje(ax2, n, -1.5, 2.3, izq=-9)
    fig.tight_layout(h_pad=3)
    return _guardar(fig, carpeta, "2_transcripcion.png")


# --------------------------------------------------------------------------
# 3. TRADUCCIÓN
# --------------------------------------------------------------------------
def _aa_circulo(ax, x, y, letra, r=0.62, fs=7.5, alpha=1.0):
    nombre_clase, col = CLASE_AA.get(letra, ("", "#999999"))
    ax.add_patch(Circle((x, y), r, fc=col, ec="white", lw=1.5, zorder=6, alpha=alpha))
    ax.text(x, y, AMINOACIDOS[letra][0], ha="center", va="center", fontsize=fs,
            color="white", fontweight="bold", zorder=7)


def figura_traduccion(trad, carpeta="resultados"):
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(15, 10.5),
                                   gridspec_kw={"height_ratios": [1.5, 1]})
    if trad.inicio < 0:
        for ax in (ax1, ax2):
            ax.axis("off")
        ax1.text(0.5, 0.5, "El ARNm no contiene ningún codón AUG:\nno se sintetiza proteína.",
                 ha="center", va="center", fontsize=16, transform=ax1.transAxes)
        return _guardar(fig, carpeta, "3_traduccion.png")

    # instantánea: un ARNt entrando en el sitio A con al menos 2 aa en la cadena
    candidatos = [e for e in trad.eventos
                  if e.sitio_a and e.fase == "Elongación" and "entra ARNt" in e.accion]
    e = next((c for c in candidatos if len(c.peptido) >= 3), None) \
        or (candidatos[len(candidatos) // 2] if candidatos else None)
    if e is None:  # proteína de un solo aminoácido: se muestra la iniciación
        e = next(c for c in trad.eventos if c.sitio_p and c.indice_codon == 0)

    # ---------------- Panel A: ribosoma ----------------
    k = e.indice_codon
    # ventana de codones a mostrar alrededor del ribosoma
    desde = max(0, k - 4)
    hasta = min(len(trad.codones), k + 5)
    x0 = 0
    fs = 10
    posiciones = {}
    for j in range(desde, hasta):
        c = trad.codones[j]
        for m, b in enumerate(c):
            _base(ax1, x0 + m, 0, b, fs)
        posiciones[j] = x0 + 1
        ax1.text(x0 + 1, -0.8, f"codón {j + 1}", ha="center", fontsize=7.5,
                 color=COLOR_SUAVE)
        x0 += 3.4
    ancho = x0
    ax1.plot([-1, ancho], [0, 0], color=COLOR_BASE["U"], lw=2, zorder=1)
    ax1.text(-1.3, 0, "5'" if desde == 0 else "5' …", ha="right", va="center",
             fontsize=10, fontweight="bold")
    ax1.text(ancho, 0, "3'" if hasta == len(trad.codones) else "… 3'", ha="left",
             va="center", fontsize=10, fontweight="bold")

    sitios = {}
    if e.sitio_a:
        sitios = {k: "A", k - 1: "P"}
        if e.sitio_e:
            sitios[k - 2] = "E"
    else:
        sitios = {k: "P"}
    x_min = posiciones[min(sitios)] - 2.2
    x_max = posiciones[max(sitios)] + 2.2
    xc = (x_min + x_max) / 2
    ax1.add_patch(Ellipse((xc, -0.9), x_max - x_min + 1.5, 2.0,
                          fc=ENZIMA_COLOR["ribosoma"], alpha=0.35, zorder=0))
    ax1.add_patch(Ellipse((xc, 2.6), x_max - x_min + 2.5, 4.6,
                          fc=ENZIMA_COLOR["ribosoma"], alpha=0.22, zorder=0))
    ax1.text(xc, -1.75, "Subunidad pequeña", ha="center", fontsize=8.5,
             color="#35587F", fontweight="bold")
    ax1.text(x_max + 1.6, 3.6, "Subunidad grande\n(peptidil transferasa)",
             ha="left", va="center", fontsize=8.5, color="#35587F", fontweight="bold")

    peptido = e.peptido
    for j, sitio in sitios.items():
        x = posiciones[j]
        codon = trad.codones[j]
        vacio = sitio == "E"
        alpha = 0.45 if vacio else 1.0
        # ARNt: cuerpo + anticodón
        ax1.add_patch(FancyBboxPatch((x - 1.1, 1.1), 2.2, 1.9,
                                     boxstyle="round,pad=0.05,rounding_size=0.4",
                                     fc="#EDE3F3", ec=COLOR_BASE["U"], lw=1.5,
                                     alpha=alpha, zorder=2))
        for m, b in enumerate(anticodon(codon)):
            _base(ax1, x - 1 + m, 0.78, b, 9, alpha=alpha)
        ax1.text(x, 2.1, f"ARNt\n3'-{anticodon(codon)}-5'", ha="center",
                 va="center", fontsize=7.5, color="#5A2A6E", alpha=alpha, zorder=3)
        ax1.text(x, -2.35, f"Sitio {sitio}", ha="center", fontsize=10,
                 fontweight="bold", color=COLOR_TEXTO)
        if sitio == "A":
            _aa_circulo(ax1, x, 3.6, CODIGO_GENETICO[codon])
        if sitio == "P":
            # cadena peptídica unida al ARNt del sitio P, sale hacia arriba-izquierda
            for m, letra in enumerate(reversed(peptido)):
                _aa_circulo(ax1, x - 1.25 * m, 3.6 + 0.85 * m, letra)
            ax1.text(x - 1.25 * (len(peptido) - 1) - 0.9,
                     3.6 + 0.85 * (len(peptido) - 1), "N-", ha="right",
                     va="center", fontsize=10, fontweight="bold")
    ax1.annotate("", xy=(ancho + 0.2, -3.0), xytext=(-1, -3.0),
                 arrowprops=dict(arrowstyle="-|>", color=COLOR_SUAVE, lw=1.5))
    ax1.text(ancho / 2, -3.35, "El ribosoma avanza 5'→3' de codón en codón",
             ha="center", fontsize=9, color=COLOR_SUAVE)
    ax1.set_xlim(-4, ancho + 2)
    ax1.set_ylim(-3.8, 3.6 + 0.85 * max(3, len(peptido)) + 0.6)
    ax1.set_aspect("equal")
    ax1.axis("off")
    ax1.set_title(f"Traducción — ribosoma durante la elongación (paso {e.numero} de "
                  f"{len(trad.eventos)}: {e.accion})", loc="left")

    # ---------------- Panel B: proteína final ----------------
    prot = trad.proteina
    for m, letra in enumerate(prot):
        if m:
            ax2.plot([m - 1 + 0.45, m - 0.45], [0, 0], color=COLOR_TEXTO, lw=1.5)
        _aa_circulo(ax2, m, 0, letra, r=0.42, fs=7 if len(prot) < 30 else 5.5)
    ax2.text(-0.8, 0, "N-", ha="right", va="center", fontsize=11, fontweight="bold")
    ax2.text(len(prot) - 0.2, 0, "-C", ha="left", va="center", fontsize=11,
             fontweight="bold")
    ax2.text(-0.8, -1.2, f"Secuencia (1 letra): {prot}", fontsize=10,
             color=COLOR_TEXTO, family="monospace")
    ax2.set_xlim(-2.5, max(len(prot), 12) + 1.5)
    ax2.set_ylim(-1.8, 1.2)
    ax2.set_aspect("equal")
    ax2.axis("off")
    ax2.set_title(f"Resultado: proteína de {len(prot)} aminoácidos, "
                  f"{trad.enlaces_peptidicos} enlaces peptídicos, "
                  f"{trad.gtp_consumidos} GTP consumidos", loc="left")
    clases = []
    for nombre, col in dict(CLASE_AA.values()).items():
        clases.append(Patch(fc=col, label=nombre))
    ax2.legend(handles=clases, loc="upper center", bbox_to_anchor=(0.5, 0.05),
               ncol=4, frameon=False, fontsize=9, title="Tipo de aminoácido")
    fig.tight_layout(h_pad=2)
    return _guardar(fig, carpeta, "3_traduccion.png")


# --------------------------------------------------------------------------
# 4. FLUJO COMPLETO
# --------------------------------------------------------------------------
def figura_flujo(rep, tra, trad, carpeta="resultados"):
    n = len(rep.hebra_superior)
    fs = _tamano_letra(n)
    fig, ax = plt.subplots(figsize=(_ancho_figura(n), 7))
    Y = {"cod": 6.4, "mol": 5.6, "arn": 2.8, "prot": -0.4}
    for i in range(n):
        _base(ax, i, Y["cod"], rep.hebra_superior[i], fs)
        _base(ax, i, Y["mol"], rep.hebra_inferior[i], fs)
        _base(ax, i, Y["arn"], tra.arnm[i], fs)
        ax.plot([i, i], [Y["mol"] + 0.34, Y["cod"] - 0.34], color="#C9CDD3", lw=1)
    _extremos(ax, -0.7, n - 0.3, Y["cod"], "5'", "3'")
    _extremos(ax, -0.7, n - 0.3, Y["mol"], "3'", "5'")
    _extremos(ax, -0.7, n - 0.3, Y["arn"], "5'", "3'")
    ax.text(-3.0, (Y["cod"] + Y["mol"]) / 2, "ADN", fontsize=11, fontweight="bold",
            va="center")
    ax.text(-3.0, Y["arn"], "ARNm", fontsize=11, fontweight="bold", va="center")
    ax.text(-3.0, Y["prot"], "Proteína", fontsize=11, fontweight="bold", va="center")

    if trad.inicio >= 0:
        for k, c in enumerate(trad.codones):
            p = trad.inicio + 3 * k
            ax.add_patch(FancyBboxPatch((p - 0.48, Y["arn"] - 0.45), 2.96, 0.9,
                                        boxstyle="round,pad=0,rounding_size=0.15",
                                        fc="none", ec=COLOR_TEXTO, lw=1, zorder=5))
            aa = CODIGO_GENETICO[c]
            if aa == "*":
                ax.text(p + 1, Y["prot"], "STOP", ha="center", va="center",
                        fontsize=8, fontweight="bold", color=COLOR_BASE["T"])
            else:
                _aa_circulo(ax, p + 1, Y["prot"], aa, r=0.75,
                            fs=max(5.5, fs - 1))
            ax.plot([p + 1, p + 1], [Y["arn"] - 0.5, Y["prot"] + 0.8],
                    color="#D5D8DD", lw=1, zorder=0)
    # flechas del dogma
    xf = n + 1.2
    for (y1, y2, txt) in [(Y["cod"] + 0.9, Y["cod"] + 0.9, None),
                          (Y["mol"] - 0.5, Y["arn"] + 0.5, "Transcripción"),
                          (Y["arn"] - 0.5, Y["prot"] + 0.8, "Traducción")]:
        if txt:
            ax.annotate("", xy=(xf, y2), xytext=(xf, y1),
                        arrowprops=dict(arrowstyle="-|>", color=COLOR_NUEVA, lw=2.5))
            ax.text(xf + 0.4, (y1 + y2) / 2, txt, va="center", fontsize=10,
                    color=COLOR_NUEVA, fontweight="bold")
    ax.text(n / 2, Y["cod"] + 1.0,
            f"Replicación: ADN → 2 moléculas de ADN idénticas ({len(rep.fragmentos)} "
            "fragmentos de Okazaki en la cadena rezagada)",
            ha="center", fontsize=10, color=COLOR_NUEVA, fontweight="bold")
    ax.set_xlim(-3.2, n + 6)
    ax.set_ylim(Y["prot"] - 1.6, Y["cod"] + 1.5)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title("Flujo de la información genética: ADN → ADN, ADN → ARN, ARN → proteína",
                 loc="left")
    prot_txt = trad.proteina_3_letras or "sin proteína"
    ax.text(-3.0, Y["prot"] - 1.3, f"Proteína: N- {prot_txt} -C", fontsize=10,
            color=COLOR_TEXTO)
    fig.tight_layout()
    return _guardar(fig, carpeta, "4_flujo_completo.png")


def generar_imagenes(rep, tra, trad, carpeta="resultados"):
    """Genera las cuatro figuras y devuelve la lista de rutas."""
    return [
        figura_replicacion(rep, carpeta),
        figura_transcripcion(tra, carpeta),
        figura_traduccion(trad, carpeta),
        figura_flujo(rep, tra, trad, carpeta),
    ]
