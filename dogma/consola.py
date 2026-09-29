"""Visualización en modo texto (consola) con colores.

Colores de las bases:
    A verde · T rojo · U magenta · G amarillo · C cian
    Cebadores de ARN (minúsculas en el modelo) -> fondo magenta
    Posiciones sin sintetizar -> '·' en gris

Los colores se generan con códigos ANSI sin depender de librerías externas.
Se pueden desactivar con la variable de entorno NO_COLOR.
"""

import os
import shutil
import sys
import textwrap


def _activar_colores():
    """Activa los códigos de color ANSI y devuelve si se pueden usar.

    En Windows 10/11 hay que activar el modo "virtual terminal" de la consola;
    se hace directamente con la API de Windows, sin librerías externas.
    Se desactivan si la salida no es una terminal o si existe NO_COLOR.
    """
    if os.environ.get("NO_COLOR") or not hasattr(sys.stdout, "isatty") \
            or not sys.stdout.isatty():
        return False
    if os.name == "nt":
        try:
            import ctypes

            kernel32 = ctypes.windll.kernel32
            handle = kernel32.GetStdHandle(-11)  # STD_OUTPUT_HANDLE
            modo = ctypes.c_uint32()
            if kernel32.GetConsoleMode(handle, ctypes.byref(modo)):
                # ENABLE_VIRTUAL_TERMINAL_PROCESSING = 0x0004
                kernel32.SetConsoleMode(handle, modo.value | 0x0004)
        except Exception:
            pass
        os.system("")  # truco adicional que también activa ANSI en cmd
    return True


COLORES = _activar_colores()


class _Ansi:
    """Códigos de escape ANSI (vacíos si no hay colores)."""

    def __init__(self, codigos):
        for nombre, codigo in codigos.items():
            setattr(self, nombre, f"\033[{codigo}m" if COLORES else "")


Fore = _Ansi({"BLACK": 30, "RED": 31, "GREEN": 32, "YELLOW": 33, "BLUE": 34,
              "MAGENTA": 35, "CYAN": 36, "WHITE": 37})
Back = _Ansi({"BLACK": 40, "RED": 41, "GREEN": 42, "YELLOW": 43, "BLUE": 44,
              "MAGENTA": 45, "CYAN": 46, "WHITE": 47})
Style = _Ansi({"BRIGHT": 1, "DIM": 2, "RESET_ALL": 0})

from .replicacion import ENZIMAS as ENZIMAS_REPLICACION
from .traduccion import AMINOACIDOS, CODIGO_GENETICO
from .traduccion import ENZIMAS as ENZIMAS_TRADUCCION
from .transcripcion import ENZIMAS as ENZIMAS_TRANSCRIPCION

ANCHO_BLOQUE = 60   # bases por línea en las secuencias largas
ANCHO_ETIQUETA = 22

COLOR_BASE = {
    "A": Fore.GREEN,
    "T": Fore.RED,
    "U": Fore.MAGENTA,
    "G": Fore.YELLOW,
    "C": Fore.CYAN,
}
COLOR_FASE = {
    "Iniciación": Fore.GREEN,
    "Elongación": Fore.CYAN,
    "Terminación": Fore.MAGENTA,
}


# --------------------------------------------------------------------------
# Utilidades generales
# --------------------------------------------------------------------------
def ancho_terminal():
    return max(60, min(110, shutil.get_terminal_size((100, 30)).columns))


def color(texto, *estilos):
    return "".join(estilos) + texto + Style.RESET_ALL if COLORES else texto


def colorear_secuencia(seq):
    """Colorea una secuencia base a base (ADN, ARN y cebadores)."""
    if not COLORES:
        return seq
    partes = []
    for b in seq:
        if b.islower():  # nucleótido de ARN de un cebador
            partes.append(Back.MAGENTA + Fore.WHITE + Style.BRIGHT + b.upper())
        elif b in COLOR_BASE:
            partes.append(COLOR_BASE[b] + Style.BRIGHT + b)
        elif b == "·":
            partes.append(Style.DIM + b)
        else:
            partes.append(b)
        partes.append(Style.RESET_ALL)
    return "".join(partes)


def titulo(texto):
    w = ancho_terminal()
    print()
    print(color("═" * w, Fore.BLUE, Style.BRIGHT))
    print(color(texto.center(w), Fore.WHITE, Style.BRIGHT))
    print(color("═" * w, Fore.BLUE, Style.BRIGHT))


def subtitulo(texto):
    print()
    print(color(f"── {texto} ", Fore.YELLOW, Style.BRIGHT)
          + color("─" * max(0, ancho_terminal() - len(texto) - 4), Fore.YELLOW))


def parrafo(texto, sangria="   "):
    w = ancho_terminal() - len(sangria)
    for linea in textwrap.wrap(texto, w):
        print(sangria + linea)


def aviso(texto):
    print()
    for i, linea in enumerate(textwrap.wrap(texto, ancho_terminal() - 6)):
        prefijo = " (!) " if i == 0 else "     "
        print(color(prefijo + linea, Fore.YELLOW))


def ok(texto):
    print(color("  [OK] ", Fore.GREEN, Style.BRIGHT) + texto)


def pausa(modo):
    """Espera a que el usuario pulse Enter en el modo paso a paso.

    Devuelve el nuevo modo: si el usuario escribe 'c' se pasa a modo continuo.
    """
    if modo != "paso":
        return modo
    r = input(color("   [Enter] siguiente paso · [c] mostrar el resto seguido > ",
                    Style.DIM)).strip().lower()
    return "continuo" if r == "c" else modo


def cabecera_evento(e, total):
    fase = color(f"{e.fase:<11}", COLOR_FASE.get(e.fase, ""), Style.BRIGHT)
    enzima = color(e.enzima, Fore.WHITE, Style.BRIGHT)
    print()
    print(f" [{e.numero:>2}/{total}] {fase} │ {enzima}")
    print("   " + color(e.accion, Style.BRIGHT))
    parrafo(e.detalle)


# --------------------------------------------------------------------------
# Secuencias alineadas en bloques
# --------------------------------------------------------------------------
def _regla(inicio, fin):
    """Regla de posiciones (1-based) para las posiciones [inicio, fin)."""
    chars = [" "] * (fin - inicio)
    for p in range(inicio + 1, fin + 1):
        if p % 10 == 0:
            s = str(p)
            ini = p - inicio - len(s)
            if ini >= 0:
                chars[ini:ini + len(s)] = s
    return "".join(chars)


def imprimir_alineado(filas, n, ancho=ANCHO_BLOQUE, regla=True):
    """Imprime varias secuencias alineadas, partidas en bloques de ``ancho``.

    ``filas``: lista de tuplas (etiqueta, extremo_izq, texto, extremo_der, tipo)
    donde ``tipo`` es 'seq' (se colorea) o 'marca' (texto plano en un color).
    """
    for ini in range(0, n, ancho):
        fin = min(n, ini + ancho)
        if regla:
            print(" " * (ANCHO_ETIQUETA + 7)
                  + color(_regla(ini, fin), Style.DIM))
        for etiqueta, izq, texto, der, tipo in filas:
            trozo = texto[ini:fin]
            if tipo == "seq":
                cuerpo = colorear_secuencia(trozo)
            else:
                cuerpo = color(trozo, tipo) if tipo else trozo
            izq_txt = izq if ini == 0 else "  "
            der_txt = der if fin == n else ""
            print(f"   {etiqueta:<{ANCHO_ETIQUETA}} {izq_txt:>2} {cuerpo} {der_txt}")
        if fin < n:
            print()


# --------------------------------------------------------------------------
# Leyendas / glosario
# --------------------------------------------------------------------------
def mostrar_leyenda_colores():
    print("   Bases: " + "  ".join(
        colorear_secuencia(b) + f" {n}" for b, n in
        [("A", "Adenina"), ("T", "Timina"), ("G", "Guanina"),
         ("C", "Citosina"), ("U", "Uracilo")]))
    print("   " + colorear_secuencia("aug") + " nucleótidos de ARN de un cebador   "
          + colorear_secuencia("···") + " aún sin sintetizar")


def mostrar_enzimas():
    for nombre_etapa, tabla in (
        ("Replicación (ADN -> ADN)", ENZIMAS_REPLICACION),
        ("Transcripción (ADN -> ARN)", ENZIMAS_TRANSCRIPCION),
        ("Traducción (ARN -> proteína)", ENZIMAS_TRADUCCION),
    ):
        subtitulo(nombre_etapa)
        for nombre, (equivalente, funcion) in tabla.items():
            print("   " + color(nombre, Fore.WHITE, Style.BRIGHT)
                  + color(f"  [{equivalente}]", Style.DIM))
            parrafo(funcion, "      ")


def mostrar_codigo_genetico():
    subtitulo("Código genético estándar (codón de ARNm -> aminoácido)")
    bases = "UCAG"
    for b1 in bases:
        for b3 in bases:
            celdas = []
            for b2 in bases:
                c = b1 + b2 + b3
                aa = CODIGO_GENETICO[c]
                abrev = AMINOACIDOS[aa][0]
                txt = f"{colorear_secuencia(c)} {abrev:<4}"
                if aa == "*":
                    txt = f"{colorear_secuencia(c)} " + color(f"{abrev:<4}", Fore.RED, Style.BRIGHT)
                elif c == "AUG":
                    txt = f"{colorear_secuencia(c)} " + color(f"{abrev:<4}", Fore.GREEN, Style.BRIGHT)
                celdas.append(txt)
            print("   " + "   ".join(celdas))
        print()
    print("   " + color("AUG", Fore.GREEN, Style.BRIGHT) + " = inicio (Metionina)   "
          + color("Stop", Fore.RED, Style.BRIGHT) + " = UAA, UAG, UGA")


# --------------------------------------------------------------------------
# REPLICACIÓN
# --------------------------------------------------------------------------
def _linea_horquilla(n, abierto):
    """Marca de la horquilla: zona abierta, helicasa y zona aún unida."""
    chars = []
    for i in range(n):
        if i < abierto:
            chars.append(" ")
        elif i == abierto:
            chars.append("▲")
        else:
            chars.append("═")
    return "".join(chars)


def _linea_fragmentos(resultado, rezagada):
    chars = [" "] * len(rezagada)
    for f in resultado.fragmentos:
        if any(rezagada[i] != "·" for i in range(f.inicio, f.fin)):
            for i in range(f.inicio, f.fin):
                chars[i] = str(f.numero % 10)
    return "".join(chars)


def dibujar_replicacion(resultado, evento):
    n = len(resultado.hebra_superior)
    filas = [
        ("Hebra parental sup.", "5'", resultado.hebra_superior, "3'", "seq"),
        ("Nueva (rezagada)", "3'", evento.rezagada, "5'", "seq"),
        ("  fragm. de Okazaki", "", _linea_fragmentos(resultado, evento.rezagada), "", Fore.MAGENTA),
        ("  horquilla", "", _linea_horquilla(n, evento.abierto_hasta), "",
         Fore.BLUE + Style.BRIGHT),
        ("Nueva (líder)", "5'", evento.lider, "3'", "seq"),
        ("Hebra parental inf.", "3'", resultado.hebra_inferior, "5'", "seq"),
    ]
    print()
    imprimir_alineado(filas, n)


def mostrar_replicacion(resultado, modo="paso"):
    titulo("1. Replicación del ADN  (ADN -> ADN)")
    parrafo("La horquilla parte del origen de replicación (izquierda) y avanza "
            "hacia la derecha. '▲' marca la helicasa y '═' la doble hélice que "
            "aún no se ha abierto. La cadena líder se copia de forma continua; "
            "la rezagada, en fragmentos de Okazaki numerados.")
    mostrar_leyenda_colores()
    print(f"   Tamaño de fragmento de Okazaki: {resultado.tamano_fragmento} nt · "
          f"cebador: {resultado.longitud_cebador} nt (valores reducidos para verlo)")
    if modo != "resumen":
        total = len(resultado.eventos)
        for e in resultado.eventos:
            cabecera_evento(e, total)
            dibujar_replicacion(resultado, e)
            modo = pausa(modo)
    else:
        dibujar_replicacion(resultado, resultado.eventos[-1])

    subtitulo("Resultado de la replicación")
    h1a, h1b = resultado.molecula_hija_1
    h2a, h2b = resultado.molecula_hija_2
    n = len(h1a)
    print("   Molécula hija 1 (hebra parental superior + nueva rezagada):")
    imprimir_alineado([("parental", "5'", h1a, "3'", "seq"),
                       ("nueva", "3'", h1b, "5'", "seq")], n, regla=False)
    print("   Molécula hija 2 (nueva líder + hebra parental inferior):")
    imprimir_alineado([("nueva", "5'", h2a, "3'", "seq"),
                       ("parental", "3'", h2b, "5'", "seq")], n, regla=False)
    print()
    print(f"   Puentes de hidrógeno rotos por la helicasa: {resultado.puentes_rotos}")
    print(f"   Cebadores de ARN: {1 + len(resultado.fragmentos)} "
          f"(1 en la líder + {len(resultado.fragmentos)} en la rezagada)")
    print(f"   Fragmentos de Okazaki: {len(resultado.fragmentos)} · "
          f"uniones de la ligasa: {resultado.ligaciones}")
    if resultado.es_correcta:
        ok("Las dos moléculas hijas son idénticas a la original "
           "(replicación semiconservativa).")
    for a in resultado.avisos:
        aviso(a)
    return modo


# --------------------------------------------------------------------------
# TRANSCRIPCIÓN
# --------------------------------------------------------------------------
def _linea_burbuja(n, burbuja):
    ini, fin = burbuja
    chars = [" "] * n
    for i in range(ini, fin):
        chars[i] = "~"
    if fin > ini:
        chars[fin - 1] = "▲"
    return "".join(chars)


def dibujar_transcripcion(resultado, evento):
    n = len(resultado.hebra_codificante)
    filas = [
        ("Hebra codificante", "5'", resultado.hebra_codificante, "3'", "seq"),
        ("  ARN polimerasa", "", _linea_burbuja(n, evento.burbuja), "",
         Fore.BLUE + Style.BRIGHT),
        ("Hebra molde", "3'", resultado.hebra_molde, "5'", "seq"),
        ("ARN (nuevo)", "5'", evento.arn, "3'", "seq"),
    ]
    print()
    imprimir_alineado(filas, n)


def mostrar_transcripcion(resultado, modo="paso"):
    titulo("2. Transcripción  (ADN -> ARNm)")
    parrafo("La ARN polimerasa lee la hebra molde (3'->5') y sintetiza el ARNm "
            "(5'->3') con las reglas A->U, T->A, G->C, C->G. '~' marca la "
            "burbuja de transcripción y '▲' la posición de la polimerasa.")
    if modo != "resumen":
        total = len(resultado.eventos)
        for e in resultado.eventos:
            cabecera_evento(e, total)
            dibujar_transcripcion(resultado, e)
            modo = pausa(modo)

    subtitulo("Resultado de la transcripción")
    n = len(resultado.arnm)
    imprimir_alineado([
        ("Hebra molde", "3'", resultado.hebra_molde, "5'", "seq"),
        ("ARNm", "5'", resultado.arnm, "3'", "seq"),
        ("Hebra codificante", "5'", resultado.hebra_codificante, "3'", "seq"),
    ], n)
    comp = resultado.composicion
    print()
    print(f"   Longitud: {n} nt · enlaces fosfodiéster: {resultado.enlaces_fosfodiester} · "
          + " ".join(f"{b}={comp[b]}" for b in "AUGC"))
    if resultado.es_correcta:
        ok("El ARNm tiene la misma secuencia que la hebra codificante con U en "
           "lugar de T.")
    for a in resultado.avisos:
        aviso(a)
    return modo


# --------------------------------------------------------------------------
# TRADUCCIÓN
# --------------------------------------------------------------------------
FONDO_SITIO = {"E": Back.WHITE + Fore.BLACK, "P": Back.BLUE + Fore.WHITE,
               "A": Back.GREEN + Fore.BLACK}


def _sitios_evento(e):
    """Índices de codón que ocupan los sitios E, P y A en un evento."""
    k = e.indice_codon
    if k < 0:
        return {}
    sitios = {}
    if e.sitio_a:
        sitios[k] = "A"
        if e.sitio_p:
            sitios[k - 1] = "P"
        if e.sitio_e:
            sitios[k - 2] = "E"
    else:
        if e.sitio_p:
            sitios[k] = "P"
        if e.sitio_e:
            sitios[k - 1] = "E"
    return sitios


def dibujar_traduccion(resultado, evento, codones_por_linea=12):
    sitios = _sitios_evento(evento)
    codones = resultado.codones
    utr5 = resultado.utr5
    print()
    if utr5:
        print("   5'-UTR: " + color(utr5, Style.DIM))
    for ini in range(0, len(codones), codones_por_linea):
        linea_c, linea_a, linea_s = [], [], []
        for k in range(ini, min(len(codones), ini + codones_por_linea)):
            c = codones[k]
            aa = CODIGO_GENETICO[c]
            abrev = AMINOACIDOS[aa][0]
            if k in sitios:
                linea_c.append(color(f" {c} ", FONDO_SITIO[sitios[k]], Style.BRIGHT))
                linea_s.append(color(f"  {sitios[k]}  ", Style.BRIGHT))
            else:
                linea_c.append(" " + colorear_secuencia(c) + " ")
                linea_s.append("     ")
            linea_a.append(f" {abrev:<4}")
        izq = "5' " if ini == 0 else "   "
        print("   " + izq + "".join(linea_c))
        print("      " + "".join(linea_a))
        if sitios:
            print("      " + "".join(linea_s))
    if resultado.utr3:
        print("   3'-UTR: " + color(resultado.utr3, Style.DIM))
    pep = "-".join(AMINOACIDOS[a][0] for a in evento.peptido)
    print("   Ribosoma: "
          + "  ".join(color(f" {s}: {getattr(evento, 'sitio_' + s.lower()) or '---'} ",
                            FONDO_SITIO[s]) for s in "EPA"))
    print("   Cadena polipeptídica: " + color(f"N- {pep or '(vacía)'} -C", Fore.GREEN, Style.BRIGHT))


def mostrar_traduccion(resultado, modo="paso"):
    titulo("3. Traducción  (ARNm -> proteína)")
    parrafo("El ribosoma lee el ARNm en codones desde el primer AUG. Colores de "
            "los sitios del ribosoma: E (salida), P (peptidil) y A (aminoacil).")
    if resultado.inicio < 0:
        for a in resultado.avisos:
            aviso(a)
        return modo
    if modo != "resumen":
        total = len(resultado.eventos)
        for e in resultado.eventos:
            cabecera_evento(e, total)
            dibujar_traduccion(resultado, e)
            modo = pausa(modo)

    subtitulo("Resultado de la traducción")
    print(f"   Codón de inicio AUG en la posición {resultado.inicio + 1} del ARNm")
    print("   Codones leídos: " + " ".join(colorear_secuencia(c) for c in resultado.codones))
    print("   Proteína (3 letras): " + color(resultado.proteina_3_letras, Fore.GREEN, Style.BRIGHT))
    print("   Proteína (1 letra):  " + color(resultado.proteina, Fore.GREEN, Style.BRIGHT))
    print(f"   Aminoácidos: {len(resultado.proteina)} · enlaces peptídicos: "
          f"{resultado.enlaces_peptidicos} · GTP consumidos: {resultado.gtp_consumidos}")
    for a in resultado.avisos:
        aviso(a)
    return modo


# --------------------------------------------------------------------------
# RESUMEN DEL FLUJO COMPLETO
# --------------------------------------------------------------------------
def mostrar_resumen(rep, tra, trad):
    titulo("Flujo completo de la información genética")
    n = len(rep.hebra_superior)
    linea_codones = [" "] * n
    linea_aa = [" "] * n
    if trad.inicio >= 0:
        for k, c in enumerate(trad.codones):
            p = trad.inicio + 3 * k
            abrev = AMINOACIDOS[CODIGO_GENETICO[c]][0]
            linea_codones[p:p + 3] = "[" + "=" + "]"
            linea_aa[p:p + 3] = abrev[:3] if abrev != "Stop" else "***"
    imprimir_alineado([
        ("ADN hebra codificante", "5'", rep.hebra_superior, "3'", "seq"),
        ("ADN hebra molde", "3'", rep.hebra_inferior, "5'", "seq"),
        ("ARNm", "5'", tra.arnm, "3'", "seq"),
        ("  codones", "", "".join(linea_codones), "", Style.DIM),
        ("Proteína", "N-", "".join(linea_aa), "-C", Fore.GREEN + Style.BRIGHT),
    ], n)
    print()
    print("   " + color("ADN -> ADN", Fore.BLUE, Style.BRIGHT)
          + f"   Replicación: 2 moléculas hijas de {n} pb, "
          f"{len(rep.fragmentos)} fragmentos de Okazaki")
    print("   " + color("ADN -> ARN", Fore.BLUE, Style.BRIGHT)
          + f"   Transcripción: ARNm de {len(tra.arnm)} nt")
    print("   " + color("ARN -> Prot", Fore.BLUE, Style.BRIGHT)
          + f"  Traducción: {trad.proteina_3_letras or 'sin proteína'}")
