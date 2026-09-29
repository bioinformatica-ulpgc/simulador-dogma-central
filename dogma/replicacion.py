"""Replicación del ADN (ADN -> ADN).

Simula una horquilla de replicación que parte de un origen situado en el
extremo izquierdo de la molécula y avanza hacia la derecha.

Modelo
------
Hebra superior (parental, codificante):  5' ----------------> 3'
Hebra inferior (parental, molde):         3' <---------------- 5'
Horquilla:                                       avanza -->

* Las ADN polimerasas solo sintetizan en dirección 5' -> 3' y necesitan un
  extremo 3'-OH libre (un cebador) para empezar.
* **Cadena líder**: se copia sobre la hebra inferior (3'->5'). La nueva hebra
  crece 5'->3' en el mismo sentido en que avanza la horquilla, por lo que se
  sintetiza de forma CONTINUA a partir de un único cebador.
* **Cadena rezagada**: se copia sobre la hebra superior (5'->3'). La nueva
  hebra debe crecer en sentido contrario al avance de la horquilla, así que se
  sintetiza de forma DISCONTINUA en **fragmentos de Okazaki**, cada uno con su
  propio cebador. Después la ADN polimerasa I sustituye los cebadores de ARN
  por ADN y la ADN ligasa une los fragmentos.

Representación del estado
-------------------------
Las cadenas nuevas se guardan como textos de la misma longitud que el ADN:
    '·'          -> posición aún no sintetizada
    minúsculas   -> nucleótidos de ARN del cebador (a, u, g, c)
    MAYÚSCULAS   -> nucleótidos de ADN (A, T, G, C)
Ambas cadenas se escriben alineadas con la hebra superior (de izquierda a
derecha), así que la cadena rezagada nueva se lee 3' -> 5'.

Los nombres de las enzimas son los de procariotas (E. coli), indicando entre
paréntesis el equivalente en eucariotas.
"""

from dataclasses import dataclass, field

from .secuencias import (
    complementaria_adn,
    puentes_hidrogeno,
    validar_adn,
)

VACIO = "·"

# --------------------------------------------------------------------------
# Moléculas y enzimas que intervienen (se usan también como leyenda)
# --------------------------------------------------------------------------
ENZIMAS = {
    "Proteínas iniciadoras": (
        "DnaA (ORC en eucariotas)",
        "Reconocen el origen de replicación (ori) y separan las primeras bases.",
    ),
    "Helicasa": (
        "DnaB (MCM en eucariotas)",
        "Rompe los puentes de hidrógeno entre las bases y abre la doble hélice.",
    ),
    "Topoisomerasa": (
        "ADN girasa (topoisomerasa I/II en eucariotas)",
        "Corta y vuelve a unir el ADN por delante de la horquilla para aliviar "
        "el superenrollamiento que genera la apertura.",
    ),
    "Proteínas SSB": (
        "SSB (RPA en eucariotas)",
        "Se unen al ADN de cadena sencilla y evitan que las hebras se vuelvan "
        "a emparejar.",
    ),
    "Primasa": (
        "DnaG (ADN pol α-primasa en eucariotas)",
        "ARN polimerasa que sintetiza un cebador corto de ARN con un extremo "
        "3'-OH libre.",
    ),
    "ADN polimerasa III": (
        "Pol III holoenzima (ADN pol ε en la líder y δ en la rezagada)",
        "Añade desoxirribonucleótidos (dNTP) en dirección 5'->3' siguiendo la "
        "complementariedad A-T y G-C. Tiene actividad correctora 3'->5'.",
    ),
    "Pinza deslizante": (
        "Pinza β (PCNA en eucariotas)",
        "Anillo que mantiene a la polimerasa unida al ADN (procesividad).",
    ),
    "ADN polimerasa I": (
        "Pol I (RNasa H + FEN1 + ADN pol δ en eucariotas)",
        "Elimina los cebadores de ARN (exonucleasa 5'->3') y rellena el hueco "
        "con ADN.",
    ),
    "ADN ligasa": (
        "Ligasa (ADN ligasa I en eucariotas)",
        "Forma el enlace fosfodiéster que une fragmentos contiguos (sella las "
        "mellas).",
    ),
}


@dataclass
class Evento:
    """Un paso de la simulación, con una foto del estado en ese momento."""

    numero: int
    fase: str          # "Iniciación", "Elongación" o "Terminación"
    enzima: str        # enzima o molécula protagonista
    accion: str        # resumen en una línea
    detalle: str       # explicación
    abierto_hasta: int # nº de pares de bases ya abiertos por la helicasa
    lider: str         # estado de la cadena líder nueva
    rezagada: str      # estado de la cadena rezagada nueva
    mellas: tuple      # posiciones i donde hay una mella entre i-1 e i (rezagada)


@dataclass
class FragmentoOkazaki:
    numero: int
    inicio: int            # posición izquierda (extremo 3' del fragmento)
    fin: int               # posición derecha exclusiva (extremo 5')
    inicio_cebador: int    # el cebador ocupa [inicio_cebador, fin)
    cebador_arn: str       # secuencia del cebador escrita 5'->3'
    ligado: bool = False


@dataclass
class ResultadoReplicacion:
    hebra_superior: str            # parental 5'->3'
    hebra_inferior: str            # parental 3'->5'
    cadena_lider: str              # nueva, 5'->3'
    cadena_rezagada: str           # nueva, alineada, se lee 3'->5'
    cebador_lider: str             # cebador de la cadena líder (5'->3')
    fragmentos: list
    eventos: list
    longitud_cebador: int
    tamano_fragmento: int
    puentes_rotos: int
    ligaciones: int
    avisos: list = field(default_factory=list)

    @property
    def molecula_hija_1(self):
        """(hebra parental superior, nueva cadena rezagada)."""
        return (self.hebra_superior, self.cadena_rezagada)

    @property
    def molecula_hija_2(self):
        """(nueva cadena líder, hebra parental inferior)."""
        return (self.cadena_lider, self.hebra_inferior)

    @property
    def es_correcta(self):
        """Las dos moléculas hijas son idénticas a la original."""
        return (
            self.cadena_lider == self.hebra_superior
            and self.cadena_rezagada == self.hebra_inferior
        )


def _a_arn(base_adn):
    """Nucleótido de ARN en minúscula equivalente a una base de ADN."""
    return "u" if base_adn == "T" else base_adn.lower()


def tamano_fragmento_automatico(n, longitud_cebador):
    """Tamaño de fragmento de Okazaki para obtener unos 4 fragmentos."""
    return max(longitud_cebador + 3, -(-n // 4))


def replicar(adn, tamano_fragmento=None, longitud_cebador=3):
    """Simula la replicación semiconservativa de ``adn`` (hebra 5'->3').

    * ``tamano_fragmento``: bases que abre la helicasa en cada avance y, por
      tanto, tamaño de los fragmentos de Okazaki. En la célula son ~1000-2000
      nt en procariotas y ~100-200 nt en eucariotas; aquí se reduce para que
      se pueda ver. Si es ``None`` se calcula automáticamente.
    * ``longitud_cebador``: longitud de los cebadores de ARN (~10 nt en la
      célula).
    """
    validar_adn(adn)
    n = len(adn)
    if longitud_cebador < 1:
        raise ValueError("El cebador debe tener al menos 1 nucleótido.")
    if n < 2 * longitud_cebador:
        raise ValueError(
            f"La secuencia es demasiado corta: necesita al menos "
            f"{2 * longitud_cebador} pb para colocar los cebadores."
        )
    if tamano_fragmento is None:
        tamano_fragmento = tamano_fragmento_automatico(n, longitud_cebador)
    if tamano_fragmento <= longitud_cebador:
        raise ValueError("El fragmento de Okazaki debe ser más largo que el cebador.")

    superior = adn
    inferior = complementaria_adn(adn)

    # Estado de las cadenas nuevas (listas mutables de caracteres)
    lider = [VACIO] * n
    rezagada = [VACIO] * n
    mellas = set()
    eventos = []
    fragmentos = []
    avisos = []
    abierto = 0
    puentes_rotos = 0
    ligaciones = 0

    def registrar(fase, enzima, accion, detalle):
        eventos.append(
            Evento(
                numero=len(eventos) + 1,
                fase=fase,
                enzima=enzima,
                accion=accion,
                detalle=detalle,
                abierto_hasta=abierto,
                lider="".join(lider),
                rezagada="".join(rezagada),
                mellas=tuple(sorted(mellas)),
            )
        )

    # ------------------------------------------------------------------
    # INICIACIÓN
    # ------------------------------------------------------------------
    registrar(
        "Iniciación", "Proteínas iniciadoras",
        "Reconocimiento del origen de replicación (ori)",
        "Las proteínas iniciadoras se unen al origen, situado en el extremo "
        "izquierdo, y desestabilizan la doble hélice para que entre la helicasa.",
    )

    # ------------------------------------------------------------------
    # ELONGACIÓN: la horquilla avanza en tramos del tamaño de un fragmento
    # ------------------------------------------------------------------
    L = longitud_cebador
    fin_anterior = 0
    while fin_anterior < n:
        fin = min(n, fin_anterior + tamano_fragmento)

        # 1) Helicasa abre el tramo [fin_anterior, fin)
        tramo = superior[fin_anterior:fin]
        rotos = puentes_hidrogeno(tramo)
        puentes_rotos += rotos
        abierto = fin
        registrar(
            "Elongación", "Helicasa",
            f"Abre la doble hélice hasta la posición {fin}",
            f"Rompe {rotos} puentes de hidrógeno entre las bases "
            f"{fin_anterior + 1}-{fin} (2 por cada par A-T y 3 por cada G-C).",
        )

        # 2) Topoisomerasa y SSB
        if fin < n:
            registrar(
                "Elongación", "Topoisomerasa",
                "Alivia la tensión por delante de la horquilla",
                "Al desenrollar la hélice se acumulan superenrollamientos "
                "positivos delante; la topoisomerasa los elimina cortando y "
                "volviendo a unir el ADN.",
            )
        registrar(
            "Elongación", "Proteínas SSB",
            "Estabilizan las hebras sencillas expuestas",
            f"Se unen a las posiciones {fin_anterior + 1}-{fin} de ambas hebras "
            "molde para que no vuelvan a emparejarse.",
        )

        # 3) Cadena líder
        if fin_anterior == 0:
            for i in range(L):
                lider[i] = _a_arn(superior[i])
            registrar(
                "Elongación", "Primasa",
                "Sintetiza el ÚNICO cebador de la cadena líder",
                f"Cebador de ARN 5'-{''.join(lider[:L]).upper()}-3' "
                f"complementario a la hebra inferior (A-U, T-A, G-C, C-G). "
                "Deja un extremo 3'-OH libre.",
            )
            inicio_sintesis = L
        else:
            inicio_sintesis = fin_anterior
        for i in range(inicio_sintesis, fin):
            lider[i] = superior[i]
        registrar(
            "Elongación", "ADN polimerasa III",
            f"Cadena líder: síntesis continua hasta la posición {fin}",
            f"Añade {fin - inicio_sintesis} dNTP en dirección 5'->3' "
            "(mismo sentido que la horquilla), leyendo la hebra inferior "
            "3'->5'. La pinza deslizante la mantiene unida: no necesita "
            "nuevos cebadores.",
        )

        # 4) Cadena rezagada: nuevo fragmento de Okazaki [fin_anterior, fin)
        num = len(fragmentos) + 1
        lf = min(L, fin - fin_anterior)
        ini_ceb = fin - lf
        for i in range(ini_ceb, fin):
            rezagada[i] = _a_arn(inferior[i])
        cebador = "".join(rezagada[ini_ceb:fin])[::-1].upper()  # 5'->3'
        frag = FragmentoOkazaki(num, fin_anterior, fin, ini_ceb, cebador)
        fragmentos.append(frag)
        registrar(
            "Elongación", "Primasa",
            f"Cadena rezagada: cebador del fragmento de Okazaki {num}",
            f"Coloca el cebador 5'-{cebador}-3' junto a la horquilla "
            f"(posiciones {ini_ceb + 1}-{fin}) sobre la hebra superior.",
        )
        for i in range(fin_anterior, ini_ceb):
            rezagada[i] = inferior[i]
        registrar(
            "Elongación", "ADN polimerasa III",
            f"Cadena rezagada: sintetiza el fragmento de Okazaki {num}",
            f"Añade {ini_ceb - fin_anterior} dNTP en dirección 5'->3', es "
            "decir, ALEJÁNDOSE de la horquilla, hasta llegar al fragmento "
            "anterior. Por eso la síntesis es discontinua.",
        )

        # 5) Maduración del fragmento anterior: Pol I + ligasa
        if num > 1:
            previo = fragmentos[-2]
            _sustituir_cebador(rezagada, inferior, previo.inicio_cebador, previo.fin)
            mellas.add(previo.fin)
            registrar(
                "Elongación", "ADN polimerasa I",
                f"Sustituye el cebador del fragmento {previo.numero} por ADN",
                f"Desde el extremo 3'-OH del fragmento {num}, elimina el ARN "
                f"5'-{previo.cebador_arn}-3' y lo rellena con ADN. Queda una "
                "mella (enlace fosfodiéster sin formar) entre ambos fragmentos.",
            )
            mellas.discard(previo.fin)
            previo.ligado = True
            ligaciones += 1
            registrar(
                "Elongación", "ADN ligasa",
                f"Une los fragmentos {previo.numero} y {num}",
                f"Forma el enlace fosfodiéster en la posición {previo.fin} "
                "(consume ATP en eucariotas o NAD+ en E. coli).",
            )

        fin_anterior = fin

    # ------------------------------------------------------------------
    # TERMINACIÓN
    # ------------------------------------------------------------------
    _sustituir_cebador(lider, superior, 0, L)
    registrar(
        "Terminación", "ADN polimerasa I",
        "Sustituye el cebador de la cadena líder",
        "En la célula este hueco lo rellena la horquilla que avanza en "
        "sentido contrario desde el mismo origen (replicación bidireccional).",
    )
    ultimo = fragmentos[-1]
    _sustituir_cebador(rezagada, inferior, ultimo.inicio_cebador, ultimo.fin)
    ultimo.ligado = True
    avisos.append(
        "Problema de la replicación de los extremos: el cebador del último "
        "fragmento de Okazaki no tiene ningún 3'-OH por delante desde el que "
        "rellenar. En un cromosoma lineal quedaría un hueco y el cromosoma se "
        "acortaría; en eucariotas lo compensa la TELOMERASA alargando los "
        "telómeros. En la simulación se rellena para completar la molécula."
    )
    registrar(
        "Terminación", "ADN polimerasa I",
        f"Sustituye el cebador del último fragmento ({ultimo.numero})",
        avisos[-1],
    )
    registrar(
        "Terminación", "Resultado",
        "Replicación semiconservativa completada",
        "Se obtienen dos moléculas de ADN idénticas a la original; cada una "
        "conserva una hebra parental y tiene una hebra nueva.",
    )

    resultado = ResultadoReplicacion(
        hebra_superior=superior,
        hebra_inferior=inferior,
        cadena_lider="".join(lider),
        cadena_rezagada="".join(rezagada),
        cebador_lider="".join(_a_arn(b) for b in superior[:L]).upper(),
        fragmentos=fragmentos,
        eventos=eventos,
        longitud_cebador=L,
        tamano_fragmento=tamano_fragmento,
        puentes_rotos=puentes_rotos,
        ligaciones=ligaciones,
        avisos=avisos,
    )
    if not resultado.es_correcta:  # comprobación interna de coherencia
        raise RuntimeError("Error interno: las moléculas hijas no coinciden.")
    return resultado


def _sustituir_cebador(cadena, molde_o_copia, inicio, fin):
    """Cambia los nucleótidos de ARN de [inicio, fin) por ADN.

    Para la cadena rezagada ``molde_o_copia`` es la hebra inferior (la nueva
    cadena tiene su misma secuencia). Para la líder es la hebra superior.
    """
    for i in range(inicio, fin):
        cadena[i] = molde_o_copia[i]


if __name__ == "__main__":
    # Demostración rápida: python -m dogma.replicacion
    from .secuencias import generar_gen_aleatorio

    adn = generar_gen_aleatorio(n_codones=8, semilla=1)
    r = replicar(adn)
    print("Superior 5'->3':", r.hebra_superior)
    print("Inferior 3'->5':", r.hebra_inferior)
    for e in r.eventos:
        print(f"\n[{e.numero:02}] {e.fase} | {e.enzima}: {e.accion}")
        print("   líder     5'", e.lider, "3'")
        print("   rezagada  3'", e.rezagada, "5'")
    print("\nFragmentos de Okazaki:", [(f.numero, f.inicio + 1, f.fin) for f in r.fragmentos])
    print("¿Moléculas hijas idénticas a la original?", r.es_correcta)
