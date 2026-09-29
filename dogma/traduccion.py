"""Traducción y síntesis de proteínas (ARNm -> proteína).

El ribosoma lee el ARNm en dirección 5' -> 3' de tres en tres nucleótidos
(codones). Cada codón se empareja con el **anticodón** complementario y
antiparalelo de un ARN de transferencia (ARNt) que transporta el aminoácido
correspondiente según el **código genético**.

    Codón (ARNm)    5'- A U G -3'
    Anticodón (ARNt) 3'- U A C -5'     ->  Metionina (Met, M)

Fases representadas:
    1. Iniciación  -> la subunidad pequeña se une al ARNm y busca el primer
                      AUG; entra el ARNt iniciador en el sitio P y se une la
                      subunidad grande.
    2. Elongación  -> por cada codón: entra el aminoacil-ARNt en el sitio A,
                      la peptidil transferasa forma el enlace peptídico y el
                      ribosoma se transloca un codón (el ARNt vacío sale por E).
    3. Terminación -> un codón de parada (UAA, UAG, UGA) entra en el sitio A,
                      un factor de liberación libera la proteína y el ribosoma
                      se disocia.
"""

from dataclasses import dataclass, field

from .secuencias import complementaria_arn

# --------------------------------------------------------------------------
# Código genético estándar
# --------------------------------------------------------------------------
_BASES = "UCAG"
_AMINOACIDOS_64 = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"

CODIGO_GENETICO = {
    b1 + b2 + b3: _AMINOACIDOS_64[16 * i + 4 * j + k]
    for i, b1 in enumerate(_BASES)
    for j, b2 in enumerate(_BASES)
    for k, b3 in enumerate(_BASES)
}

AMINOACIDOS = {
    # letra: (abreviatura de 3 letras, nombre)
    "A": ("Ala", "Alanina"),
    "R": ("Arg", "Arginina"),
    "N": ("Asn", "Asparagina"),
    "D": ("Asp", "Ácido aspártico"),
    "C": ("Cys", "Cisteína"),
    "Q": ("Gln", "Glutamina"),
    "E": ("Glu", "Ácido glutámico"),
    "G": ("Gly", "Glicina"),
    "H": ("His", "Histidina"),
    "I": ("Ile", "Isoleucina"),
    "L": ("Leu", "Leucina"),
    "K": ("Lys", "Lisina"),
    "M": ("Met", "Metionina"),
    "F": ("Phe", "Fenilalanina"),
    "P": ("Pro", "Prolina"),
    "S": ("Ser", "Serina"),
    "T": ("Thr", "Treonina"),
    "W": ("Trp", "Triptófano"),
    "Y": ("Tyr", "Tirosina"),
    "V": ("Val", "Valina"),
    "*": ("Stop", "Codón de parada"),
}

CODON_INICIO = "AUG"
CODONES_STOP = ("UAA", "UAG", "UGA")

ENZIMAS = {
    "Ribosoma": (
        "70S = 30S + 50S (80S = 40S + 60S en eucariotas)",
        "Formado por ARNr y proteínas. Tiene tres sitios: A (entra el "
        "aminoacil-ARNt), P (ARNt con la cadena peptídica) y E (salida).",
    ),
    "ARNt": (
        "ARN de transferencia",
        "Lleva un aminoácido en su extremo 3' (CCA) y reconoce el codón con "
        "su anticodón mediante complementariedad de bases.",
    ),
    "Aminoacil-ARNt sintetasa": (
        "Una por aminoácido",
        "Une cada aminoácido a su ARNt correcto gastando ATP ('carga' el ARNt).",
    ),
    "Factores de iniciación": (
        "IF1, IF2, IF3 (eIF en eucariotas)",
        "Ayudan a montar el ribosoma sobre el AUG con el ARNt iniciador.",
    ),
    "Factores de elongación": (
        "EF-Tu y EF-G (eEF1A y eEF2 en eucariotas)",
        "EF-Tu lleva el aminoacil-ARNt al sitio A y EF-G mueve el ribosoma "
        "un codón (translocación). Cada uno gasta 1 GTP.",
    ),
    "Peptidil transferasa": (
        "Actividad del ARNr 23S (28S en eucariotas): es una ribozima",
        "Forma el enlace peptídico entre el aminoácido del sitio A y la "
        "cadena del sitio P.",
    ),
    "Factores de liberación": (
        "RF1/RF2 y RF3 (eRF1 y eRF3 en eucariotas)",
        "Reconocen el codón de parada y liberan la proteína.",
    ),
}


def aminoacido(codon):
    """Letra del aminoácido codificado por un codón de ARN."""
    return CODIGO_GENETICO[codon]


def anticodon(codon):
    """Anticodón del ARNt, alineado con el codón (se lee 3' -> 5')."""
    return complementaria_arn(codon)


@dataclass
class EventoTraduccion:
    numero: int
    fase: str
    enzima: str
    accion: str
    detalle: str
    indice_codon: int     # nº de codón (desde 0) que ocupa el sitio A; -1 si no hay
    sitio_e: str
    sitio_p: str
    sitio_a: str
    peptido: str          # cadena de aminoácidos (1 letra) hasta el momento


@dataclass
class ResultadoTraduccion:
    arnm: str
    inicio: int           # índice del AUG en el ARNm, -1 si no hay
    fin: int              # índice justo después del codón stop (o del último codón leído)
    codones: list         # codones leídos desde el AUG, stop incluido si existe
    proteina: str         # 1 letra, sin el stop
    tiene_stop: bool
    eventos: list
    avisos: list = field(default_factory=list)

    @property
    def utr5(self):
        return self.arnm[: self.inicio] if self.inicio >= 0 else self.arnm

    @property
    def utr3(self):
        return self.arnm[self.fin:] if self.inicio >= 0 else ""

    @property
    def proteina_3_letras(self):
        return "-".join(AMINOACIDOS[a][0] for a in self.proteina)

    @property
    def enlaces_peptidicos(self):
        return max(0, len(self.proteina) - 1)

    @property
    def gtp_consumidos(self):
        """1 (iniciación) + 2 por cada elongación + 1 (terminación, RF3)."""
        if not self.proteina:
            return 0
        return 1 + 2 * self.enlaces_peptidicos + (1 if self.tiene_stop else 0)


def traducir(arnm):
    """Simula la traducción de un ARNm (5'->3') empezando en el primer AUG."""
    arnm = arnm.upper().replace("T", "U")
    for pos, b in enumerate(arnm, start=1):
        if b not in "AUGC":
            raise ValueError(f"Posición {pos}: '{b}' no es un nucleótido de ARN.")

    eventos = []
    avisos = []
    peptido = []
    sitios = {"E": "", "P": "", "A": ""}

    def registrar(fase, enzima, accion, detalle, indice=-1):
        eventos.append(
            EventoTraduccion(
                numero=len(eventos) + 1,
                fase=fase,
                enzima=enzima,
                accion=accion,
                detalle=detalle,
                indice_codon=indice,
                sitio_e=sitios["E"],
                sitio_p=sitios["P"],
                sitio_a=sitios["A"],
                peptido="".join(peptido),
            )
        )

    inicio = arnm.find(CODON_INICIO)
    if inicio == -1:
        avisos.append(
            "El ARNm no contiene ningún codón de inicio AUG: el ribosoma no "
            "puede empezar la traducción y no se sintetiza proteína."
        )
        registrar(
            "Iniciación", "Ribosoma",
            "No se encuentra el codón de inicio AUG",
            avisos[-1],
        )
        return ResultadoTraduccion(arnm, -1, -1, [], "", False, eventos, avisos)

    # Codones en el marco de lectura que define el AUG
    codones = []
    tiene_stop = False
    for i in range(inicio, len(arnm) - 2, 3):
        c = arnm[i:i + 3]
        codones.append(c)
        if c in CODONES_STOP:
            tiene_stop = True
            break
    fin = inicio + 3 * len(codones)

    # ------------------------------------------------------------------
    # INICIACIÓN
    # ------------------------------------------------------------------
    registrar(
        "Iniciación", "Factores de iniciación",
        "La subunidad pequeña se une al ARNm",
        "En procariotas se une a la secuencia Shine-Dalgarno; en eucariotas "
        "se une a la caperuza 5' y recorre el ARNm (escaneo) buscando el "
        f"primer AUG. La región 5' no traducida (5'-UTR) tiene {inicio} nt.",
    )
    sitios["P"] = CODON_INICIO
    peptido.append("M")
    registrar(
        "Iniciación", "ARNt",
        f"El ARNt iniciador reconoce el AUG (posición {inicio + 1})",
        f"Anticodón 3'-{anticodon(CODON_INICIO)}-5' emparejado con el codón "
        "5'-AUG-3'. Lleva Metionina (formil-metionina, fMet, en procariotas) "
        "y se coloca directamente en el sitio P. Se gasta 1 GTP (IF2).",
        0,
    )
    registrar(
        "Iniciación", "Ribosoma",
        "Se une la subunidad grande: ribosoma completo",
        "Queda el sitio P ocupado por el ARNt iniciador y el sitio A libre "
        "para el siguiente codón. El marco de lectura queda fijado por el AUG.",
        0,
    )

    # ------------------------------------------------------------------
    # ELONGACIÓN
    # ------------------------------------------------------------------
    for k in range(1, len(codones)):
        c = codones[k]
        if c in CODONES_STOP:
            break
        aa = aminoacido(c)
        abrev, nombre = AMINOACIDOS[aa]
        sitios["A"] = c
        registrar(
            "Elongación", "Factores de elongación",
            f"Codón {k + 1}: {c} -> entra ARNt con {nombre} ({abrev})",
            f"EF-Tu lleva al sitio A el aminoacil-ARNt con anticodón "
            f"3'-{anticodon(c)}-5', complementario del codón 5'-{c}-3' "
            "(A-U, G-C). Gasta 1 GTP.",
            k,
        )
        peptido.append(aa)
        registrar(
            "Elongación", "Peptidil transferasa",
            f"Enlace peptídico nº {len(peptido) - 1}",
            f"La cadena del sitio P se transfiere al {nombre} del sitio A. "
            f"Péptido: {'-'.join(AMINOACIDOS[x][0] for x in peptido)}.",
            k,
        )
        sitios["E"], sitios["P"], sitios["A"] = sitios["P"], c, ""
        registrar(
            "Elongación", "Factores de elongación",
            "Translocación del ribosoma un codón (5'->3')",
            f"EF-G mueve el ribosoma: el ARNt con la cadena pasa al sitio P y "
            f"el ARNt vacío ({sitios['E']}) pasa al sitio E y sale. Gasta 1 GTP.",
            k,
        )

    # ------------------------------------------------------------------
    # TERMINACIÓN
    # ------------------------------------------------------------------
    proteina = "".join(peptido)
    if tiene_stop:
        stop = codones[-1]
        sitios["E"], sitios["A"] = "", stop
        registrar(
            "Terminación", "Factores de liberación",
            f"Codón de parada {stop} en el sitio A",
            f"Ningún ARNt reconoce {stop}. Un factor de liberación (RF1/RF2, "
            "eRF1) ocupa el sitio A y hace que se hidrolice el enlace entre la "
            "proteína y el último ARNt.",
            len(codones) - 1,
        )
        sitios["P"], sitios["A"] = "", ""
        registrar(
            "Terminación", "Ribosoma",
            "Se libera la proteína y el ribosoma se disocia",
            f"Proteína de {len(proteina)} aminoácidos unidos por "
            f"{len(proteina) - 1} enlaces peptídicos. Extremo amino (N) = Met "
            "inicial; extremo carboxilo (C) = último aminoácido.",
        )
    else:
        avisos.append(
            "No se ha encontrado ningún codón de parada en el marco de lectura: "
            "el ribosoma llega al final del ARNm sin terminar. La proteína "
            "quedaría incompleta (en bacterias la rescata el ARNtm)."
        )
        registrar(
            "Terminación", "Ribosoma",
            "Fin del ARNm sin codón de parada",
            avisos[-1],
        )

    return ResultadoTraduccion(
        arnm=arnm,
        inicio=inicio,
        fin=fin,
        codones=codones,
        proteina=proteina,
        tiene_stop=tiene_stop,
        eventos=eventos,
        avisos=avisos,
    )


if __name__ == "__main__":
    # Demostración rápida: python -m dogma.traduccion
    from .secuencias import generar_gen_aleatorio

    arnm = generar_gen_aleatorio(n_codones=8, semilla=1).replace("T", "U")
    r = traducir(arnm)
    print("ARNm 5'->3':", r.arnm)
    print("5'-UTR:", r.utr5, "| codones:", " ".join(r.codones), "| 3'-UTR:", r.utr3)
    for e in r.eventos:
        print(f"[{e.numero:02}] {e.fase} | {e.enzima}: {e.accion}"
              f"   E={e.sitio_e or '---'} P={e.sitio_p or '---'} A={e.sitio_a or '---'}"
              f"   péptido={e.peptido}")
    print("\nProteína (1 letra):", r.proteina)
    print("Proteína (3 letras):", r.proteina_3_letras)
    print("Enlaces peptídicos:", r.enlaces_peptidicos, "| GTP:", r.gtp_consumidos)
