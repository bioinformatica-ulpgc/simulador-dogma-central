"""Transcripción (ADN -> ARNm).

La ARN polimerasa usa la **hebra molde** (3' -> 5') como plantilla y sintetiza
una molécula de ARN mensajero en dirección 5' -> 3', aplicando la
complementariedad ADN -> ARN:

    ADN molde:  A   T   G   C
    ARN:        U   A   C   G

Como consecuencia, el ARNm tiene la misma secuencia que la **hebra
codificante** (5' -> 3'), cambiando T por U. El simulador lo comprueba.

Fases representadas:
    1. Iniciación  -> el factor sigma reconoce el promotor, la ARN polimerasa
                      se une y abre la burbuja de transcripción. No necesita
                      cebador (a diferencia de la ADN polimerasa).
    2. Elongación  -> la burbuja (~12-14 pb en la célula) avanza, se añaden
                      ribonucleótidos (NTP) y la hélice se vuelve a cerrar
                      por detrás.
    3. Terminación -> la polimerasa encuentra la señal de terminación y
                      libera el ARNm.

La unidad de transcripción es la secuencia completa introducida.
"""

from dataclasses import dataclass, field

from .secuencias import (
    COMPLEMENTO_ADN_A_ARN,
    complementaria_adn,
    complementaria_adn_a_arn,
    validar_adn,
)

VACIO = "·"

ENZIMAS = {
    "Factor sigma (σ)": (
        "σ70 en E. coli (factores generales TFIID/TBP, TFIIB... en eucariotas)",
        "Reconoce el promotor (cajas -35 TTGACA y -10 TATAAT; caja TATA en "
        "eucariotas) y coloca a la ARN polimerasa en el inicio.",
    ),
    "ARN polimerasa": (
        "Núcleo α2ββ'ω (ARN polimerasa II para ARNm en eucariotas)",
        "Abre la doble hélice, lee la hebra molde 3'->5' y une ribonucleótidos "
        "(ATP, UTP, GTP, CTP) en dirección 5'->3'. No necesita cebador.",
    ),
    "Terminador": (
        "Horquilla rica en GC + cola de U, o factor Rho (señal poli-A en eucariotas)",
        "Señal que hace que la ARN polimerasa se detenga y libere el ARNm.",
    ),
}


@dataclass
class EventoTranscripcion:
    numero: int
    fase: str
    enzima: str
    accion: str
    detalle: str
    burbuja: tuple   # (inicio, fin) de la región abierta; (0, 0) si cerrada
    arn: str         # ARN sintetizado hasta el momento, alineado con el ADN


@dataclass
class ResultadoTranscripcion:
    hebra_codificante: str   # 5' -> 3'
    hebra_molde: str         # 3' -> 5' (alineada debajo de la codificante)
    arnm: str                # 5' -> 3'
    eventos: list
    tamano_burbuja: int
    avisos: list = field(default_factory=list)

    @property
    def enlaces_fosfodiester(self):
        return max(0, len(self.arnm) - 1)

    @property
    def composicion(self):
        return {b: self.arnm.count(b) for b in "AUGC"}

    @property
    def es_correcta(self):
        """El ARNm coincide con la hebra codificante cambiando T por U."""
        return self.arnm == self.hebra_codificante.replace("T", "U")


def transcribir(adn, paso=None, tamano_burbuja=12):
    """Simula la transcripción de ``adn`` (hebra codificante 5'->3').

    * ``paso``: nucleótidos añadidos en cada evento de elongación (para la
      visualización). Si es ``None`` se ajusta para obtener unos 6 pasos.
    * ``tamano_burbuja``: pares de bases abiertos alrededor de la polimerasa.
    """
    validar_adn(adn)
    n = len(adn)
    if paso is None:
        paso = max(3, -(-n // 6))
    if paso < 1:
        raise ValueError("El paso de elongación debe ser al menos 1.")

    codificante = adn
    molde = complementaria_adn(adn)
    arn = [VACIO] * n
    eventos = []

    def burbuja_en(pos):
        """Región abierta cuando la polimerasa está en ``pos``."""
        mitad = tamano_burbuja // 2
        return (max(0, pos - mitad), min(n, pos + mitad))

    def registrar(fase, enzima, accion, detalle, burbuja):
        eventos.append(
            EventoTranscripcion(
                numero=len(eventos) + 1,
                fase=fase,
                enzima=enzima,
                accion=accion,
                detalle=detalle,
                burbuja=burbuja,
                arn="".join(arn),
            )
        )

    # ------------------------------------------------------------------
    # INICIACIÓN
    # ------------------------------------------------------------------
    registrar(
        "Iniciación", "Factor sigma (σ)",
        "Reconoce el promotor",
        "El factor σ se une a las secuencias del promotor situadas antes del "
        "inicio (-35 y -10) y recluta a la ARN polimerasa. Se usa como molde "
        "la hebra inferior (3'->5'); la superior es la hebra codificante.",
        (0, 0),
    )
    registrar(
        "Iniciación", "ARN polimerasa",
        "Abre la burbuja de transcripción",
        f"Separa unos {min(tamano_burbuja, n)} pb de la doble hélice para "
        "dejar expuesta la hebra molde. No necesita cebador: puede empezar "
        "una cadena de ARN desde cero.",
        (0, min(n, tamano_burbuja)),
    )

    # ------------------------------------------------------------------
    # ELONGACIÓN
    # ------------------------------------------------------------------
    pos = 0
    primero = True
    while pos < n:
        fin = min(n, pos + paso)
        for i in range(pos, fin):
            arn[i] = COMPLEMENTO_ADN_A_ARN[molde[i]]
        tramo_molde = molde[pos:fin]
        tramo_arn = "".join(arn[pos:fin])
        detalle = (
            f"Lee el molde 3'-{tramo_molde}-5' y añade 5'-{tramo_arn}-3' "
            f"(posiciones {pos + 1}-{fin}). "
        )
        if primero:
            detalle += (
                "Tras los primeros nucleótidos el factor σ se libera y la "
                "polimerasa pasa a la fase de elongación. "
            )
            primero = False
        detalle += "Por detrás, la doble hélice de ADN se vuelve a cerrar."
        registrar(
            "Elongación", "ARN polimerasa",
            f"Añade {fin - pos} ribonucleótidos (5'->3')",
            detalle,
            burbuja_en(fin),
        )
        pos = fin

    # ------------------------------------------------------------------
    # TERMINACIÓN
    # ------------------------------------------------------------------
    registrar(
        "Terminación", "Terminador",
        "La ARN polimerasa llega a la señal de terminación",
        "En procariotas, una horquilla rica en GC seguida de U en el ARN (o el "
        "factor Rho) desprende la polimerasa. En eucariotas, el corte ocurre "
        "tras la señal de poliadenilación AAUAAA.",
        burbuja_en(n),
    )
    registrar(
        "Terminación", "ARN polimerasa",
        "Libera el ARNm y el ADN recupera la doble hélice",
        f"Se obtiene un ARNm de {n} nucleótidos unidos por {max(0, n - 1)} "
        "enlaces fosfodiéster.",
        (0, 0),
    )

    avisos = [
        "En eucariotas el ARN recién sintetizado (pre-ARNm) se procesa en el "
        "núcleo antes de salir al citoplasma: se añade la caperuza 5' (7-metil-"
        "guanosina), la cola poli-A en el extremo 3' y se eliminan los intrones "
        "(splicing). La simulación trabaja con la secuencia codificante ya "
        "madura, como ocurre directamente en procariotas."
    ]

    arnm = "".join(arn)
    resultado = ResultadoTranscripcion(
        hebra_codificante=codificante,
        hebra_molde=molde,
        arnm=arnm,
        eventos=eventos,
        tamano_burbuja=tamano_burbuja,
        avisos=avisos,
    )
    # Comprobaciones de coherencia: dos formas distintas de obtener el ARNm
    if arnm != complementaria_adn_a_arn(molde) or not resultado.es_correcta:
        raise RuntimeError("Error interno: el ARNm no es coherente.")
    return resultado


if __name__ == "__main__":
    # Demostración rápida: python -m dogma.transcripcion
    from .secuencias import generar_gen_aleatorio

    adn = generar_gen_aleatorio(n_codones=8, semilla=1)
    r = transcribir(adn)
    print("Codificante 5'->3':", r.hebra_codificante)
    print("Molde       3'->5':", r.hebra_molde)
    for e in r.eventos:
        print(f"\n[{e.numero:02}] {e.fase} | {e.enzima}: {e.accion}  burbuja={e.burbuja}")
        print("   ARN      5'", e.arn, "3'")
    print("\nARNm 5'->3':", r.arnm)
    print("Composición:", r.composicion)
    print("¿ARNm = codificante con U en lugar de T?", r.es_correcta)
