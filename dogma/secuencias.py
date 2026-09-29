"""Reglas básicas sobre secuencias de ácidos nucleicos.

Este módulo contiene todo lo que usan el resto de etapas del simulador:

* Las reglas de complementariedad de bases:
    - ADN -> ADN (replicación):     A-T, T-A, G-C, C-G
    - ADN -> ARN (transcripción):   A-U, T-A, G-C, C-G
* Limpieza y validación de la secuencia que introduce el usuario.
* Generación de un gen aleatorio "válido" (con codón de inicio y de parada)
  para poder probar el simulador sin tener que escribir una secuencia.

Convención usada en todo el proyecto
-----------------------------------
La secuencia de ADN con la que trabaja el usuario es la **hebra codificante**
(o hebra sentido), escrita en dirección 5' -> 3'. Su hebra complementaria,
antiparalela (3' -> 5'), es la **hebra molde** que usará la ARN polimerasa en
la transcripción. Así, el ARNm resultante tiene la misma secuencia que la
hebra codificante, cambiando T por U.
"""

import random

# --------------------------------------------------------------------------
# Alfabetos
# --------------------------------------------------------------------------
BASES_ADN = "ATGC"
BASES_ARN = "AUGC"

NOMBRES_BASES = {
    "A": "Adenina",
    "T": "Timina",
    "G": "Guanina",
    "C": "Citosina",
    "U": "Uracilo",
}

# --------------------------------------------------------------------------
# Reglas de complementariedad
# --------------------------------------------------------------------------
# ADN -> ADN: la adenina se empareja con la timina (2 puentes de hidrógeno)
# y la guanina con la citosina (3 puentes de hidrógeno).
COMPLEMENTO_ADN = {"A": "T", "T": "A", "G": "C", "C": "G"}

# ADN (hebra molde) -> ARN: igual, pero el ARN usa uracilo en lugar de timina.
COMPLEMENTO_ADN_A_ARN = {"A": "U", "T": "A", "G": "C", "C": "G"}

# ARN -> ARN: usado para calcular los anticodones de los ARNt.
COMPLEMENTO_ARN = {"A": "U", "U": "A", "G": "C", "C": "G"}

PUENTES_HIDROGENO = {"A": 2, "T": 2, "U": 2, "G": 3, "C": 3}

CODON_INICIO = "ATG"                  # AUG en el ARNm (Metionina)
CODONES_STOP = ("TAA", "TAG", "TGA")  # UAA, UAG, UGA en el ARNm


# --------------------------------------------------------------------------
# Limpieza y validación
# --------------------------------------------------------------------------
def limpiar_secuencia(texto):
    """Pasa a mayúsculas y elimina espacios, saltos de línea, números y guiones.

    Permite que el usuario pegue secuencias en formato "bonito", por ejemplo
    ``"atg gcc 5'-TTA-3'"`` o secuencias numeradas tipo GenBank.
    """
    texto = texto.upper().replace("5'", "").replace("3'", "")
    return "".join(c for c in texto if c.isalpha())


def validar_adn(secuencia):
    """Comprueba que la secuencia solo contiene A, T, G y C.

    Lanza ``ValueError`` indicando la primera posición no válida
    (numerada desde 1, como en biología).
    """
    if not secuencia:
        raise ValueError("La secuencia está vacía.")
    for posicion, base in enumerate(secuencia, start=1):
        if base not in BASES_ADN:
            if base == "U":
                raise ValueError(
                    f"Posición {posicion}: 'U' (uracilo) es una base del ARN, "
                    "no del ADN. Introduce la secuencia con T (timina)."
                )
            raise ValueError(
                f"Posición {posicion}: '{base}' no es una base válida del ADN "
                "(solo se admiten A, T, G y C)."
            )
    return True


# --------------------------------------------------------------------------
# Operaciones de complementariedad
# --------------------------------------------------------------------------
def complementaria_adn(secuencia):
    """Devuelve la hebra complementaria base a base.

    Si la entrada va 5' -> 3', la salida se lee alineada debajo, es decir,
    en dirección 3' -> 5' (las dos hebras son antiparalelas).
    """
    return "".join(COMPLEMENTO_ADN[b] for b in secuencia)


def reversa_complementaria(secuencia):
    """Hebra complementaria escrita en su propio sentido 5' -> 3'."""
    return complementaria_adn(secuencia)[::-1]


def complementaria_adn_a_arn(hebra_molde):
    """Aplica las reglas de la transcripción a una hebra molde de ADN."""
    return "".join(COMPLEMENTO_ADN_A_ARN[b] for b in hebra_molde)


def complementaria_arn(secuencia_arn):
    """Complementaria de un ARN (para obtener anticodones)."""
    return "".join(COMPLEMENTO_ARN[b] for b in secuencia_arn)


# --------------------------------------------------------------------------
# Estadísticas sencillas
# --------------------------------------------------------------------------
def contenido_gc(secuencia):
    """Porcentaje de G + C de una secuencia."""
    if not secuencia:
        return 0.0
    gc = sum(1 for b in secuencia if b in "GC")
    return 100 * gc / len(secuencia)


def puentes_hidrogeno(secuencia):
    """Número total de puentes de hidrógeno de la doble hélice."""
    return sum(PUENTES_HIDROGENO[b] for b in secuencia)


def buscar_codon_inicio(secuencia, desde=0):
    """Posición (índice 0) del primer ATG a partir de ``desde``, o -1."""
    return secuencia.find(CODON_INICIO, desde)


# --------------------------------------------------------------------------
# Generación de secuencias de prueba
# --------------------------------------------------------------------------
def generar_gen_aleatorio(n_codones=10, flanco=6, semilla=None):
    """Genera una hebra codificante 5'->3' que contiene un gen completo.

    Estructura:  [flanco 5'] ATG [n_codones aleatorios sin stop] STOP [flanco 3']

    * ``n_codones``: número de codones entre el de inicio y el de parada.
    * ``flanco``: bases aleatorias antes y después del gen (regiones no
      traducidas), para que la búsqueda del AUG tenga sentido.
    * ``semilla``: para obtener siempre la misma secuencia (reproducible).

    El flanco 5' se construye sin ningún ATG para que el primer codón de
    inicio sea el del gen.
    """
    rng = random.Random(semilla)

    def flanco_sin_atg(n):
        while True:
            s = "".join(rng.choice(BASES_ADN) for _ in range(n))
            if CODON_INICIO not in s:
                return s

    def codon_sentido():
        while True:
            c = "".join(rng.choice(BASES_ADN) for _ in range(3))
            if c not in CODONES_STOP:
                return c

    cuerpo = "".join(codon_sentido() for _ in range(n_codones))
    stop = rng.choice(CODONES_STOP)
    inicio_5 = flanco_sin_atg(flanco)
    final_3 = "".join(rng.choice(BASES_ADN) for _ in range(flanco))

    gen = inicio_5 + CODON_INICIO + cuerpo + stop + final_3
    # Si al unir el flanco con el ATG apareciera un ATG antes (p. ej. "..A" + "TG"),
    # se vuelve a generar el flanco.
    while buscar_codon_inicio(gen) != len(inicio_5):
        inicio_5 = flanco_sin_atg(flanco)
        gen = inicio_5 + CODON_INICIO + cuerpo + stop + final_3
    return gen


if __name__ == "__main__":
    # Pequeña demostración: python -m dogma.secuencias
    adn = generar_gen_aleatorio(semilla=1)
    print("Hebra codificante 5'->3':", adn)
    print("Hebra molde       3'->5':", complementaria_adn(adn))
    print(f"Longitud: {len(adn)} pb | GC: {contenido_gc(adn):.1f}% | "
          f"Puentes de H: {puentes_hidrogeno(adn)}")
