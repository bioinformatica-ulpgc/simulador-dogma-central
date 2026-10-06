"""Simulador del dogma central de la biología molecular.

Uso interactivo (menú):
    python main.py

Uso directo desde la línea de comandos (ejemplos):
    python main.py --ejemplo
    python main.py --secuencia ATGGCTAGCAAATAA --modo continuo
    python main.py --aleatoria 12 --semilla 7 --modo resumen
    python main.py --archivo ejemplos/egfp.fasta
    python main.py --archivo ejemplos/lacZ_Ecoli_NC_000913.3.gb --modo resumen
    python main.py --ejemplo --modo resumen --sin-imagenes
"""

import argparse
import os
import sys

from dogma import consola as ui
from dogma.replicacion import replicar
from dogma.secuencias import (
    contenido_gc,
    generar_gen_aleatorio,
    limpiar_secuencia,
    validar_adn,
)
from dogma.traduccion import traducir
from dogma.transcripcion import transcribir

# Inicio del gen de la proteína verde fluorescente (EGFP) con una región 5'
# (secuencia Kozak GCCACC) y un codón de parada añadido al final.
SECUENCIA_EJEMPLO = "GCCACC" "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGG" "TAA" "GCGGCC"
DESCRIPCION_EJEMPLO = ("Inicio del gen EGFP (proteína verde fluorescente) con "
                       "secuencia Kozak y codón de parada añadido")

MODOS = {"paso": "Paso a paso (Enter entre pasos)",
         "continuo": "Todo seguido, sin pausas",
         "resumen": "Solo resultados de cada etapa"}


# --------------------------------------------------------------------------
# Simulación completa
# --------------------------------------------------------------------------
def generar_y_mostrar_imagenes(rep, tra, trad, carpeta="resultados", abrir=None):
    """Genera las figuras PNG y, si se pide, abre la carpeta en Windows."""
    try:
        from dogma.visualizacion import generar_imagenes
    except ImportError:
        ui.aviso("No se pueden generar imágenes: falta matplotlib "
                 "(pip install matplotlib).")
        return []
    ui.subtitulo("Imágenes")
    rutas = generar_imagenes(rep, tra, trad, carpeta)
    for r in rutas:
        ui.ok(os.path.abspath(r))
    if abrir is None and hasattr(os, "startfile"):
        abrir = input("   ¿Abrir la carpeta de imágenes? [s/N]: ").strip().lower() == "s"
    if abrir and hasattr(os, "startfile"):
        os.startfile(os.path.abspath(carpeta))
    return rutas


def simular(adn, modo="paso", descripcion=None, imagenes=None, carpeta="resultados"):
    """Ejecuta replicación, transcripción y traducción y las muestra.

    ``imagenes``: True genera las figuras, False no, None pregunta al usuario.
    """
    validar_adn(adn)
    ui.titulo("Simulador del dogma central de la biología molecular")
    if descripcion:
        ui.parrafo(descripcion)
    print(f"   Secuencia de ADN (hebra codificante 5'->3'): {len(adn)} pb, "
          f"GC = {contenido_gc(adn):.1f}%")
    ui.imprimir_alineado([("ADN", "5'", adn, "3'", "seq")], len(adn))

    rep = replicar(adn)
    tra = transcribir(rep.cadena_lider)   # se transcribe una molécula hija
    trad = traducir(tra.arnm)

    modo = ui.mostrar_replicacion(rep, modo)
    modo = ui.pausa(modo)
    modo = ui.mostrar_transcripcion(tra, modo)
    modo = ui.pausa(modo)
    ui.mostrar_traduccion(trad, modo)
    ui.mostrar_resumen(rep, tra, trad)

    if imagenes is None:
        imagenes = input("\n   ¿Generar las imágenes de cada etapa (PNG)? [S/n]: "
                         ).strip().lower() != "n"
        abrir = None
    else:
        abrir = False
    if imagenes:
        generar_y_mostrar_imagenes(rep, tra, trad, carpeta, abrir)
    return rep, tra, trad


# --------------------------------------------------------------------------
# Entrada de datos
# --------------------------------------------------------------------------
def leer_archivo(ruta):
    """Lee una secuencia de un archivo de texto, FASTA o GenBank.

    - FASTA / texto: se ignoran las líneas de cabecera ('>' o ';').
    - GenBank (.gb): la secuencia está entre la línea ORIGIN y '//'.
    """
    with open(ruta, encoding="utf-8") as f:
        lineas = f.readlines()
    if any(l.startswith("ORIGIN") for l in lineas):
        inicio = next(i for i, l in enumerate(lineas) if l.startswith("ORIGIN"))
        lineas = [l for l in lineas[inicio + 1:] if not l.startswith("//")]
    else:
        lineas = [l for l in lineas if not l.startswith(">") and not l.startswith(";")]
    return limpiar_secuencia("".join(lineas))


def pedir_entero(mensaje, defecto, minimo=1, maximo=None):
    while True:
        r = input(f"   {mensaje} [{defecto}]: ").strip()
        if not r:
            return defecto
        try:
            v = int(r)
            if v < minimo or (maximo is not None and v > maximo):
                raise ValueError
            return v
        except ValueError:
            rango = f"entre {minimo} y {maximo}" if maximo else f"mayor o igual que {minimo}"
            print(f"   Introduce un número {rango}.")


def pedir_secuencia_manual():
    print("   Escribe o pega la hebra codificante de ADN en dirección 5'->3'.")
    print("   (Se ignoran espacios, números y minúsculas. Deja vacío para volver.)")
    while True:
        texto = input("   ADN > ")
        if not texto.strip():
            return None
        adn = limpiar_secuencia(texto)
        try:
            validar_adn(adn)
            if len(adn) < 6:
                raise ValueError("La secuencia debe tener al menos 6 pb.")
            return adn
        except ValueError as e:
            print(ui.color(f"   Error: {e}", ui.Fore.RED))


def pedir_modo():
    print()
    claves = list(MODOS)
    for i, k in enumerate(claves, start=1):
        print(f"   {i}. {MODOS[k]}")
    r = input("   Modo de visualización [1]: ").strip()
    try:
        return claves[int(r) - 1]
    except (ValueError, IndexError):
        return "paso"


# --------------------------------------------------------------------------
# Menú principal
# --------------------------------------------------------------------------
def menu():
    while True:
        ui.titulo("Simulador del dogma central de la biología molecular")
        print("   ADN --(replicación)--> ADN --(transcripción)--> ARNm "
              "--(traducción)--> proteína\n")
        print("   1. Introducir una secuencia de ADN")
        print("   2. Generar un gen aleatorio")
        print("   3. Usar la secuencia de ejemplo (inicio del gen EGFP)")
        print("   4. Cargar una secuencia desde un archivo (.txt / FASTA / GenBank)")
        print("   5. Ver las enzimas y moléculas que intervienen")
        print("   6. Ver el código genético")
        print("   0. Salir")
        opcion = input("\n   Elige una opción: ").strip()

        adn, descripcion = None, None
        if opcion == "1":
            adn = pedir_secuencia_manual()
        elif opcion == "2":
            n = pedir_entero("Número de codones entre el AUG y el stop", 8, 1, 200)
            semilla_txt = input("   Semilla (vacío = aleatoria): ").strip()
            semilla = int(semilla_txt) if semilla_txt.lstrip("-").isdigit() else None
            adn = generar_gen_aleatorio(n_codones=n, semilla=semilla)
            descripcion = f"Gen aleatorio de {n} codones (semilla: {semilla})"
        elif opcion == "3":
            adn, descripcion = SECUENCIA_EJEMPLO, DESCRIPCION_EJEMPLO
        elif opcion == "4":
            ruta = input("   Ruta del archivo: ").strip().strip('"')
            try:
                adn = leer_archivo(ruta)
                validar_adn(adn)
                descripcion = f"Secuencia cargada de {ruta}"
            except (OSError, ValueError) as e:
                print(ui.color(f"   Error: {e}", ui.Fore.RED))
                adn = None
        elif opcion == "5":
            ui.mostrar_enzimas()
        elif opcion == "6":
            ui.mostrar_codigo_genetico()
        elif opcion == "0":
            print("   ¡Hasta luego!")
            return
        else:
            print("   Opción no válida.")

        if adn:
            if len(adn) > 300:
                ui.aviso(f"La secuencia tiene {len(adn)} pb: el modo paso a paso "
                         "será muy largo. Se recomienda el modo resumen.")
            try:
                simular(adn, pedir_modo(), descripcion)
            except ValueError as e:
                print(ui.color(f"   Error: {e}", ui.Fore.RED))
        input(ui.color("\n   Pulsa Enter para volver al menú...", ui.Style.DIM))


def main(argv=None):
    # Asegura que los caracteres especiales (▲, σ, ═...) se ven bien en Windows
    for flujo in (sys.stdout, sys.stderr):
        if hasattr(flujo, "reconfigure"):
            flujo.reconfigure(encoding="utf-8")

    p = argparse.ArgumentParser(
        description="Simulador del dogma central de la biología molecular "
                    "(replicación, transcripción y traducción).")
    g = p.add_mutually_exclusive_group()
    g.add_argument("--secuencia", help="hebra codificante de ADN (5'->3')")
    g.add_argument("--aleatoria", type=int, metavar="N",
                   help="genera un gen aleatorio con N codones")
    g.add_argument("--ejemplo", action="store_true",
                   help="usa la secuencia de ejemplo (inicio del gen EGFP)")
    g.add_argument("--archivo", help="lee la secuencia de un archivo .txt, FASTA o GenBank")
    p.add_argument("--semilla", type=int, help="semilla para --aleatoria")
    p.add_argument("--modo", choices=MODOS, default="paso",
                   help="paso (por defecto), continuo o resumen")
    p.add_argument("--sin-imagenes", action="store_true",
                   help="no genera las imágenes PNG")
    p.add_argument("--carpeta", default="resultados",
                   help="carpeta donde se guardan las imágenes (por defecto: resultados)")
    args = p.parse_args(argv)
    img = {"imagenes": not args.sin_imagenes, "carpeta": args.carpeta}

    try:
        if args.secuencia:
            simular(limpiar_secuencia(args.secuencia), args.modo, **img)
        elif args.aleatoria:
            simular(generar_gen_aleatorio(args.aleatoria, semilla=args.semilla),
                    args.modo, f"Gen aleatorio de {args.aleatoria} codones", **img)
        elif args.ejemplo:
            simular(SECUENCIA_EJEMPLO, args.modo, DESCRIPCION_EJEMPLO, **img)
        elif args.archivo:
            simular(leer_archivo(args.archivo), args.modo, f"Archivo: {args.archivo}",
                    **img)
        else:
            menu()
    except ValueError as e:
        print(ui.color(f"Error: {e}", ui.Fore.RED))
        sys.exit(1)
    except (KeyboardInterrupt, EOFError):
        print("\n   Simulación interrumpida.")


if __name__ == "__main__":
    main()
