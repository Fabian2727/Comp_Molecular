from collections import defaultdict
import math

# Clase para almacenar las probabilidades de emisión
class ProbabilidadesEmision:
    def __init__(self):
        self.estadoCoincidencia = defaultdict(lambda: defaultdict(float))

# Clase para almacenar las probabilidades de transición
class ProbabilidadesTransicion:
    def __init__(self):
        self.transiciones = defaultdict(float)

# Inicializa las probabilidades de transición
def inicializar_probabilidades(transition, num_columnas):
    for i in range(1, num_columnas):
        transition.transiciones[f"M{i}toM{i+1}"] = 0.0
        transition.transiciones[f"M{i}toI{i}"] = 0.0
        transition.transiciones[f"M{i}toD{i}"] = 0.0
        transition.transiciones[f"I{i}toM{i+1}"] = 0.0
        transition.transiciones[f"I{i}toI{i}"] = 0.0
        transition.transiciones[f"D{i}toM{i+1}"] = 0.0
        transition.transiciones[f"D{i}toD{i+1}"] = 0.0

# Selecciona las columnas de interés en las secuencias
def seleccionar_columnas(secuencias):
    columnas_seleccionadas = []
    num_secuencias = len(secuencias)
    longitud_secuencia = len(secuencias[0])
    umbral = num_secuencias // 2

    for i in range(longitud_secuencia):
        conteo_simbolos = sum(1 for seq in secuencias if seq[i] != '-')
        if conteo_simbolos >= umbral:
            columnas_seleccionadas.append(i)

    return columnas_seleccionadas

# Calcula las probabilidades de emisión
def calcular_probabilidades_emision(secuencias, columnas_seleccionadas, emision):
    num_secuencias = len(secuencias)
    tamaño_alfabeto = 20  # 20 aminoácidos
    pseudo_conteo = 1  # Para evitar problemas con frecuencias cero

    for col in columnas_seleccionadas:
        conteos = defaultdict(int)
        conteo_total = 0

        # Contar la ocurrencia de cada símbolo
        for seq in secuencias:
            simbolo = seq[col]
            if simbolo != '-':  # No contar gaps
                conteos[simbolo] += 1
                conteo_total += 1

        # Asignar probabilidades, incluyendo pseudo conteo para evitar frecuencia cero
        for c in "ACDEFGHIKLMNPQRSTVWY":
            probabilidad = (conteos[c] + pseudo_conteo) / float(conteo_total + pseudo_conteo * tamaño_alfabeto)
            emision.estadoCoincidencia[col][c] = probabilidad

# Calcula las probabilidades de transición
def calcular_probabilidades_transicion(secuencias, columnas_seleccionadas, transicion):
    num_secuencias = len(secuencias)
    pseudo_conteo = 1
    conteos_transicion = defaultdict(int)
    conteos_estados = defaultdict(int)

    for i in range(1, len(columnas_seleccionadas)):
        conteos_transicion[f"M{i}toM{i+1}"] = 0
        conteos_transicion[f"M{i}toI{i}"] = 0
        conteos_transicion[f"M{i}toD{i}"] = 0
        conteos_estados[f"M{i}"] = 0

    # Contar transiciones
    for seq_index in range(num_secuencias):
        for i in range(1, len(columnas_seleccionadas)):
            col_index = columnas_seleccionadas[i]
            prev_col_index = columnas_seleccionadas[i - 1]

            if secuencias[seq_index][prev_col_index] != '-':
                conteos_estados[f"M{i}"] += 1
                if secuencias[seq_index][col_index] != '-':
                    conteos_transicion[f"M{i}toM{i+1}"] += 1
                else:
                    conteos_transicion[f"M{i}toD{i}"] += 1
            else:
                conteos_transicion[f"M{i}toI{i}"] += 1

    # Calcular las probabilidades de transición
    for i in range(1, len(columnas_seleccionadas)):
        m_to_m = f"M{i}toM{i+1}"
        m_to_i = f"M{i}toI{i}"
        m_to_d = f"M{i}toD{i}"

        total_transiciones = sum([conteos_transicion[m_to_m], conteos_transicion[m_to_i], conteos_transicion[m_to_d]]) + pseudo_conteo * 3
        transicion.transiciones[m_to_m] = (conteos_transicion[m_to_m] + pseudo_conteo) / total_transiciones
        transicion.transiciones[m_to_i] = (conteos_transicion[m_to_i] + pseudo_conteo) / total_transiciones
        transicion.transiciones[m_to_d] = (conteos_transicion[m_to_d] + pseudo_conteo) / total_transiciones

# Calcula la probabilidad logarítmica de una secuencia
def calcular_log_probabilidad(secuencia, columnas_seleccionadas, emision, transicion):
    log_prob = 0
    estado_anterior = "M1"  # Comienza en el primer estado de coincidencia

    for i, col in enumerate(columnas_seleccionadas):
        if i == 0:
            continue  # El estado inicial ya está definido

        simbolo = secuencia[col]
        if simbolo != '-':
            estado_actual = f"M{i+1}"
            log_prob += math.log(emision.estadoCoincidencia[col].get(simbolo, 1e-6))  # Prob. emisión
            log_prob += math.log(transicion.transiciones.get(f"{estado_anterior}to{estado_actual}", 1e-6))  # Prob. transición
            estado_anterior = estado_actual
        else:
            estado_actual = f"D{i+1}"
            log_prob += math.log(transicion.transiciones.get(f"{estado_anterior}to{estado_actual}", 1e-6))  # Prob. transición
            estado_anterior = estado_actual

    return log_prob

# Verifica si una secuencia pertenece a la familia
def pertenece_a_familia(secuencia, columnas_seleccionadas, emision, transicion, umbral):
    log_prob = calcular_log_probabilidad(secuencia, columnas_seleccionadas, emision, transicion)
    return log_prob >= umbral, log_prob

# Función principal
def HMM(secuencias, secuencia_nueva, umbral):
    emision = ProbabilidadesEmision()
    transicion = ProbabilidadesTransicion()

    columnas_seleccionadas = seleccionar_columnas(secuencias)
    inicializar_probabilidades(transicion, len(columnas_seleccionadas))
    calcular_probabilidades_emision(secuencias, columnas_seleccionadas, emision)
    calcular_probabilidades_transicion(secuencias, columnas_seleccionadas, transicion)

    pertenece, log_prob = pertenece_a_familia(secuencia_nueva, columnas_seleccionadas, emision, transicion, umbral)

    print(f"Log-Probabilidad de la secuencia nueva: {log_prob}")
    print("¿Pertenece a la familia? ", "Sí" if pertenece else "No")

# Ejemplo de uso
if __name__ == "__main__":
    secuencias = [
        "VGA--HAGEY",
        "V----NVDEV",
        "VEA--DVAGH",
        "VKG------D",
        "VYS--TYETS",
        "FNA--NIPKH",
        "IAGADNGAGY"
    ]

    secuencia_nueva = "VGA--AGGFK"
    umbral = -30  # Umbral definido para pertenencia

    HMM(secuencias, secuencia_nueva, umbral)
