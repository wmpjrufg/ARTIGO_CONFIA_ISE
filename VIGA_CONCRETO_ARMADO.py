################################################################################
# UNIVERSIDADE FEDERAL DE CATALÃO (UFCAT)
# WANDERLEI MALAQUIAS PEREIRA JUNIOR,                 ENG. CIVIL / PROF. (UFCAT)
# ANDRÉ TEÓFILO BECK,                                   ENG. MEC. / PROF. (EESC)
# DANIEL LIMA ARAÚJO,                                   ENG. CIVIL / PROF. (UFG)
# MAURO SOUSA,                                         ENG. CIVIL / PROF. (UFPA)
# MATHEUS HENRIQUE MORATO DE MORAES,                          ENG. CIVIL (UFCAT)
################################################################################

################################################################################
# DESCRIÇÃO ALGORITMO:
# BIBLIOTECA DE VERIFICAÇÃO DE ESTADO LIMITE EM VIGAS DE CONCRETO ARMADO DESEN-
# VOLVIDA PELO GRUPO DE PESQUISAS E ESTUDOS EM ENGENHARIA (GPEE)
################################################################################

################################################################################
# BIBLIOTECAS PYTHON
import numpy as np

################################################################################
# BIBLIOTECAS DESENVOLVEDORES GPEE

def MODULO_ELASTICIDADE_CONCRETO(AGREGADO, F_CK, F_CKJ):
    """
    Esta função calcula os módulos de elasticidade do concreto.  

    Entrada:
    AGREGADO    | Tipo de agragado usado no traço do cimento       |        | string    
                |   'BAS' - Agregado de Basalto                    |        | 
                |   'GRA' - Agregado de Granito                    |        |              
                |   'CAL' - Agregado de Calcário                   |        |
                |   'ARE' - Agregado de Arenito                    |        | 
    F_CK        | Resistência característica à compressão          | kN/m²  | float   
    F_CKJ       | Resistência característica à compressão idade J  | kN/m²  | float
    
    Saída:
    E_CIJ       | Módulo de elasticidade tangente                  | kN/m²  | float
    E_CSJ       | Módulo de elasticidade do secante                | kN/m²  | float   
    """
    # Determinação do módulo tangente E_CI idade T
    if AGREGADO == 'BAS':         
        ALFA_E = 1.2
    elif AGREGADO == 'GRA':         
        ALFA_E = 1.0
    elif AGREGADO == 'CAL':       
        ALFA_E = 0.9
    elif AGREGADO == 'ARE':       
        ALFA_E = 0.7
    F_CK /= 1E3
    if F_CK <= 50:        
        E_CI = ALFA_E * 5600 * np.sqrt(F_CK)
    elif F_CK > 50:   
        E_CI = 21.5 * (10 ** 3) * ALFA_E * (F_CK / 10 + 1.25) ** (1 / 3)
    ALFA_I = 0.8 + 0.2 * F_CK / 80
    if ALFA_I > 1:        
        ALFA_I = 1
    # Determinação do módulo secante E_CS idade T
    E_CS = E_CI * ALFA_I
    if F_CK <= 45 :
        F_CK *= 1E3
        E_CIJ = E_CI * (F_CKJ / F_CK) ** 0.5  
    elif  F_CK > 45 : 
        F_CK *= 1E3
        E_CIJ = E_CI * (F_CKJ / F_CK) ** 0.3  
    E_CSJ = E_CIJ * ALFA_I
    E_CIJ *= 1E3 
    E_CSJ *= 1E3 
    return E_CIJ, E_CSJ
    
def MOMENTO_RESISTENE_MRD(B_W, D, A_S, F_Y, F_CK, GAMMA_M_ACO, GAMMA_M_CONCRETO):
    """
    Esta função determina o momento resistente M_RD de uma viga de concreto armado.

    Entrada:
    B_W                | Largura da viga                               | m     | float  
    D                  | Altura útil da seção                          | m     | float
    A_S                | Área de aço necessária na seção               | m²    | float
    F_Y                | Tensão de escoamento do aço                   | kN/m² | float
    GAMMA_M_ACO        | Coeficiente parcial de segurança do aço       |       | float
    GAMMA_M_CONCRETO   | Coeficiente parcial de segurança do concreto  |       | float

    Saída:
    M_RD               | Momento resistente da peça                    | kN.m  | Float 

    """
    # Definição dos parâmetros que dependem do F_CK
    F_CK /= 1E3
    if F_CK >  50:
        LAMBDA = 0.80 - ((F_CK - 50) / 400)
        ALPHA_C = (1.00 - ((F_CK - 50) / 200)) * 0.85
    else:
        LAMBDA = 0.80
        ALPHA_C = 0.85
    F_CK *= 1E3    
    F_CD = F_CK / GAMMA_M_CONCRETO
    # Linha neutra real e momento resistente
    X_III = (A_S * (F_Y / GAMMA_M_ACO)) / (ALPHA_C * F_CD * B_W * LAMBDA)
    M_RD = A_S * (F_Y / GAMMA_M_ACO) * (D - 0.50 * LAMBDA * X_III)
    return M_RD

def CORTANTE_RESISTENE_VRD2(F_CK, B_W, D, GAMMA_M_CONCRETO):
    """
    Esta função verifica o valor da resistência da biela comprimida V_RD2.

    Entrada:
    F_CK               | Resistência característica à compressão         | kN/m² | float
    B_W                | Largura da viga                                 | m     | float  
    D                  | Altura útil da seção                            | m     | float
    GAMMA_M_CONCRETO   | Coeficiente parcial de segurança do concreto    |       | float

    Saída:
    V_RD2              | Força resistente da biela comprimida            | kN    | float
    """
    V_RD2 = RESISTENCIA_BIELA_COMPRIMIDA(F_CK, B_W, D, GAMMA_M_CONCRETO)
    return V_RD

def RESISTENCIA_BIELA_COMPRIMIDA(F_CK, B_W, D, GAMMA_M_CONCRETO):
    """
    Esta função verifica o valor da resistência da biela comprimida V_RD2.

    Entrada:
    F_CK               | Resistência característica à compressão         | kN/m² | float
    B_W                | Largura da viga                                 | m     | float  
    D                  | Altura útil da seção                            | m     | float
    GAMMA_M_CONCRETO   | Coeficiente parcial de segurança do concreto    |       | float
    
    Saída:
    V_RD2              | Resitência da biela comprimida                  | kN    | float 
    """
    # Força resistente da biela de compressão
    F_CK /= 1E3 
    ALFA_V2 = (1 - (F_CK / 250))
    F_CK *= 1E3 
    F_CD = F_CK / GAMMA_M_CONCRETO
    V_RD2 = 0.27 * ALFA_V2 * F_CD * B_W * D
    return V_RD2

def FLECHA_LIMITE_DELTA_RD(L):
    """
    Esta função verifica o valor da flecha limite para uma viga.

    Entrada:
    L         | Vão da peça analisada           | m     | float
    
    Saída:
    DELTA_RD  | Flecha resistente               | m     | float 
    """
    # Flecha Limite
    DELTA_RD = (L / 250)
    return DELTA_RD

def VERIFICA_PILAR(TIPO_PILAR, A_C, RHO, E_S, F_CK, GAMMA_M_CONCRETO):
    """
    Esta função verifica o valor da carga normal máxima em um pilar de concreto armado.
    
    Entrada:
    TIPO_PILAR        | Tipo de pilar da estrutura                    |       | string
                      |    "INTERMEDIARIO"                            |       | 
                      |    "EXTREMIDADE"                              |       | 
                      |    "CANTO"                                    |       | 
    A_C               | Área da seção de concreto do pilar            | m²    | float
    RHO               | Taxa mecânica de armadura                     | %     | float
    E_S               | Módulo de elasticidade do aço                 | kN/m² | float
    F_CK              | Resistência característica à compressão       | kN/m² | float
    GAMMA_M_CONCRETO  | Coeficiente parcial de segurança do concreto  |       | float
    
    Saída:
    N_RD              | Compressão resistente do pilar                | m     | float 
    """
    # Definição do coeficiente de conversão de flexo-compressão para compressão
    if TIPO_PILAR == 'CANTO':
        ALPHA = 2.5
    elif TIPO_PILAR == 'EXTREMIDADE':
        ALPHA = 2.2
    elif TIPO_PILAR == 'INTERMEDIARIO':
        ALPHA = 1.8
    # Tensão de escoamento do aço a 0.20% de deformação de compressão
    SGIMA_SE = 0.2 / 100 * E_S
    # Carga máxima resistente
    F_CD = F_CK / GAMMA_M_CONCRETO
    RHO /= 100
    N_RD = (A_C * (0.85 * F_CD + RHO * SGIMA_SE)) / ALPHA
    return N_RD

import pandas as pd
import math
import numpy as np

def CALCULO_KV(TIPO_SOLO, NSPT, N_GK, N_QK):
    """
    Esta função determina o valor do K_V (Mola vertical) de uma estrutura. qwerqwerqwerwqe
    
    Entrada:
    TIPO_SOLO         | Tipo de solo encontrado na sondagem           |       | string
                      |    "INTERMEDIARIO"                            |       | 
                      |    "EXTREMIDADE"                              |       | 
                      |    "CANTO"                                    |       | 
    N_SPT             | Área da seção de concreto do pilar            | m²    | float
    N_GK              | Carga Normal Permanente característica do pilar  | kN    | float
    N_QK              | Carga Normal Variável característica do pilar    | kN    | float
    E_S               | Módulo de elasticidade do aço                 | kN/m² | float
    F_CK              | Resistência característica à compressão       | kN/m² | float
    GAMMA_M_CONCRETO  | Coeficiente parcial de segurança do concreto  |       | float
    
    Saída:
    N_RD              | Compressão resistente do pilar                | m     | float 
    """
    # Calculo deformabilidade do solo Es

    # Es = alfa * k * nspt 

  
    carga_pilar = SOLO["CARGA_PILAR"] 
    base_pilar = SOLO["BASE_PILAR"] 
    largura_pilar = SOLO["LARGURA_PILAR"]

    AUX = [('AREIA', '3', '900', '0.3'), ('AREIA ARGILO SILTOSA', '3', '600', '0.3'),
                ('AREIA ARGILOSA', '3', '550', '0.3'), ('AREIA SILTO ARGILOSA', '3', '575', '0.3'),
                ('AREIA SILTOSA', '3', '700', '0.3'), ('ARGILA', '7', '200', '0.45'),
                ('ARGILA ARENO SILTOSA', '6', '250', '0.45'), ('ARGILA ARENOSA', '6', '300', '0.45'),
                ('ARGILA SILTO sa', '6', '250', '0.45'), ('ARGILA SILTOSA', '6', '200', '0.45'),
                ('SILTE', '5', '350', '0.4'), ('SILTE  ARGILOSO', '5', '450', '0.4'),
                ('SILTE ARENOSO', '5', '450', '0.4'), ('SILTE ARGILOSO ARENOSO', '5', '325', '0.4'),
                ('SILTE ARGILOSO', '5', '250', '0.4')]          
    TABELA_SOLOS = pd.DataFrame(AUX)
    TABELA_SOLOS.columns = ['Solo', 'alfa', 'k(kPa)', 'n']
    DADO_SOLO = TABELA_SOLOS[(TABELA_SOLOS['Solo'] == TIPO_SOLO)]
    # Módulo de elasticidade do solo
    ALPHA = DADO_SOLO['alfa']
    K = DADO_SOLO['k(kPa)']
    N = DADO_SOLO['n']
    E_SOLO = float(ALPHA) * float(K) * float(N)
    # Tensão no solo
    B = 1
    N_SPTMEDIO = N_SPT / (3 * B) 
    SIGMA_SOLO = N_SPTMEDIO /  5     # Saída em kgf/cm²
    SIGMA_SOLO *= 100                # Conversão kgf/cm² para kN/m²
    # Área da sapata
    N_SK = 1.10 * N_GK + N_QK
    A_SAPATA = N_SK / SIGMA_SOLO
    L = A_SAPATA / B
    AUX = ((B ** 2) * L) ** (1 / 3) 
    # Coeficiente de mola
    K_VOLUME = 1.33 * E_SOLO / AUX
    K_V = A_SAPATA * K_VOLUME
    return K_V