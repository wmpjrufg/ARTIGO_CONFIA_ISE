################################################################################
# UNIVERSIDADE FEDERAL DE CATALÃO (UFCAT)
# WANDERLEI MALAQUIAS PEREIRA JUNIOR,                 ENG. CIVIL / PROF. (UFCAT)
# ANDRÉ TEÓFILO BECK,                                   ENG. MEC. / PROF. (EESC)
# DANIEL LIMA ARAÚJO,                                   ENG. CIVIL / PROF. (UFG)
# MAURO,                                               ENG. CIVIL / PROF. (UFPA)
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

def MOMENTO_RESISTENE_MRD(B_W, D, A_S, F_Y, F_CK, GAMMA_M_ACO, GAMMA_M_CONCRETO, THETA_R):
    """
    Esta função determina o momento resistente M_RD de uma viga de concreto armado para análise de confiabilidade.

    Entrada:
    B_W                | Largura da viga                               | m     | float  
    D                  | Altura útil da seção                          | m     | float
    A_S                | Área de aço necessária na seção               | m²    | float
    F_Y                | Tensão de escoamento do aço                   | kN/m² | float
    GAMMA_M_ACO        | Coeficiente parcial de segurança do aço       |       | float
    GAMMA_M_CONCRETO   | Coeficiente parcial de segurança do concreto  |       | float
    THETA_R            | Erro de modelo da capacidade                  |       | float

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
    # Linha neutra real e momento resistente
    X_III = (A_S * (F_Y / GAMMA_M_ACO)) / (ALPHA_C * (F_CK * GAMMA_M_CONCRETO) * B_W * LAMBDA)
    M_RD = A_S * (F_Y / GAMMA_M_ACO) * (D - 0.50 * LAMBDA * X_III)
    # Capacidade
    M_RD *= THETA_R
    return M_RD

def CORTANTE_RESISTENE_VRD2(F_CK, B_W, D, GAMMA_M_CONCRETO, THETA_R):
    """
    Esta função verifica o valor da resistência da biela comprimida V_RD2.

    Entrada:
    F_CK               | Resistência característica à compressão         | kN/m² | float
    B_W                | Largura da viga                                 | m     | float  
    D                  | Altura útil da seção                            | m     | float
    GAMMA_M_CONCRETO   | Coeficiente parcial de segurança do concreto    |       | float
    THETA_R            | Erro de modelo da capacidade                    |       | float
    """
    # Capacidade
    V_RD2 = RESISTENCIA_BIELA_COMPRIMIDA(F_CK, B_W, D, GAMMA_M_CONCRETO)
    V_RD = V_RD2 * THETA_R
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

def FLECHA_LIMITE_DELTA_RD(L, THETA_R):
    """
    Esta função verifica o valor da flecha limite para uma análise de confiabilidade.

    Entrada:
    L         | Vão da peça analisada           | m     | float
    THETA_R   | Erro de modelo da capacidade    |       | float
    
    Saída:
    DELTA_RD  | Flecha resistente               | m     | float 
    """
    # Flecha Limite
    DELTA_RD = (L / 250) * THETA_R 
    return DELTA_RD