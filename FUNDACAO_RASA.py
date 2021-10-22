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
# BIBLIOTECA QUE DETERMINA A CONSTANTE DE MOLA E ÁREA DE UMA SAPATA DESENVOL-
# VIDA PELO GRUPO DE PESQUISAS E ESTUDOS EM ENGENHARIA (GPEE)
################################################################################

################################################################################
# BIBLIOTECAS DESENVOLVEDORES GPEE

def AREA_SAPATA(N_GK, N_QK, N_SPT, COTA_SAPATA, A_P, B_P):
	"""
	Esta função determina a área de uma sapata em função do carregamento e da tensão admíssivel do solo.

	Entrada:
	N_GK         | Carga normal permanente característica do pilar  | kN    | float
    N_QK         | Carga normal variável característica do pilar    | kN    | float
	N_SPT        | Valor do ensaio de sondagem tipo N_SPT           |       | list
	COTA_SAPATA  | Cota de assentamento da sapata                   | m     | int
	A_P          | Maior dimensão do pilar                          | m     | float
	B_P          | Menor dimensão do pilar                          | m     | float

	Saída:
	A_SAPATA     | Área da sapata                                   | m²    | float
    """
   	# Carga de projeto
    N_SK = 1.10 * N_GK + N_QK
    # Área da sapata
    B = 1
    SIGMA_SOLOADM = TENSAO_ADMISSIVEL(N_SPT, COTA_SAPATA, B)
    A = (N_SK / SIGMA_SOLO) / B
	A_SAPATA = A * B
    return A_SAPATA

def TENSAO_ADMISSIVEL(N_SPT, COTA_SAPATA, B):
	"""
	Esta função determina a tensão admissível no solo a partir de dados de uma sondagem.

	Entrada:
	N_SPT        | Valor do ensaio de sondagem tipo N_SPT   |       | list
	COTA_SAPATA  | Cota de assentamento da sapata           | m     | int
	B            | Menor dimensão da sapata                 | m     | int

	Saída:
	SIGMA_ADM    | Valor da tensão admissível média do solo | kN/m² | float 
	"""
	COTA_BULBO = 2 * B
	COTA_FINALCOMBULBO = COTA_SAPATA + round(COTA_BULBO)
	PERFIL_BULBO = []
	PERFIL_INICIAL = []
	for I_CONT in range(COTA_FINALCOMBULBO):
	    PERFIL_BULBO.append(I_CONT) 
	for J_CONT in range(COTA_SAPATA):
	    PERFIL_INICIAL.append(J_CONT ) 
	PROFUNDIDADES = np.setdiff1d(PERFIL_BULBO, PERFIL_INICIAL)
	N_SPTBULBO = []
	for K_CONT, VALOR_NSPT in enumerate(PROFUNDIDADES):
	    N_SPTBULBO.append(N_SPT[VALOR_NSPT])
	SIGMA_ADM = np.mean(N_SPTBULBO) / 0.05
	return SIGMA_ADM

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

    L = A_SAPATA / B
    AUX = ((B ** 2) * L) ** (1 / 3) 
    # Coeficiente de mola
    K_VOLUME = 1.33 * E_SOLO / AUX
    K_V = A_SAPATA * K_VOLUME
    return K_V