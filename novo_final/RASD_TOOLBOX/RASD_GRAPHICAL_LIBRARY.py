#  /$$$$$$$   /$$$$$$   /$$$$$$  /$$$$$$$        /$$$$$$$$ /$$$$$$   /$$$$$$  /$$       /$$$$$$$   /$$$$$$  /$$   /$$
# | $$__  $$ /$$__  $$ /$$__  $$| $$__  $$      |__  $$__//$$__  $$ /$$__  $$| $$      | $$__  $$ /$$__  $$| $$  / $$
# | $$  \ $$| $$  \ $$| $$  \__/| $$  \ $$         | $$  | $$  \ $$| $$  \ $$| $$      | $$  \ $$| $$  \ $$|  $$/ $$/
# | $$$$$$$/| $$$$$$$$|  $$$$$$ | $$  | $$         | $$  | $$  | $$| $$  | $$| $$      | $$$$$$$ | $$  | $$ \  $$$$/ 
# | $$__  $$| $$__  $$ \____  $$| $$  | $$         | $$  | $$  | $$| $$  | $$| $$      | $$__  $$| $$  | $$  >$$  $$ 
# | $$  \ $$| $$  | $$ /$$  \ $$| $$  | $$         | $$  | $$  | $$| $$  | $$| $$      | $$  \ $$| $$  | $$ /$$/\  $$
# | $$  | $$| $$  | $$|  $$$$$$/| $$$$$$$/         | $$  |  $$$$$$/|  $$$$$$/| $$$$$$$$| $$$$$$$/|  $$$$$$/| $$  \ $$
# |__/  |__/|__/  |__/ \______/ |_______/          |__/   \______/  \______/ |________/|_______/  \______/ |__/  |__/

################################################################################
# UNIVERSIDADE FEDERAL DE CATALÃO (UFCAT)
# WANDERLEI MALAQUIAS PEREIRA JUNIOR                   ENG. CIVIL / PROF (UFCAT)
# ROMES ANTÔNIO BORGES                                       MAT. / PROF (UFCAT)
# DONIZETTI A. DE SOUZA JÚNIOR                                ENG. CIVIL (UFCAT)
################################################################################

################################################################################
# DESCRIÇÃO ALGORITMO:
# BIBLIOTECA GRÁFICA PARA ANÁLISE DE CONFIABILIDADE ESTRUTURAL DESENVOLVIDA 
# PELO GRUPO DE PESQUISAS E ESTUDOS EM ENGENHARIA (GPEE)
################################################################################

################################################################################
# BIBLIOTECAS NATIVAS PYTHON
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.cm import ScalarMappable
import numpy as np
import pandas as pd

################################################################################
# BIBLIOTECAS DESENVOLVEDORES GPEE 

# CONVERTE SI PARA POLEGADAS NO TAMANHO DO GRÁFICO
def CONVERT_SI_TO_INCHES(WIDTH, HEIGHT):
    """ 
    This function convert figure size meters to inches.
    
    Input:
    WIDTH    |  Figure width in SI units       | Float
    HEIGHT   |  Figure height in SI units      | Float
    
    Output:
    WIDTH    |  Figure width in INCHES units   | Float
    HEIGHT   |  Figure height in INCHES units  | Float
    """
    WIDTH = WIDTH / 0.0254
    HEIGHT = HEIGHT / 0.0254
    return WIDTH, HEIGHT

# SALVA A FIGURA
def SAVE_GRAPHIC(NAME, EXT, DPI):
    """ 
    This function save graphics on a specific path extensions options.

    - 'svg'
    - 'png'
    - 'eps'
    - 'pdf'

    Input: 
    NAME  | Path + name figure               | String
    EXT   | File extension                   | String
    DPI   | The resolution in dots per inch  | Integer
    
    Output:
    N/A
    """
    plt.savefig(NAME + EXT, dpi = DPI, bbox_inches = 'tight', transparent = True)


# PLOTAGEM 1
def RASD_PLOT_1(DATASET, PLOT_SETUP):
    """
    This function shows a boxplot and histograms in a single chart.
    
    Input: 
    DATASET     | Results from a RASD Toolboox                              | Py dataframe or Py Numpy array[N_POP x ??]
                |    Dictionary tags                                        |
                |    'DATASET'       == Complete data                       | Py dataframe or Py Numpy array[N_POP x ??]
                |    'COLUMN'        == Dataframe column                    | String
    PLOT_SETUP  | Contains specifications of each model of chart            | Py dictionary
                |    Dictionary tags                                        |
                |    'NAME'          == Filename output file                | String 
                |    'WIDTH'         == Width figure                        | Float
                |    'HEIGHT         == Height figure                       | Float
                |    'X AXIS SIZE'   == X axis size                         | Float
                |    'Y AXIS SIZE'   == Y axis size                         | Float
                |    'AXISES COLOR'  == Axis color                          | String
                |    'X AXIS LABEL'  == X label name                        | String
                |    'LABELS SIZE'   == Labels size                         | Float
                |    'LABELS COLOR'  == Labels color                        | String
                |    'CHART COLOR'   == Boxplot and histogram color         | String
                |    'BINS'          == Range representing the width of     | Float
                |                       a single bar                        | 
                |    'KDE'           == Smooth of the random distribution   | Boolean      
                |    'DPI'           == Dots Per Inch - Image quality       | Integer   
                |    'EXTENSION'     == Extension output file               | String ('.svg, '.png', '.eps' or '.pdf')
    
    Output:
    N/A

    Example:
    GRAPH_SETUP = { 'NAME': 'WANDER',
                    'WIDTH': 0.40, 
                    'HEIGHT': 0.20,
                    'DPI': 600, 
                    'BINS' : 10,
                    'EXTENSION': '.svg'
                    'X AXIS SIZE': 20,
                    'Y AXIS SIZE': 20,
                    'AXISES COLOR': '#000000',
                    'X AXIS LABEL': '$G(kN)$',
                    'LABELS SIZE': 16,
                    'LABELS COLOR': '#000000',  
                    'CHART COLOR': '#FEB625',
                    'KDE': False}
    DATA = RESULTS_TEST[0]['TOTAL RESULTS']
    DATA_SETUP = {'DATASET': DADOS, 'COLUMN': 'S_0'}
    RASD_PLOT_1(DATA_SETUP, GRAPH_SETUP)
    """
    # Setup
    NAME = PLOT_SETUP['NAME']
    W = PLOT_SETUP['WIDTH']
    H = PLOT_SETUP['HEIGHT']
    X_AXIS_LABEL = PLOT_SETUP['X AXIS LABEL']
    X_AXIS_SIZE = PLOT_SETUP['X AXIS SIZE']
    Y_AXIS_SIZE = PLOT_SETUP['Y AXIS SIZE']
    AXISES_COLOR = PLOT_SETUP['AXISES COLOR']
    LABELS_SIZE = PLOT_SETUP['LABELS SIZE']     
    LABELS_COLOR = PLOT_SETUP['LABELS COLOR']
    CHART_COLOR = PLOT_SETUP['CHART COLOR']
    BINS = int(PLOT_SETUP['BINS'])
    KDE = PLOT_SETUP['KDE']
    DPI = PLOT_SETUP['DPI']
    EXT = PLOT_SETUP['EXTENSION']
    AUX = DATASET['DATASET']
    COLUMN = DATASET['COLUMN']
    DATA = AUX[COLUMN]
    
    # Plot
    [W, H] = CONVERT_SI_TO_INCHES(W, H)
    sns.set(style = 'ticks')
    FIG, (AX_BOX, AX_HIST) = plt.subplots(2, figsize = (W, H), sharex = True, gridspec_kw = {'height_ratios': (.15, .85)})
    sns.boxplot(DATA, ax = AX_BOX, color = CHART_COLOR)
    sns.histplot(DATA, ax = AX_HIST, kde = KDE, color = CHART_COLOR, bins = BINS)
    AX_BOX.set(yticks = [])
    AX_BOX.set(xlabel='')
    font = {'fontname': 'Arial',
            'color':  LABELS_COLOR,
            'weight': 'normal',
            'size': LABELS_SIZE}
    AX_HIST.set_xlabel(X_AXIS_LABEL, fontdict = font)
    AX_HIST.set_ylabel('Frequência', fontdict = font)
    AX_HIST.tick_params(axis = 'x', labelsize = X_AXIS_SIZE, colors = AXISES_COLOR)
    AX_HIST.tick_params(axis = 'y', labelsize = Y_AXIS_SIZE, colors = AXISES_COLOR)
    sns.despine(ax = AX_HIST)
    sns.despine(ax = AX_BOX, left = True)
    
    # Save figure
    SAVE_GRAPHIC(NAME, EXT, DPI)

# PLOTAGEM 2
def RASD_PLOT_2(DATASET, PLOT_SETUP):
    """

    This function shows a scatter chart with results between limits and resistances demands.

    Input: 
    DATASET     | Results from a RASD Toolboox                             | Py dataframe or Py Numpy array[N_POP x ??]
                |    Dictionary tags                                       |
                |    'DATASET'       == Complete data                      | Py dataframe or Py Numpy array[N_POP x ??]
                |    'X DATA'        == Dataframe column plots in X axis   | String
                |    'Y DATA'        == Dataframe column plots in Y axis   | String
                |    'HUE VALUE'     == Dataframe column plots in Hue      | String
    PLOT_SETUP  | Contains specifications of each model of chart           | Py dictionary
                |    Dictionary tags                                       |
                |    'NAME'          == Filename output file               | String 
                |    'DPI'           == Dots Per Inch - Image quality      | Integer   
                |    'EXTENSION'     == Extension output file              | String ('.svg, '.png', '.eps' or '.pdf')
                |    'WIDTH'         == Width figure                       | Float
                |    'HEIGHT         == Height figure                      | Float
                |    'X AXIS SIZE'   == X axis size                        | Float
                |    'Y AXIS SIZE'   == Y axis size                        | Float
                |    'AXISES COLOR'  == Axis color                         | String
                |    'X AXIS LABEL'  == X label name                       | String
                |    'Y AXIS LABEL'  == Y label name                       | String             
                |    'LABELS SIZE'   == Labels size                        | Float
                |    'LABELS COLOR'  == labels color                       | Float
                |    'LOC LEGEND'    == Legend position                    | String
                |    'TITLE LEGEND'  == Text in legend                     | String
    
    Output:
    N/A

    Example:
    GRAPH_SETUP = { 'NAME': '5',
                    'EXTENSION': '.svg',
                    'DPI': 600,
                    'WIDTH': 0.300, 
                    'HEIGHT': 0.150,              
                    'X AXIS SIZE': 16,
                    'Y AXIS SIZE': 16,
                    'AXISES COLOR': '#000000',
                    'X AXIS LABEL': '$S_0 (N.m)$',
                    'Y AXIS LABEL': '$R_0 (N.m)$',
                    'LABELS SIZE': 18,
                    'LABELS COLOR': '#000000',
                    'LOC LEGEND': 'lower right',
                    'TITLE LEGEND': 'Índice de falha (I):'}
    DATA = RESULTS_TEST[0]['TOTAL RESULTS']
    DATA_SETUP = {'DATASET': DATA, 'X DATA': 'S_0', 'Y DATA': 'R_0', 'HUE VALUE': 'I_0'}    
    RASD_PLOT_2(DATA_SETUP, GRAPH_SETUP)
    """
    # Setup
    NAME = PLOT_SETUP['NAME']
    EXT = PLOT_SETUP['EXTENSION']
    DPI = PLOT_SETUP['DPI']
    W = PLOT_SETUP['WIDTH']
    H = PLOT_SETUP['HEIGHT']
    X_AXIS_SIZE = PLOT_SETUP['X AXIS SIZE']
    Y_AXIS_SIZE = PLOT_SETUP['Y AXIS SIZE']
    AXISES_COLOR = PLOT_SETUP['AXISES COLOR']
    X_AXIS_LABEL = PLOT_SETUP['X AXIS LABEL']
    Y_AXIS_LABEL = PLOT_SETUP['Y AXIS LABEL']
    LABELS_SIZE = PLOT_SETUP['LABELS SIZE']
    LABELS_COLOR = PLOT_SETUP['LABELS COLOR']
    LOC_LEGEND = PLOT_SETUP['LOC LEGEND']
    TITLE_LEGEND = PLOT_SETUP['TITLE LEGEND']
    DATA = DATASET['DATASET']
    X_DATA = DATASET['X DATA']
    Y_DATA = DATASET['Y DATA']
    HUE_VALUE = DATASET['HUE VALUE']
    
    # Plot
    sns.set(style = 'ticks')
    [W, H] = CONVERT_SI_TO_INCHES(W, H)
    FIG, AX = plt.subplots(figsize = (W, H))
    sns.scatterplot(data = DATA, x = X_DATA, y = Y_DATA, hue = HUE_VALUE)
    font = {'fontname': 'Arial',
        'color':  LABELS_COLOR,
        'weight': 'bold',
        'size': LABELS_SIZE}
    AX.set_xlabel(X_AXIS_LABEL, fontdict = font)
    AX.set_ylabel(Y_AXIS_LABEL, fontdict = font)
    AX.tick_params(axis = 'x', labelsize = X_AXIS_SIZE, colors = AXISES_COLOR)
    AX.tick_params(axis = 'y', labelsize = Y_AXIS_SIZE, colors = AXISES_COLOR)
    AX.legend(loc = LOC_LEGEND, title = TITLE_LEGEND)
    
    # Save figure
    SAVE_GRAPHIC(NAME, EXT, DPI)

# PLOTAGEM 3
def RASD_PLOT_3(DATASET, PLOT_SETUP):
    """
    This functions plots a scatter chart with variables X and Z in function of G value.
    
    Input: 
    DATASET     | Results from a RASD Toolboox                             | Py dataframe or Py Numpy array[N_POP x ??]
                |    Dictionary tags                                       |
                |    'DATASET'       == Complete data                      | Py dataframe or Py Numpy array[N_POP x ??]
                |    'X DATA'        == Dataframe name column plots in X   | String
                |    'Y DATA'        == Dataframe name column plots in Y   | String
                |    'G VALUE'       == Dataframe column plots in Hue      | String
    PLOT_SETUP  | Contains specifications of each model of chart           | Py dictionary
                |    Dictionary tags                                       |
                |    'NAME'          == Filename output file               | String 
                |    'EXTENSION'     == Extension output file              | String ('.svg, '.png', '.eps' or '.pdf')
                |    'DPI'           == Dots Per Inch - Image quality      | Integer   
                |    'WIDTH'         == Width figure                       | Float
                |    'HEIGHT         == Height figure                      | Float
                |    'X AXIS SIZE'   == X axis size                        | Float
                |    'Y AXIS SIZE'   == Y axis size                        | Float
                |    'AXISES COLOR'  == Axis color                         | String
                |    'X AXIS LABEL'  == X label name                       | String
                |    'Y AXIS LABEL'  == Y label name                       | String             
                |    'LABELS SIZE'   == Labels size                        | Float
                |    'LABELS COLOR'  == Labels color                       | Float
                |    'TRANSPARENCY'  == Blending value                     | Float
                |    'COLOR MAP'     == Colormap instance, Registered name | String

    Output:
    N/A

    Example:
    GRAPH_SETUP = { 'NAME': 'WANDER',
                    'EXTENSION': '.svg',
                    'DPI': 600,
                    'WIDTH': 0.20, 
                    'HEIGHT': 0.10,              
                    'X AXIS SIZE': 20,
                    'Y AXIS SIZE': 20,
                    'AXISES COLOR': '#000000',
                    'X AXIS LABEL': '$S_0$',
                    'Y AXIS LABEL': '$R_0$',
                    'LABELS SIZE': 16,
                    'LABELS COLOR': '#000000',
                    'TRANSPARENCY': 0.8,
                    'COLOR MAP': 'viridis'}
    DATA_SETUP = {'DATASET': DADOS, 'X DATA': 'S_0', 'Y DATA': 'R_0', 'G VALUE': 'G_0'}       
    RASD_PLOT_3(DATA_SETUP, GRAPH_SETUP)
    """
    # Setup
    NAME = PLOT_SETUP['NAME']
    EXT = PLOT_SETUP['EXTENSION']
    DPI = PLOT_SETUP['DPI']
    W = PLOT_SETUP['WIDTH']
    H = PLOT_SETUP['HEIGHT']
    X_AXIS_SIZE = PLOT_SETUP['X AXIS SIZE']
    Y_AXIS_SIZE = PLOT_SETUP['Y AXIS SIZE']
    AXISES_COLOR = PLOT_SETUP['AXISES COLOR']
    X_AXIS_LABEL = PLOT_SETUP['X AXIS LABEL']
    Y_AXIS_LABEL = PLOT_SETUP['Y AXIS LABEL']
    LABELS_SIZE = PLOT_SETUP['LABELS SIZE']
    LABELS_COLOR = PLOT_SETUP['LABELS COLOR']
    TRANSPARENCY = PLOT_SETUP['TRANSPARENCY']
    COLOR_MAP = PLOT_SETUP['COLOR MAP']
    A_UX = DATASET['DATASET']
    X_DATA = DATASET['X DATA']
    Y_DATA = DATASET['Y DATA']
    C_VALUE = DATASET['G VALUE']

    # Plot
    [W, H] = CONVERT_SI_TO_INCHES(W, H)
    AUX = plt.Normalize(A_UX[C_VALUE].min(), A_UX[C_VALUE].max())
    FIG, AX = plt.subplots(figsize = (W, H))
    plt.scatter(x = A_UX[X_DATA], y = A_UX[Y_DATA], c = A_UX[C_VALUE], cmap = COLOR_MAP, alpha = TRANSPARENCY)
    font = {'fontname': 'Arial',
        'color':  LABELS_COLOR,
        'weight': 'bold',
        'size': LABELS_SIZE}
    AX.set_xlabel(X_AXIS_LABEL, fontdict = font)
    AX.set_ylabel(Y_AXIS_LABEL, fontdict = font)
    AX.tick_params(axis = 'x', labelsize = X_AXIS_SIZE, colors = AXISES_COLOR)
    AX.tick_params(axis = 'y', labelsize = Y_AXIS_SIZE, colors = AXISES_COLOR)
    AUX1 =  ScalarMappable(norm = AUX, cmap = COLOR_MAP)
    FIG.colorbar(AUX1, ax = AX)
    
    # Save figure
    SAVE_GRAPHIC(NAME, EXT, DPI)

# PLOTAGEM 4
def RASD_PLOT_4(DATASET, PLOT_SETUP):
    """
    This function plots two histograms in a single one chart.

    Input: 
    DATASET     | Results from a RASD Toolboox                             | Py dataframe or Py Numpy array[N_POP x ??]
                |    Dictionary tags                                       |
                |    'DATASET'       == Complete data                      | Py dataframe or Py Numpy array[N_POP x ??]
                |    'HIST_1'        == Dataframe hist 1                   | String
                |    'HIST_2'        == Dataframe hist 2                   | String
    PLOT_SETUP  | Contains specifications of each model of chart           | Py dictionary
                |    Dictionary tags                                       |
                |    'NAME'          == Filename output file               | String 
                |    'EXTENSION'     == Extension output file              | String ('.svg, '.png', '.eps' or '.pdf')
                |    'DPI'           == Dots Per Inch - Image quality      | Integer   
                |    'WIDTH'         == Width figure                       | Float
                |    'HEIGHT         == Height figure                      | Float
                |    'X AXIS SIZE'   == X axis size                        | Float
                |    'Y AXIS SIZE'   == Y axis size                        | Float
                |    'HIST 1 LABEL'  == Histogram 1 label                  | String
                |    'HIST 2 LABEL'  == Histogram 2 label                  | String
                |    'HIST 1 COLOR'  == Histogram 1 color                  | String
                |    'HIST 2 COLOR'  == Histogram 2 color                  | String
                |    'X AXIS LABEL'  == X label name                       | String
                |    'Y AXIS LABEL'  == Y label name                       | String             
                |    'LABELS SIZE'   == Labels size                        | Float
                |    'LABELS COLOR'  == Labels color                       | String
                |    'TRANSPARENCY'  == Blending value                     | Float
                |    'BINS'          == Equal width bins in the range      | Integer

    Output:
    N/A

    Example:
    GRAPH_SETUP = { 'NAME': 'WANDER',
                        'EXTENSION': '.svg',
                        'DPI': 600,
                        'WIDTH': 0.20,
                        'HEIGHT': 0.10,
                        'X AXIS SIZE': 20,
                        'Y AXIS SIZE': 20,
                        'AXISES COLOR': '#00000',
                        'HIST 1 LABEL': '$R_0$',
                        'HIST 2 LABEL': '$S_0$',
                        'HIST 1 COLOR': '#0000FF',
                        'HIST 2 COLOR': '#00FF7F',
                        'X AXIS LABEL': '$Z_0$',
                        'Y AXIS LABEL': 'Frequência',
                        'LABELS SIZE': 16,
                        'LABELS COLOR': '#000000',
                        'BINS': 10,
                        'TRANSPARENCY': 0.80}
    DATA_SETUP = {'DATASET': DADOS, 'HIST_1': 'R_0', 'HIST_2': 'S_0'}  
    RASD_PLOT_4(DATA_SETUP, GRAPH_SETUP)
    """
    # Setup
    NAME = PLOT_SETUP['NAME']
    EXT = PLOT_SETUP['EXTENSION']
    DPI = PLOT_SETUP['DPI']
    W = PLOT_SETUP['WIDTH']
    H = PLOT_SETUP['HEIGHT']
    X_AXIS_SIZE = PLOT_SETUP['X AXIS SIZE']
    Y_AXIS_SIZE = PLOT_SETUP['Y AXIS SIZE']
    AXISES_COLOR = PLOT_SETUP['AXISES COLOR']
    X_AXIS_LABEL = PLOT_SETUP['X AXIS LABEL']
    Y_AXIS_LABEL = PLOT_SETUP['Y AXIS LABEL']
    HIST_1_LABEL = PLOT_SETUP['HIST 1 LABEL']
    HIST_2_LABEL = PLOT_SETUP['HIST 2 LABEL']
    HIST_1_COLOR = PLOT_SETUP['HIST 1 COLOR']
    HIST_2_COLOR = PLOT_SETUP['HIST 2 COLOR']
    LABELS_SIZE = PLOT_SETUP['LABELS SIZE']
    LABELS_COLOR = PLOT_SETUP['LABELS COLOR']
    BINS = int(PLOT_SETUP['BINS'])
    ALPHA = float(PLOT_SETUP['TRANSPARENCY'])
    A_UX = DATASET['DATASET']
    X_DATA = DATASET['HIST_1']
    Y_DATA = DATASET['HIST_2']

    # Plot
    [W, H] = CONVERT_SI_TO_INCHES(W, H)
    plt.subplots(figsize = (W, H))
    plt.hist(A_UX[X_DATA], bins = BINS, label = HIST_1_LABEL, facecolor = HIST_1_COLOR, alpha = ALPHA)
    plt.hist(A_UX[Y_DATA], bins = BINS, label = HIST_2_LABEL, facecolor = HIST_2_COLOR, alpha = ALPHA)
    plt.legend()
    plt.xlabel(X_AXIS_LABEL)
    plt.ylabel(Y_AXIS_LABEL)

    # Save figure
    SAVE_GRAPHIC(NAME, EXT, DPI)

    # PLOTAGEM 5
def RASD_PLOT_5(DATASET, PLOT_SETUP):
    """
    This function plots a chart with values of number of simulations versus Beta/Failure Probability.

    Input: 
    DATASET     | Results from a RASD Toolboox                             | Py dataframe or Py Numpy array[N_POP x ??]
                |    Dictionary tags                                       |
                |    'DATASET'       == Complete data                      | Py dataframe or Py Numpy array[N_POP x ??]
                |    'G NUMBER'      == State Limit function plot          | Int
    PLOT_SETUP  | Contains specifications of each model of chart           | Py dictionary
                |    Dictionary tags                                       |
                |    'NAME'          == Filename output file               | String 
                |    'EXTENSION'     == Extension output file              | String ('.svg, '.png', '.eps' or '.pdf')
                |    'DPI'           == Dots Per Inch - Image quality      | Integer 
                |    'WIDTH'         == Width figure                       | Float
                |    'HEIGHT         == Height figure                      | Float
                |    'X AXIS SIZE'   == X axis size                        | Float
                |    'Y AXIS SIZE'   == Y axis size                        | Float
                |    'AXISES COLOR'  == Axis color                         | String
                |    'X AXIS LABEL'  == X label name                       | String
                |    'Y AXIS LABEL'  == Y label name                       | String             
                |    'LABELS SIZE'   == Labels size                        | Float
                |    'LABELS COLOR'  == Labels color                       | Float
                |    'CHART COLOR'   == Chart color                        | String
                |    'POPULATION'    == Population array                   | List
                |    'TYPE'          == Chart type                         | String ('Beta', 'Pf', 'N_fails')

    Output:
    N/A

    Example:

    """
    # Setup 
    NAME = PLOT_SETUP['NAME']
    EXT = PLOT_SETUP['EXTENSION']
    DPI = PLOT_SETUP['DPI']
    W = PLOT_SETUP['WIDTH']
    H = PLOT_SETUP['HEIGHT']
    X_AXIS_SIZE = PLOT_SETUP['X AXIS SIZE']
    Y_AXIS_SIZE = PLOT_SETUP['Y AXIS SIZE']
    AXISES_COLOR = PLOT_SETUP['AXISES COLOR']
    X_AXIS_LABEL = PLOT_SETUP['X AXIS LABEL']
    Y_AXIS_LABEL = PLOT_SETUP['Y AXIS LABEL']
    LABELS_SIZE = PLOT_SETUP['LABELS SIZE']
    LABELS_COLOR = PLOT_SETUP['LABELS COLOR']
    CHART_COLOR = PLOT_SETUP['CHART COLOR']
    POP_SIZE = len(PLOT_SETUP['POPULATION'])
    POPULATION = PLOT_SETUP['POPULATION']
    CHART_TYPE = PLOT_SETUP['TYPE']
    DATA = DATASET['DATASET']
    G_NUMBER = DATASET['G NUMBER']
    PF_AUX = []
    BETA_AUX = []  
    N_FAILS = [] 
    # Plot
    [W, H] = CONVERT_SI_TO_INCHES(W, H)
    for I_COUNT in range(POP_SIZE):
        PF_VALUE = DATA[I_COUNT]['PROBABILITY OF FAILURE'][G_NUMBER]
        BETA_VALUE = DATA[I_COUNT]['BETA INDEX'][G_NUMBER]
        N_FAILSVALUE = DATA[I_COUNT]['NUMBER OF FAILURES'][G_NUMBER]
        PF_AUX.append(PF_VALUE)
        BETA_AUX.append(BETA_VALUE)
        N_FAILS.append(N_FAILSVALUE)
    plt.subplots(figsize = (W, H))    
    if CHART_TYPE.upper() == 'PF':
        plt.plot(POPULATION, PF_AUX, color = CHART_COLOR)
    elif CHART_TYPE.upper() == 'BETA':
        plt.plot(POPULATION, BETA_AUX, color = CHART_COLOR)
    elif CHART_TYPE.upper() == 'N_FAILS':
        plt.plot(POPULATION, N_FAILS, color = CHART_COLOR)
    plt.xlabel(X_AXIS_LABEL)
    plt.ylabel(Y_AXIS_LABEL)
    
    # Save figure
    SAVE_GRAPHIC(NAME, EXT, DPI)

#   /$$$$$$  /$$$$$$$  /$$$$$$$$ /$$$$$$$$       /$$$$$$$$ /$$$$$$$$  /$$$$$$  /$$   /$$ /$$   /$$  /$$$$$$  /$$        /$$$$$$   /$$$$$$  /$$$$$$ /$$$$$$$$  /$$$$$$ 
#  /$$__  $$| $$__  $$| $$_____/| $$_____/      |__  $$__/| $$_____/ /$$__  $$| $$  | $$| $$$ | $$ /$$__  $$| $$       /$$__  $$ /$$__  $$|_  $$_/| $$_____/ /$$__  $$
# | $$  \__/| $$  \ $$| $$      | $$               | $$   | $$      | $$  \__/| $$  | $$| $$$$| $$| $$  \ $$| $$      | $$  \ $$| $$  \__/  | $$  | $$      | $$  \__/
# | $$ /$$$$| $$$$$$$/| $$$$$   | $$$$$            | $$   | $$$$$   | $$      | $$$$$$$$| $$ $$ $$| $$  | $$| $$      | $$  | $$| $$ /$$$$  | $$  | $$$$$   |  $$$$$$ 
# | $$|_  $$| $$____/ | $$__/   | $$__/            | $$   | $$__/   | $$      | $$__  $$| $$  $$$$| $$  | $$| $$      | $$  | $$| $$|_  $$  | $$  | $$__/    \____  $$
# | $$  \ $$| $$      | $$      | $$               | $$   | $$      | $$    $$| $$  | $$| $$\  $$$| $$  | $$| $$      | $$  | $$| $$  \ $$  | $$  | $$       /$$  \ $$
# |  $$$$$$/| $$      | $$$$$$$$| $$$$$$$$         | $$   | $$$$$$$$|  $$$$$$/| $$  | $$| $$ \  $$|  $$$$$$/| $$$$$$$$|  $$$$$$/|  $$$$$$/ /$$$$$$| $$$$$$$$|  $$$$$$/
#  \______/ |__/      |________/|________/         |__/   |________/ \______/ |__/  |__/|__/  \__/ \______/ |________/ \______/  \______/ |______/|________/ \______/ 