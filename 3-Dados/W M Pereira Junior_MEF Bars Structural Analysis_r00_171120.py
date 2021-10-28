#!/usr/bin/env python
# coding: utf-8

# In[1]:


################################################################################
# FEDERAL UNIVERSITY OF CATALÃO (UFCAT)
# DEVELOPERS:
# JOSÉ VITOR CARVALHO SILVA, CIVIL ENGINEER (UFCAT / FENG)
# PROF. WANDERLEI MALAQUIAS PEREIRA JUNIOR, Ph.D. (UFCAT / FENG)
# PROF. DAVIDSON FRANÇA JUNIOR, MSc. (UFCAT / FENG)
# GITHUB: wwww.linklink.com.br
# LICENSE: aqui verificar
################################################################################

################################################################################
# DESCRIPTION SCRIPT:

# LIBRARY FOR DETERMINING INERNAL LOADS IN BARS STRUCTURES BY FINITE ELEMENT 
# METHOD 
################################################################################

################################################################################
# LIBRARIES
# NATIVE PYTHON
import numpy as np

# GPEE DEVELOPERS
def INDEX_ASSEMBLY(TYPE_ELEMENT):
    if TYPE_ELEMENT == "BAR_2D_1DOF":  
        TOTAL_DOF_NODE = 1
        TOTAL_NODES_ELEMENT = 2
        ARRAY_DOF_ACTIVE = np.array([[1, 0, 0]])
        ARRAY_LOCAL_DOF = [0, 1]
    elif TYPE_ELEMENT == "TRUSS_2D_4DOF":   
        TOTAL_DOF_NODE = 2
        TOTAL_NODES_ELEMENT = 2
        ARRAY_DOF_ACTIVE = np.array([[1, 1, 0]])
        ARRAY_LOCAL_DOF = [0, 1, 2, 3]
    elif TYPE_ELEMENT == "BEAM_2D_4DOF":   
        TOTAL_DOF_NODE = 2
        TOTAL_NODES_ELEMENT = 2
        ARRAY_DOF_ACTIVE = np.array([[0, 1, 1]])
        ARRAY_LOCAL_DOF = [0, 1, 2, 3]
    elif TYPE_ELEMENT == "FRAME_2D_6DOF": 
        TOTAL_DOF_NODE = 3
        TOTAL_NODES_ELEMENT = 2
        ARRAY_DOF_ACTIVE = np.array([[1, 1, 1]])
        ARRAY_LOCAL_DOF = [0, 1, 2, 3, 4, 5]
    COLUMN_MATERIAL_INDEX = TOTAL_NODES_ELEMENT + 0
    COLUMN_SECTION_INDEX = TOTAL_NODES_ELEMENT + 1
    TOTAL_DOF_ELEMENT = TOTAL_DOF_NODE * TOTAL_NODES_ELEMENT
    return TOTAL_DOF_NODE, TOTAL_NODES_ELEMENT, TOTAL_DOF_ELEMENT, COLUMN_MATERIAL_INDEX, COLUMN_SECTION_INDEX, ARRAY_DOF_ACTIVE, ARRAY_LOCAL_DOF

def GLOBAL_DOF_ASSEMBLY(TYPE_ELEMENT, TOTAL_DOF_NODE, TOTAL_NODES, TOTAL_ELEMENTS):
    MATRIX_GLOBAL_DOF = np.zeros((TOTAL_NODES,3))
    for I_COUNT in range (TOTAL_NODES):
        if TYPE_ELEMENT == "BAR_2D_1DOF":  
            MATRIX_GLOBAL_DOF[I_COUNT, 0] = int(TOTAL_DOF_NODE * I_COUNT + 0) 
            MATRIX_GLOBAL_DOF[I_COUNT, 1] = - 1989
            MATRIX_GLOBAL_DOF[I_COUNT, 2] = - 1989
        elif TYPE_ELEMENT == "TRUSS_2D_4DOF":   
            MATRIX_GLOBAL_DOF[I_COUNT, 0] = int(TOTAL_DOF_NODE * I_COUNT + 0)
            MATRIX_GLOBAL_DOF[I_COUNT, 1] = int(TOTAL_DOF_NODE * I_COUNT + 1)
            MATRIX_GLOBAL_DOF[I_COUNT, 2] = - 1989
        elif TYPE_ELEMENT == "BEAM_2D_4DOF":   
            MATRIX_GLOBAL_DOF[I_COUNT, 0] = - 1989
            MATRIX_GLOBAL_DOF[I_COUNT, 1] = int(TOTAL_DOF_NODE * I_COUNT + 0)
            MATRIX_GLOBAL_DOF[I_COUNT, 2] = int(TOTAL_DOF_NODE * I_COUNT + 1)
        elif TYPE_ELEMENT == "FRAME_2D_6DOF": 
            MATRIX_GLOBAL_DOF[I_COUNT, 0] = int(TOTAL_DOF_NODE * I_COUNT + 0)
            MATRIX_GLOBAL_DOF[I_COUNT, 1] = int(TOTAL_DOF_NODE * I_COUNT + 1)
            MATRIX_GLOBAL_DOF[I_COUNT, 2] = int(TOTAL_DOF_NODE * I_COUNT + 2)
    return MATRIX_GLOBAL_DOF

def TOTAL_DEGREE_FREEDOM(TOTAL_DOF_NODE, TOTAL_NODES):
    ARRAY_TOTAL_DOF = []
    TOTAL_DOF = TOTAL_NODES * TOTAL_DOF_NODE
    for I_COUNT in range (TOTAL_DOF):
        ARRAY_TOTAL_DOF.append (I_COUNT)
    return ARRAY_TOTAL_DOF, TOTAL_DOF

def PRESCRIBED_DEGREE_FREEDOM(MATRIX_NODAL_PRESCRIPTIONS, MATRIX_GLOBAL_DOF, TOTAL_DOF_NODE):
    ARRAY_PRESCRIBED_DOF = []
    ARRAY_PRESCRIBED_DOF_VALUE = []
    TOTAL_PRESCRIBED_DOF = MATRIX_NODAL_PRESCRIPTIONS.shape[0]
    for I_COUNT in range(TOTAL_PRESCRIBED_DOF):
        NODE = MATRIX_NODAL_PRESCRIPTIONS[I_COUNT, 0]
        INDEX_DOF = int(MATRIX_NODAL_PRESCRIPTIONS[I_COUNT, 1])
        DOF_VALUE = int(MATRIX_GLOBAL_DOF[NODE, INDEX_DOF])
        ARRAY_PRESCRIBED_DOF.append(DOF_VALUE)
        PRESCRIBED_VALUE = MATRIX_NODAL_PRESCRIPTIONS[I_COUNT, 2]
        ARRAY_PRESCRIBED_DOF_VALUE.append(PRESCRIBED_VALUE)
    return TOTAL_PRESCRIBED_DOF, ARRAY_PRESCRIBED_DOF_VALUE, ARRAY_PRESCRIBED_DOF

def FREE_DEGREE_FREEDOM(ARRAY_PRESCRIBED_DOF, ARRAY_TOTAL_DOF):
    ARRAY_FREE_DOF = np.setdiff1d(ARRAY_TOTAL_DOF, ARRAY_PRESCRIBED_DOF)
    TOTAL_FREE_DOF = len(ARRAY_FREE_DOF)
    return TOTAL_FREE_DOF, ARRAY_FREE_DOF

def NODAL_EXTERNAL_LOAD(MATRIX_NODAL_EXTERNAL_LOAD, TOTAL_NODAL_LOADS, TOTAL_DOF, MATRIX_GLOBAL_DOF):
    ARRAY_COMPLETE_NODAL_FORCE = np.zeros((TOTAL_DOF,1))
    for I_COUNT in range(TOTAL_NODAL_LOADS):
        NODE = int(MATRIX_NODAL_EXTERNAL_LOAD[I_COUNT, 0])
        INDEX_DOF = int(MATRIX_NODAL_EXTERNAL_LOAD[I_COUNT, 1])
        DOF_VALUE = int(MATRIX_GLOBAL_DOF[NODE, INDEX_DOF])
        LOAD = int(MATRIX_NODAL_EXTERNAL_LOAD[I_COUNT, 2])
        ARRAY_COMPLETE_NODAL_FORCE[DOF_VALUE, 0] = LOAD
    return ARRAY_COMPLETE_NODAL_FORCE

def DISTRIBUTED_EXTERNAL_LOAD(TOTAL_DOF, ARRAY_ELEMENT_GEOMETRY, TYPE_ELEMENT, TYPE_EXTERNAL_LOAD, MATRIX_DISTRIBUTED_EXTERNAL_LOAD, TOTAL_NODES_ELEMENT, MATRIX_ELEMENT_PROPERTIES, MATRIX_GLOBAL_DOF):
    ARRAY_COMPLETE_EQUIVALENT_NODAL_FORCE = np.zeros((TOTAL_DOF,1))
    L = ARRAY_ELEMENT_GEOMETRY[0]
    SINN = ARRAY_ELEMENT_GEOMETRY[1]
    COSS = ARRAY_ELEMENT_GEOMETRY[2]
    if TYPE_ELEMENT == "BEAM_2D_4DOF":
        if TYPE_EXTERNAL_LOAD == 0:
            QY = MATRIX_DISTRIBUTED_EXTERNAL_LOAD[I_COUNT, 3]
            ARRAY_EQUIVALENT_LOCAL_NODAL_FORCE = np.array([[QY * L / 2], [(QY * L * L) / 12], [QY * L / 2], [- (QY * L * L) / 12]])
            MATRIX_TRANSFORM = ROTATION(TYPE_ELEMENT, ARRAY_ELEMENT_GEOMETRY)
            ARRAY_EQUIVALENT_GLOBAL_NODAL_FORCE = np.dot(np.transpose(MATRIX_TRANSFORM), ARRAY_EQUIVALENT_LOCAL_NODAL_FORCE)
    elif TYPE_ELEMENT == "FRAME_2D_6DOF":
        if TYPE_EXTERNAL_LOAD == 0:
            QX = MATRIX_DISTRIBUTED_EXTERNAL_LOAD[I_COUNT, 2]
            QY = MATRIX_DISTRIBUTED_EXTERNAL_LOAD[I_COUNT, 3]
            ARRAY_EQUIVALENT_LOCAL_NODAL_FORCE = np.array([[QX * L / 2], [QY * L / 2], [(QY * L * L) / 12], [QX * L / 2], [QY * L / 2], [- (QY * L * L) / 12]])
            MATRIX_TRANSFORM = ROTATION(TYPE_ELEMENT, ARRAY_ELEMENT_GEOMETRY)
            ARRAY_EQUIVALENT_GLOBAL_NODAL_FORCE = np.dot(np.transpose(MATRIX_TRANSFORM), ARRAY_EQUIVALENT_LOCAL_NODAL_FORCE)
    for I_COUNT in range(TOTAL_NODES_ELEMENT):
        NODE = MATRIX_ELEMENT_PROPERTIES[I_ELEMENT, I_COUNT]
        COUNT = 0
        for J_COUNT in range(TOTAL_DOF_NODE):
            DOF_VALUE = int(MATRIX_GLOBAL_DOF[NODE, J_COUNT])
            ARRAY_COMPLETE_EQUIVALENT_NODAL_FORCE[DOF_VALUE, 0] = ARRAY_EQUIVALENT_LOCAL_NODAL_FORCE[COUNT, 0]
            COUNT = COUNT + 1
    return ARRAY_COMPLETE_EQUIVALENT_NODAL_FORCE

def MATERIALS_PROPRETIES(MATRIX_ELEMENT_PROPERTIES, MATRIX_MATERIAL_PROPERTIES, I_ELEMENT, INDEX_COLUMN_MAT):
    MATERIAL_INDEX = MATRIX_ELEMENT_PROPERTIES[I_ELEMENT, INDEX_COLUMN_MAT]
    E_MODULUS_VALUE = MATRIX_MATERIAL_PROPERTIES[MATERIAL_INDEX, 0]
    POISSON_VALUE = MATRIX_MATERIAL_PROPERTIES[MATERIAL_INDEX, 1]
    THERMAL_COEFFICIENT = MATRIX_MATERIAL_PROPERTIES[MATERIAL_INDEX, 2]
    G_MODULUS_VALUE = E_MODULUS_VALUE / (2 * (1 + POISSON_VALUE))
    ARRAY_ELEMENT_MATERIAL = [E_MODULUS_VALUE, G_MODULUS_VALUE, POISSON_VALUE, THERMAL_COEFFICIENT]
    return ARRAY_ELEMENT_MATERIAL

def GEOMETRIC_PROPRETIES(MATRIX_COORDINATES_PROPERTIES, MATRIX_ELEMENT_PROPERTIES, SECTION_MATRIX_PROPERTIES, I_ELEMENT, INDEX_COLUMN_SEC):
    NODE_1 = MATRIX_ELEMENT_PROPERTIES[I_ELEMENT, 0]
    NODE_2 = MATRIX_ELEMENT_PROPERTIES[I_ELEMENT, 1]
    X_NODE_1 = MATRIX_COORDINATES_PROPERTIES[NODE_1, 0]
    Y_NODE_1 = MATRIX_COORDINATES_PROPERTIES[NODE_1, 1]
    X_NODE_2 = MATRIX_COORDINATES_PROPERTIES[NODE_2, 0]
    Y_NODE_2 = MATRIX_COORDINATES_PROPERTIES[NODE_2, 1]
    DELTA_X = X_NODE_2 - X_NODE_1
    DELTA_Y = Y_NODE_2 - Y_NODE_1
    LENGTH = ((DELTA_X) ** 2 + (DELTA_Y) ** 2) ** 0.50
    COSS = DELTA_X / LENGTH
    SINN = DELTA_Y / LENGTH
    SECTION_INDEX = MATRIX_ELEMENT_PROPERTIES[I_ELEMENT, INDEX_COLUMN_SEC]
    AREA_VALUE = SECTION_MATRIX_PROPERTIES[SECTION_INDEX, 0]
    INERTIAX_VALUE = SECTION_MATRIX_PROPERTIES[SECTION_INDEX, 1]
    INERTIAY_VALUE = SECTION_MATRIX_PROPERTIES[SECTION_INDEX, 2]
    ARRAY_ELEMENT_GEOMETRY = [LENGTH, SINN, COSS, AREA_VALUE, INERTIAX_VALUE, INERTIAY_VALUE]
    return ARRAY_ELEMENT_GEOMETRY

def ELEMENTAR_STIFFNESS(TYPE_ELEMENT, ARRAY_ELEMENT_GEOMETRY, ARRAY_ELEMENT_MATERIAL):
    if TYPE_ELEMENT == "BAR_2D_1DOF":
        L = ARRAY_ELEMENT_GEOMETRY[0]
        A = ARRAY_ELEMENT_GEOMETRY[3]
        E = ARRAY_ELEMENT_MATERIAL[0]
        MATRIX_ELEMENTAR_STIFFNESS = (A * E / L) * np.array([[1, -1],
                                                            [-1, 1]])        
    elif TYPE_ELEMENT == "TRUSS_2D_4DOF":
            L = ARRAY_ELEMENT_GEOMETRY[0]
            A = ARRAY_ELEMENT_GEOMETRY[3]
            E = ARRAY_ELEMENT_MATERIAL[0]
            MATRIX_ELEMENTAR_STIFFNESS = (A * E / L) * np.array([[1, 0,-1, 0],
                                                                [0, 0, 0, 0],
                                                                [-1, 0, 1, 0],
                                                                [0, 0, 0, 0]])
    elif TYPE_ELEMENT == "BEAM_2D_4DOF":
            L = ARRAY_ELEMENT_GEOMETRY[0]
            I = ARRAY_ELEMENT_GEOMETRY[5]
            E = ARRAY_ELEMENT_MATERIAL[0]
            MATRIX_ELEMENTAR_STIFFNESS = (E * I) * np.array([[12 / L ** 3,  6 / L ** 2, -12 / L ** 3,  6 / L ** 2],
                                                            [6 / L ** 2, 4 / L, -6 / L ** 2, 2 / L],
                                                            [-12 / L ** 3, -6 / L ** 2, 12 / L ** 3, -6 / L ** 2],
                                                            [6 / L ** 2, 2 / L, -6 / L ** 2, 4 / L]])
    elif TYPE_ELEMENT == "FRAME_2D_6DOF":
            L = ARRAY_ELEMENT_GEOMETRY[0]
            A = ARRAY_ELEMENT_GEOMETRY[3]
            I = ARRAY_ELEMENT_GEOMETRY[5]
            E = ARRAY_ELEMENT_MATERIAL[0]
            C1 = A * E / L
            C2 = E * I / (L **3)
            MATRIX_ELEMENTAR_STIFFNESS = np.array([[C1, 0, 0, -C1, 0, 0],
                                                  [0, 12 * C2, 6 * C2 * L, 0, -12 * C2, 6 * C2 * L],
                                                  [0, 6 * C2 * L, 4 * C2 * L ** 2, 0, -6 * C2 * L, 2 * C2 * L ** 2],
                                                  [-C1, 0, 0, C1, 0, 0],
                                                  [0, -12 * C2, -6 * C2 * L, 0, 12 * C2, -6 * C2 * L],
                                                  [0, 6 * C2 * L, 2 * C2 * L ** 2, 0, -6 * C2 * L, 4 * C2 * L **2]])
    return MATRIX_ELEMENTAR_STIFFNESS

def ROTATION(TYPE_ELEMENT, ARRAY_ELEMENT_GEOMETRY):
    SINN = ARRAY_ELEMENT_GEOMETRY[1]
    COSS = ARRAY_ELEMENT_GEOMETRY[2]
    if TYPE_ELEMENT == "BAR_2D_1DOF":
        MATRIX_TRANSFORM = np.array([[1 , 0],
                                    [0, 1]])
    elif TYPE_ELEMENT == "TRUSS_2D_4DOF":
        MATRIX_TRANSFORM = np.array([[COSS, SINN, 0, 0],
                                     [-SINN, COSS, 0, 0],
                                     [ 0, 0, COSS, SINN],
                                     [0, 0, -SINN, COSS]])
    elif TYPE_ELEMENT == "BEAM_2D_4DOF":
        MATRIX_TRANSFORM = np.array([[1, 0, 0, 0],
                                     [0, 1, 0, 0],
                                     [0, 0, 1, 0],
                                     [0, 0, 0, 1]])
    elif TYPE_ELEMENT == "FRAME_2D_6DOF":
        MATRIX_TRANSFORM = np.array([[COSS, SINN, 0, 0, 0, 0],
                                     [-SINN,COSS, 0, 0, 0, 0],
                                     [0, 0, 1, 0, 0, 0],
                                     [0, 0, 0, COSS, SINN, 0],
                                     [0, 0, 0, -SINN, COSS, 0],
                                     [0, 0, 0, 0, 0, 1]])
    return MATRIX_TRANSFORM

def GLOBAL_DOF_ELEMENT(TOTAL_NODES_ELEMENT, TOTAL_DOF_NODE, MATRIX_GLOBAL_DOF, MATRIX_ELEMENT_PROPERTIES, I_ELEMENT):
    ARRAY_GLOBAL_DOF = []
    for I_COUNT in range(TOTAL_NODES_ELEMENT):
        NODE = MATRIX_ELEMENT_PROPERTIES[I_ELEMENT, I_COUNT]
        for J_COUNT in range(TOTAL_DOF_NODE):
            DOF_VALUE = int(MATRIX_GLOBAL_DOF[NODE, J_COUNT])
            ARRAY_GLOBAL_DOF.append(DOF_VALUE)
    return ARRAY_GLOBAL_DOF

def GLOBAL_STIFFNESS(TOTAL_DOF, ARRAY_GLOBAL_DOF, MATRIX_ELEMENTAR_STIFFNESS):
    MATRIX_GLOBAL_STIFFNESS = np.zeros((TOTAL_DOF, TOTAL_DOF))
    for I_COUNT, I_VALUE in enumerate(ARRAY_GLOBAL_DOF):
        for J_COUNT, J_VALUE in enumerate(ARRAY_GLOBAL_DOF):
            MATRIX_GLOBAL_STIFFNESS[I_VALUE, J_VALUE] = MATRIX_GLOBAL_STIFFNESS[I_VALUE, J_VALUE] + MATRIX_ELEMENTAR_STIFFNESS[I_COUNT, J_COUNT]
    return MATRIX_GLOBAL_STIFFNESS

def CONDENSE_FREE_GLOBAL_STIFFNESS(MATRIX_GLOBAL_STIFFNESS, ARRAY_FREE_DOF, TOTAL_FREE_DOF):
    MATRIX_FREE_FREE_GLOBAL_STIFFNESS = np.zeros((TOTAL_FREE_DOF, TOTAL_FREE_DOF))
    for I_COUNT in range(TOTAL_FREE_DOF):
        FREE_DOF_LINE = ARRAY_FREE_DOF[I_COUNT]
        for J_COUNT in range(TOTAL_FREE_DOF):
            FREE_DOF_COLUMN = ARRAY_FREE_DOF[J_COUNT]
            MATRIX_FREE_FREE_GLOBAL_STIFFNESS[I_COUNT, J_COUNT] = MATRIX_GLOBAL_STIFFNESS[FREE_DOF_LINE, FREE_DOF_COLUMN]
    return MATRIX_FREE_FREE_GLOBAL_STIFFNESS

def CONDENSE_PRESCRIBED_FREE_GLOBAL_STIFFNESS(MATRIX_GLOBAL_STIFFNESS, ARRAY_FREE_DOF, TOTAL_FREE_DOF, ARRAY_PRESCRIBED_DOF, TOTAL_PRESCRIBED_DOF):
    MATRIX_PRESCRIBED_FREE_STIFFNESS = np.zeros((TOTAL_PRESCRIBED_DOF, TOTAL_FREE_DOF))
    for I_COUNT in range(TOTAL_PRESCRIBED_DOF):
        PRECRIBED_DOF_LINE = ARRAY_PRESCRIBED_DOF[I_COUNT]
        for J_COUNT in range(TOTAL_FREE_DOF):
            FREE_DOF_COLUMN = ARRAY_FREE_DOF[J_COUNT]
            MATRIX_PRESCRIBED_FREE_STIFFNESS[I_COUNT, J_COUNT] = MATRIX_GLOBAL_STIFFNESS[PRECRIBED_DOF_LINE, FREE_DOF_COLUMN]
    return MATRIX_PRESCRIBED_FREE_STIFFNESS

def CONDENSE_FREE_GLOBAL_FORCES(ARRAY_EXTERNAL_GLOBAL_FORCES, ARRAY_FREE_DOF, TOTAL_FREE_DOF):
    ARRAY_FREE_GLOBAL_FORCES = np.zeros((TOTAL_FREE_DOF, 1))
    for I_COUNT in range(TOTAL_FREE_DOF):
        FREE_DOF = ARRAY_FREE_DOF[I_COUNT]
        ARRAY_FREE_GLOBAL_FORCES[I_COUNT, 0] = ARRAY_EXTERNAL_GLOBAL_FORCES[FREE_DOF, 0]
    return ARRAY_FREE_GLOBAL_FORCES

def CONDENSE_PRESCRIBED_GLOBAL_DISPLACEMENT(ARRAY_PRESCRIBED_DOF_VALUE, TOTAL_PRESCRIBED_DOF):
    ARRAY_PRESCRIBED_GLOBAL_DISPLACEMENT = np.zeros((TOTAL_PRESCRIBED_DOF, 1))
    for I_COUNT in range(TOTAL_PRESCRIBED_DOF):
        PRESCRIBED_DOF_VALUE = ARRAY_PRESCRIBED_DOF_VALUE[I_COUNT]
        ARRAY_PRESCRIBED_GLOBAL_DISPLACEMENT[I_COUNT, 0] = PRESCRIBED_DOF_VALUE
    return ARRAY_PRESCRIBED_GLOBAL_DISPLACEMENT

def ASSEMBLY_TOTAL_DISPLACEMENT(ARRAY_FREE_GLOBAL_DISPLACEMENT, ARRAY_PRESCRIBED_GLOBAL_DISPLACEMENT, TOTAL_DOF, ARRAY_PRESCRIBED_DOF, ARRAY_FREE_DOF):
    ARRAY_GLOBAL_DISPLACEMENT = np.zeros((TOTAL_DOF, 1))
    for I_COUNT, I_VALUE in enumerate(ARRAY_PRESCRIBED_DOF):
        DOF_DISPLACEMENT_VALUE = ARRAY_PRESCRIBED_GLOBAL_DISPLACEMENT[I_COUNT, 0]
        ARRAY_GLOBAL_DISPLACEMENT[I_VALUE, 0] = DOF_DISPLACEMENT_VALUE
    for J_COUNT, J_VALUE in enumerate(ARRAY_FREE_DOF):
        DOF_DISPLACEMENT_VALUE = ARRAY_FREE_GLOBAL_DISPLACEMENT[J_COUNT, 0]
        ARRAY_GLOBAL_DISPLACEMENT[J_VALUE, 0] = DOF_DISPLACEMENT_VALUE
    return ARRAY_GLOBAL_DISPLACEMENT


# In[2]:


import numpy as np
typeElement              = "TRUSS_2D_4DOF"
typeSolution             = "CONDENSE_PROCEDURE"
totalNodes               = 6
coordinatesMatrix        = np.array([[18.28, 9.14], [18.28, 0.0], [9.14, 9.14], [9.14, 0.0], [0.0, 9.14], [0.0, 0.0]])
materialMatrix           = np.array([[69000000, 0.1, 0.2],[69000000, 0.1, 0.2]])
totalMaterials           = 2
sectionMatrix            = np.array([[0.01, 0.01, 0.02]])
totalSection             = 1
totalNodalForces         = 2
TotalElements            = 10
ElementsPropertiesMatrix       = np.array([[4, 2, 1, 0], [2, 0, 1, 0], [5, 3, 1, 0], [3, 1, 1, 0], [3, 2, 1, 0], [1, 0, 1, 0], [4, 3, 1, 0], [5, 2, 1, 0], [3, 0, 1, 0], [2, 1, 1, 0]]) 
prescribedDisplacementMatrix   = np.array([[4, 0, 0],[4, 1, 0],[5, 0, 0],[5, 1, 0]])
externalNodalForcesMatrix      = np.array([[1, 1, -450],[3, 1, -450]])
totalNodalForces         = 2

totalNodes               = 3
coordinatesMatrix        = np.array([[1.0, 0.0], [0.0, 0.0], [1.0, 1.0]])
materialMatrix           = np.array([[1.0, 1.0, 0]])
totalMaterials           = 1
sectionMatrix            = np.array([[1.0, 0, 0]])
totalSection             = 1
totalNodalForces         = 1
TotalElements            = 2
ElementsPropertiesMatrix       = np.array([[0, 2, 0, 0], [1, 2, 0, 0]]) 
prescribedDisplacementMatrix   = np.array([[0, 0, 0], [0, 1, 0], [1, 0, 0], [1, 1, 0]])
externalNodalForcesMatrix      = np.array([[2, 0, 10]])
totalNodalForces         = 1



print("     ")
print("linha")
[totalDOFNode, totalNodesElement, totalDOFElement,  columnMaterial, columnSection, arrayDOFActive, localDOF] = INDEX_ASSEMBLY(typeElement)

print("totalDOFNode", totalDOFNode)
print("totalNodesElement", totalNodesElement)
print("totalDOFElement", totalDOFElement)
print("columnMaterial", columnMaterial)
print("columnSecao", columnSection)
print("arrayDOFActive", arrayDOFActive)
print("Local DOF", localDOF)

print("     ")
print("linha")
nodesGlobalDOF = GLOBAL_DOF_ASSEMBLY(typeElement, totalDOFNode, totalNodes, TotalElements)
print("nodesGlobalDOF")
print(nodesGlobalDOF)

print("     ")
print("linha")
[totalDOFArray, totalDOF] = TOTAL_DEGREE_FREEDOM(totalDOFNode, totalNodes)
print("totalDOF", totalDOF)
print("totalDOFArray")
print(totalDOFArray)

print("     ")
print("linha")
[numberDOFPrescribed, valuePrescribedDOFArray, totalPrescribedDOFArray] = PRESCRIBED_DEGREE_FREEDOM(prescribedDisplacementMatrix, nodesGlobalDOF, totalDOFNode)
print("numberDOFPrescribed", numberDOFPrescribed)
print("valuePrescribedDOFArray")
print(valuePrescribedDOFArray)
print("totalPrescribedDOFArray")
print(totalPrescribedDOFArray)

print("     ")
print("linha")
[numberDOFFree, totalFreeDOFArray] = FREE_DEGREE_FREEDOM(totalPrescribedDOFArray, totalDOFArray)
print("numberDOFFree", numberDOFFree)
print("totalFreeDOFArray")
print(totalFreeDOFArray)

print("     ")
print("linha")
nodalForceContribuition = NODAL_EXTERNAL_LOAD(externalNodalForcesMatrix , totalNodalForces, totalDOF, nodesGlobalDOF)
print("nodalForceContribuition")
print(nodalForceContribuition)

print("     ")
print("linha")
if totalElementForces > 0
    nodalElementContribuition =  np.zeros((totalDOF, 1))
    for iElement in range(totalElementForces):
        loadElementContribuition = DISTRIBUTED_EXTERNAL_LOAD(totalDOF, ARRAY_ELEMENT_GEOMETRY, TYPE_ELEMENT, TYPE_EXTERNAL_LOAD, MATRIX_DISTRIBUTED_EXTERNAL_LOAD, TOTAL_NODES_ELEMENT, MATRIX_ELEMENT_PROPERTIES, MATRIX_GLOBAL_DOF)
        nodalElementContribuition = nodalElementContribuition + loadElementContribuition

print("     ")
print("linha")
structureStiffness = np.zeros((totalDOF, totalDOF))
for iElement in range(TotalElements):
    
    print("ielement", iElement)
    materialsElement = MATERIALS_PROPRETIES(ElementsPropertiesMatrix, materialMatrix, iElement, columnMaterial)
    
    print("material")
    print(materialsElement)
    sectionsElement = GEOMETRIC_PROPRETIES(coordinatesMatrix, ElementsPropertiesMatrix, sectionMatrix, iElement, columnSection)
    
    print("section")
    print(sectionsElement)
    
    elementarStiffnessLocalAxis = ELEMENTAR_STIFFNESS(typeElement, sectionsElement, materialsElement)
    
    print("KelLocal")
    print(elementarStiffnessLocalAxis)
    
    rotationMatrix = ROTATION(typeElement, sectionsElement)
    print("Rotation")
    print(rotationMatrix)
    
    elementarStiffnessGlobalAxis = np.dot(np.dot(rotationMatrix.T, elementarStiffnessLocalAxis), rotationMatrix)
    print("KelGlobal")
    print(elementarStiffnessGlobalAxis)
    
    globalDOF = GLOBAL_DOF_ELEMENT(totalNodesElement, totalDOFNode, nodesGlobalDOF, ElementsPropertiesMatrix, iElement)
    elementarContribuitionStructureStiffness = GLOBAL_STIFFNESS(totalDOF, globalDOF, elementarStiffnessGlobalAxis)
    structureStiffness = structureStiffness + elementarContribuitionStructureStiffness
    print("KelGlobal Contribuition")
    print(structureStiffness)    

print("Final Stiffness")
print(structureStiffness)    
if typeSolution == "CONDENSE_PROCEDURE":
    arrayGlobalForceFreeFree = CONDENSE_FREE_GLOBAL_FORCES(nodalForceContribuition, totalFreeDOFArray, numberDOFFree)
    print("Fff")
    print(arrayGlobalForceFreeFree)
    structureStiffnessFreeFree = CONDENSE_FREE_GLOBAL_STIFFNESS(structureStiffness, totalFreeDOFArray, numberDOFFree)
    print("Kff")
    print(structureStiffnessFreeFree)
    structureStiffnessPrescribedFree = CONDENSE_PRESCRIBED_FREE_GLOBAL_STIFFNESS(structureStiffness, totalFreeDOFArray, numberDOFFree, totalPrescribedDOFArray, numberDOFPrescribed)
    print("Kpf")
    print(structureStiffnessPrescribedFree)    
    arrayGlobalDisplacementPrescribed = CONDENSE_PRESCRIBED_GLOBAL_DISPLACEMENT(valuePrescribedDOFArray, numberDOFPrescribed)
    print("Up")
    print(arrayGlobalDisplacementPrescribed)
    structureStiffnessFreePrescribed = structureStiffnessPrescribedFree.T
    kurDotUr = np.dot(structureStiffnessFreePrescribed, arrayGlobalDisplacementPrescribed)
    frMinuskurDotUr = arrayGlobalForceFreeFree - kurDotUr
    inverseKff = np.linalg.inv(structureStiffnessFreeFree)
    arrayGlobalDisplacementFree = np.dot(inverseKff, frMinuskurDotUr)
    print("Uf")
    print(arrayGlobalDisplacementFree)
elif typeSolution == 'ZERO_AND_ONE_PROCEDURE':
    pass

arrayGlobalDisplacement = ASSEMBLY_TOTAL_DISPLACEMENT(arrayGlobalDisplacementFree, arrayGlobalDisplacementPrescribed, totalDOF, totalPrescribedDOFArray, totalFreeDOFArray)
print("Global displacement")
print(arrayGlobalDisplacement)
    


# In[ ]:





# In[ ]:




