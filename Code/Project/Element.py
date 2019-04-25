import math

import numpy as np
from numpy.linalg import inv
from Section import *
from log_ import *


class Element:


    def __init__(self, id, start_node, end_node, cross_section, n_sections, angle, length):
        self.id = id
        self.start_node = start_node
        self.end_node = end_node
        self.cross_section = cross_section
        self.n_sections = n_sections
        self.angle = angle
        self.length = length
        self.sections = np.empty(n_sections, dtype=Section)

        for section_id in range(n_sections):
            section = Section(section_id, cross_section)
            self.sections.put(section_id, section)

    def rotMatrix(self):
        cosVal = math.cos(self.angle)
        sinVal = math.sin(self.angle)
        return np.array([[cosVal, sinVal, 0, 0, 0, 0], [-sinVal, cosVal, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0],[0, 0, 0, cosVal, sinVal, 0], [0, 0, 0, -sinVal, cosVal, 0],[0, 0, 0, 0, 0, 1]])  # should return the rotational matrix for the element

    def rigidBodyTransMatrix(self):  # L--> length of the element
        L = self.length

        return np.array([[0, 1 / L, 1, 0, -1 / L, 0], [0, 1 / L, 0, 0, -1 / L, 1], [-1, 0, 0, 1, 0, 0]], dtype=float)

    def calInitialElement_K(self, Condition):  # Condition="LOCAL" for refering local co-ordinate system & Condition="GLOBAL" for refering global co-ordinate system

        logger.info("Call element %d initial stiffness matrix" % self.id)
        if self.n_sections == 3:
            wh = [1 / 3, 4 / 3, 1 / 3]
            x = [-1, 0, 1]
        elif self.n_sections == 4:
            wh = [5 / 6, 1 / 6, 1 / 6, 5 / 6]
            x = [-1, - 0.447214, 0.447214, 1]
        elif self.n_sections == 5:
            wh = [1 / 10, 49 / 90, 32 / 45, 49 / 90, 1 / 10]
            x = [-1, - 0.654654, 0, 0.654654, 1]
        elif self.n_sections == 6:
            wh = [0.066667, 0.378475, 0.554858, 0.554858, 0.378475, 0.066667]
            x = [-1, - 0.765055, - 0.285232, 0.285232, 0.765055, 1]

        initialElementFlexibMat = 0

        for section_ in range(self.n_sections):
            Section_K = self.sections[section_].analyze([0, 0])
            NP = np.array([[0, 0, 1], [((x[section_] + 1) / 2) - 1, (x[section_] + 1) / 2 , 0]])
            fh =inv(Section_K[1])
            mat1 = np.transpose(NP) @ fh
            mat2 = mat1 @ NP
            mat3 = mat2 * wh[section_]
            mat4=mat3 * (self.length/2)
            initialElementFlexibMat += mat4

        k_element_initial = inv(initialElementFlexibMat)
        # self.K_element_initial=np.transpose(self.rotMatrix())@np.transpose(self.rigidBodyTransMatrix())@k_element_initial@self.rigidBodyTransMatrix()@self.rotMatrix()

        if Condition == "GLOBAL":
            returnMatrix = np.transpose(self.rotMatrix()) @ np.transpose(self.rigidBodyTransMatrix()) @ k_element_initial @ self.rigidBodyTransMatrix() @ self.rotMatrix()  # 6x6 matrix refering global co-ordinate system
            logger.info("return Element %d inital stiffness matrix GLOBAL" % self.id)
            logger.debug("Element %d inital stiffness matrix GLOBAL is: " % self.id)
            logger.debug(returnMatrix)
            return returnMatrix

        elif Condition == "LOCAL":
            logger.info("return Element %d inital stiffness matrix LOCAL" % self.id)
            logger.debug("Element %d inital stiffness matrix LOCAL is: " % self.id)
            logger.debug(k_element_initial)
            return k_element_initial

    def analyze(self, tolerance):  # for the first iteration set the initial call to True

        logger.info("Element:%d Sectional level iteration running" % self.id)

        elementDefINCR = np.array([self.start_node.d_y, self.start_node.d_x, self.start_node.dm_z, self.end_node.d_y, self.end_node.d_x,self.end_node.dm_z])

        logger.info("######Element %d elementDefINCR :" %self.id)
        logger.info(elementDefINCR)

        if self.n_sections == 3:
            wh = [1 / 3, 4 / 3, 1 / 3]
            x = [-1, 0, 1]
        elif self.n_sections == 4:
            wh = [5 / 6, 1 / 6, 1 / 6, 5 / 6]
            x = [-1, - 0.447214, 0.447214, 1]
        elif self.n_sections == 5:
            wh = [1 / 10, 49 / 90, 32 / 45, 49 / 90, 1 / 10]
            x = [-1, - 0.654654, 0, 0.654654, 1]
        elif self.n_sections == 6:
            wh = [0.066667, 0.378475, 0.554858, 0.554858, 0.378475, 0.066667]
            x = [-1, - 0.765055, - 0.285232, 0.285232, 0.765055, 1]

        logger.debug("Element %d roatational matrix is:" % self.id)
        logger.debug(self.rotMatrix())
        logger.debug("Element %d rigidBody matrix is:" % self.id)
        logger.debug(self.rigidBodyTransMatrix())

        rotate = np.matmul(self.rotMatrix(), elementDefINCR)  # convert defINCR to local co-ordinate systme
        basicSystem = np.matmul(self.rigidBodyTransMatrix(), rotate)  # remove rigid body modes (basicSystem 3x1 matrix)

        #########################################################################

        K_element_intial = self.calInitialElement_K("LOCAL")  # calculate initialElement Stiffness matrix and update K_element_initial
        elementForceINCR = np.matmul(K_element_intial, basicSystem)

        logger.info("Element %d initial matrix"%self.id)
        logger.info(K_element_intial)
        logger.info("Basic system")
        logger.info(basicSystem)
        logger.info("Element %d Element force increment" % self.id)
        logger.info(elementForceINCR)

            ########################################

        # elementForceINCR ---> 3x1 matrix
        section_out_converge=True

        for section_ in range(self.n_sections):  # newton raphson iteration
            logger.info("Element %d sectional iteration running" % self.id)
            logger.info(self.sections[section_].s_h)
            NP = np.array([[0, 0, 1], [((x[section_] + 1) / 2 )- 1, (x[section_] + 1) / 2 , 0]])

            sectionForceINCR = np.matmul(NP, elementForceINCR)  # sectionForceINCR ---> 2X1 matrix

            self.sections[section_].s_h+=sectionForceINCR

            [Section_R, Section_K_Initial] = self.sections[section_].analyze([0, 0])

            section_in_converge=False

            logger.info("Element %d section %d section force before converge"%(self.id,section_))
            logger.info(self.sections[section_].s_h)
            for itr in range(1,1000):

                sectionDefINCR_=Section_K_Initial@sectionForceINCR
                self.sections[section_].e_h+=sectionDefINCR_
                cross_section_result = self.sections[section_].analyze(self.sections[section_].e_h)
                logger.info("self.sections[section_].e_h")
                logger.info(self.sections[section_].e_h)
                sectionResistingForce = cross_section_result[0]
                Section_K = cross_section_result[1]
                unbalanceForce=self.sections[section_].s_h-sectionResistingForce

                if self.conditionCheck(unbalanceForce,tolerance)==False:
                    logger.info("Element %d section %d convergence done" %(self.id,section_))
                    self.sections[section_].section_k = Section_K
                    section_in_converge=True
                    break

            if section_in_converge==False:
                logger.info("Section LEVEL convergence can not achieve ")
                section_out_converge=False
                break
        if section_out_converge==True:

            elementFlexibMat = 0  # calculate element stiffness
            logger.info("Element %d Sectional assembling running" % self.id)
            for section_ in range(self.n_sections):
                NP = np.array([[0, 0, 1], [((x[section_] + 1) / 2) - 1, (x[section_] + 1) / 2, 0]])
                fh = inv(self.sections[section_].section_k)
                mat1 = np.matmul(np.transpose(NP), fh)
                mat2 = np.matmul(mat1, NP)
                mat3 = mat2 * wh[section_]
                mat4= mat3 * (self.length/2)
                elementFlexibMat += mat4

                logger.debug("Element %d section %d section flexibility matrix is" % (self.id, section_))
                logger.debug(fh)

            K_element = inv(elementFlexibMat)

            returnMatrix = np.transpose(self.rotMatrix()) @ np.transpose(self.rigidBodyTransMatrix()) @ K_element @ self.rigidBodyTransMatrix() @ self.rotMatrix()  # 6x6 matrix refering global co-ordinate system

            return returnMatrix
        else:
            logger.info("should die")

    def conditionCheck(self, mat, value):
        max_abs_val=abs(max(mat.min(),mat.max(),key=abs))
        logger.info("Checking convergence. max_abs_val = %f"%max_abs_val)
        if max_abs_val > value:
            return True
        else:

            return False
