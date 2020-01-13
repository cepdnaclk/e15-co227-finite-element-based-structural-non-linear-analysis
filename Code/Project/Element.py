import math

import numpy as np
from numpy.linalg import inv
from Section import *


class Element:
    k_element_initial = None  # initial stiffness matrix refering local co-ordinate system
    k_element_unbalance_force=None
    element_initial_status=True
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
            # section.f_section_resist = section.analyze([0, 0])[0]
            section.analyze([0, 0])
            # initk = section.k_section
            # print(initk)
            # break
            self.sections.put(section_id, section)

        if self.n_sections == 3:
            self.wh = [1 / 3, 4 / 3, 1 / 3]
            self.x = [-1, 0, 1]
        elif self.n_sections == 4:
            self.wh = [5 / 6, 1 / 6, 1 / 6, 5 / 6]
            self.x = [-1, - 0.447214, 0.447214, 1]
        elif self.n_sections == 5:
            self.wh = [1 / 10, 49 / 90, 32 / 45, 49 / 90, 1 / 10]
            self.x = [-1, - 0.654654, 0, 0.654654, 1]
        elif self.n_sections == 6:
            self.wh = [0.066667, 0.378475, 0.554858, 0.554858, 0.378475, 0.066667]
            self.x = [-1, - 0.765055, - 0.285232, 0.285232, 0.765055, 1]

    def rotMatrix(self):
        cos = math.cos(self.angle)
        sin = math.sin(self.angle)
        return np.array([[cos, sin, 0, 0, 0, 0],
                         [-sin, cos, 0, 0, 0, 0],
                         [0, 0, 1, 0, 0, 0],
                         [0, 0, 0, cos, sin, 0],
                         [0, 0, 0, -sin, cos, 0],
                         [0, 0, 0, 0, 0, 1]])  # should return the rotational matrix for the element

    def rigidBodyTransMatrix(self):  # L--> length of the element
        L = self.length

        return np.array([[0, 1 / L, 1, 0, -1 / L, 0],
                         [0, 1 / L, 0, 0, -1 / L, 1],
                         [-1, 0, 0, 1, 0, 0]])

    def calInitialElement_K(self):  # Condition="LOCAL" for refering local co-ordinate system & Condition="GLOBAL" for refering global co-ordinate system

        initialElementFlexibMat = 0

        for section_ in range(self.n_sections):
            Section_K = self.sections[section_].k_section
            NP = np.array([[0, 0, 1],
                           [((self.x[section_] + 1) / 2) - 1, (self.x[section_] + 1) / 2, 0]])

            fh = np.linalg.inv(Section_K)
            mat1 = np.matmul(np.transpose(NP), fh)
            mat2 = np.matmul(mat1, NP)

            mat3 = mat2 * self.wh[section_]
            mat4 = mat3 * (self.length / 2)
            initialElementFlexibMat += mat4

        k_element_initial = inv(initialElementFlexibMat)
        self.k_element_initial = k_element_initial
        mat1 = np.matmul(np.transpose(self.rigidBodyTransMatrix()), k_element_initial)

        k_element_initial_local = np.matmul(mat1, self.rigidBodyTransMatrix())

        mat2 = np.matmul(np.transpose(self.rotMatrix()), k_element_initial_local)

        k_element_initial_global = np.matmul(mat2, self.rotMatrix())

        return k_element_initial_global

    def analyze(self, tolerance):

        elementDefINCR = np.array([[self.start_node.d_x],
                                   [self.start_node.d_y],
                                   [self.start_node.dm_z],
                                   [self.end_node.d_x],
                                   [self.end_node.d_y],
                                   [self.end_node.dm_z]])

        Rot =self.rotMatrix()
        rotate = np.matmul(Rot, elementDefINCR)  # convert defINCR to local co-ordinate systme
        RB=self.rigidBodyTransMatrix()

        basicSystem = np.matmul(RB, rotate)  # remove rigid body modes (basicSystem 3x1 matrix)
        k_ =self.k_element_initial
        #mat_cal = np.matmul(k_, basicSystem)

        if self.element_initial_status==True:
            elementForceINCR = np.matmul(self.k_element_initial, basicSystem)
            self.element_initial_status=False
        else:
            elementForceINCR = np.matmul(self.k_element_initial, basicSystem)-self.k_element_unbalance_force

        for section_ in range(self.n_sections):  # newton raphson iteration

            section = self.sections[section_]
            total_sectionDefINCR_ = section.total_deformation
            total_sectionForceINCR = section.total_force

            k_section = section.k_section
            NP = np.array([[0, 0, 1], [((self.x[section_] + 1) / 2) - 1, (self.x[section_] + 1) / 2, 0]])

            # ////////////starting newton-raphson iteration////////////////////

            sectionForceINCR = np.matmul(NP, elementForceINCR)


            sectionDefINCR_ = np.matmul(inv(k_section), sectionForceINCR)

            total_sectionDefINCR_ += sectionDefINCR_
            total_sectionForceINCR += sectionForceINCR

            section.analyze(np.transpose(total_sectionDefINCR_)[0])

            unbalanceForce = total_sectionForceINCR - section.f_section_resist

            print("section level iteration starts: ")
            loop_count=0
            while (self.conditionCheck(unbalanceForce, tolerance)):
                #print("section.k_section_initial",section.k_section_initial)
                #print("unbalanceForce",unbalanceForce)
                corrective_d = np.matmul(inv(section.k_section), unbalanceForce)

                total_sectionDefINCR_ += corrective_d
                section.analyze(np.transpose(total_sectionDefINCR_)[0])
                unbalanceForce = total_sectionForceINCR - section.f_section_resist
                loop_count +=1
                #print("total_sectionDefINCR_In loop:", total_sectionDefINCR_)
            # ///////////ending newton raphson iteration/////////////////

            section.total_deformation = total_sectionDefINCR_
            section.total_force = total_sectionForceINCR
            print("Section level iterations ends =",loop_count)
        K_element = 0

        for section_ in range(self.n_sections):
            section = self.sections[section_]
            k_section=section.k_section
            NP = np.array([[0, 0, 1], [((self.x[section_] + 1) / 2) - 1, (self.x[section_] + 1) / 2, 0]])
            fh = inv(k_section)
            mat1 = np.matmul(np.transpose(NP), fh)
            mat2 = np.matmul(mat1, NP)
            mat3 = mat2 * self.wh[section_]
            mat4 = mat3 * (self.length / 2)
            K_element += mat4

        K_element = inv(K_element)
        self.k_element_initial=K_element
        element_res_force=np.matmul(K_element,basicSystem)
        #self.k_element_unbalance_force=mat_cal-element_res_force
        self.k_element_unbalance_force = elementForceINCR - element_res_force

        return np.transpose(self.rotMatrix()) @ np.transpose(
            self.rigidBodyTransMatrix()) @ K_element @ self.rigidBodyTransMatrix() @ self.rotMatrix()  # 6x6 matrix refering global co-ordinate system

    def conditionCheck(self, mat, value):
        #print(mat)
        #print(value)
        max_abs_val = abs(max(mat.min(), mat.max(), key=abs))
        #print(max_abs_val)
        if max_abs_val > value:
            return True
        else:

            return False
