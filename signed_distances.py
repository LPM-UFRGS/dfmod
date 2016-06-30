#-------------------------------------------------------------------------------
#
# Signed Distance Function Calculator
# ***********************************
#
# This SGeMS plugin calculates the anisotropic signed distances for each data
# point and each rock type
#
# AUTHOR: Roberto Mentzingen Rolo
#
#-------------------------------------------------------------------------------

#!/bin/python
import sgems
import math
import numpy as np

#Calculates the distances
def dist(x1, y1, z1, x2, y2, z2):
    return math.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2)

#Defines the rotation and dilatation matrices
def rot(range1, range2, range3, azimuth, dip, rake, vetor):

    if azimuth >= 0 and azimuth <=270:
	alpha = math.radians(90-azimuth)
    else:
	alpha = math.radians(450-azimuth)
    beta = -math.radians(dip)
    phi = math.radians(rake)

    rot_matrix = np.zeros((3,3))

    rot_matrix[0,0] = math.cos(beta)*math.cos(alpha)
    rot_matrix[0,1] = math.cos(beta)*math.sin(alpha)
    rot_matrix[0,2] = -math.sin(beta)
    rot_matrix[1,0] = (range1/range2)*(-math.cos(phi)*math.sin(alpha)+math.sin(phi)*math.sin(beta)*math.cos(alpha))
    rot_matrix[1,1] = (range1/range2)*(math.cos(phi)*math.cos(alpha)+math.sin(phi)*math.sin(beta)*math.sin(alpha))
    rot_matrix[1,2] = (range1/range2)*(math.sin(phi)*math.cos(beta))
    rot_matrix[2,0] = (range1/range3)*(math.sin(phi)*math.sin(alpha)+math.cos(phi)*math.sin(beta)*math.cos(alpha))
    rot_matrix[2,1] = (range1/range3)*(-math.sin(phi)*math.cos(alpha)+math.cos(phi)*math.sin(beta)*math.sin(alpha))
    rot_matrix[2,2] = (range1/range3)*(math.cos(phi)*math.cos(beta))

    vetor = np.array(vetor)

    return np.dot(rot_matrix, vetor)

#Transform the data with the ratation/dilatation matrices
def anis_search(X, Y, Z, range1, range2, range3, azimuth, dip, rake):

    X_linha = []
    Y_linha = []
    Z_linha = []

    for i in range(len(X)):
	vet = [X[i],Y[i],Z[i]]

	vet_rot = rot(range1, range2, range3, azimuth, dip, rake, vet)

	X_linha.append(vet_rot[0])
	Y_linha.append(vet_rot[1])
	Z_linha.append(vet_rot[2])

    return X_linha, Y_linha, Z_linha

#Shows every parameter of the plugin in the command pannel
def read_params(a,j=''):
  for i in a:
    if (type(a[i])!=type({'a':1})):
      print j+"['"+str(i)+"']="+str(a[i])
    else:
      read_params(a[i],j+"['"+str(i)+"']")

class signed_distances:
    def __init__(self):
        pass

    def initialize(self, params):
        self.params = params
        return True

    def execute(self):

        '''#Execute the funtion read_params
        read_params(self.params)
        print self.params'''

        #Get the grid and rock type propery
        grid = self.params['propertyselectornoregion']['grid']
        prop = self.params['propertyselectornoregion']['property']

        #Error message
        if len(grid) == 0 or len(prop) == 0:
            print 'Select the rocktype property'
            return False

        #Get the X, Y and Z coordinates
        X = sgems.get_property(grid, '_X_')
        Y = sgems.get_property(grid, '_Y_')
        Z = sgems.get_property(grid, '_Z_')
        RT = sgems.get_property(grid, prop)

        elipsoide = self.params['ellipsoidinput']['value']
        elipsoide_split = elipsoide.split()

        range1 = float(elipsoide_split[0])
        range2 = float(elipsoide_split[1])
        range3 = float(elipsoide_split[2])

        azimuth = float(elipsoide_split[3])
        dip = float(elipsoide_split[4])
        rake = float(elipsoide_split[5])

        X, Y, Z = anis_search(X, Y, Z, range1, range2, range3, azimuth, dip, rake)

        #Creates a list of all rock types
        rt_list = []
        for i in RT:
            if i not in rt_list and not math.isnan(i):
                rt_list.append(i)

        #Sort the rock type list in crescent order
        rt_list = [int(x) for x in rt_list]
        rt_list.sort()

        #Create a empty distance matrix
        dist_matrix = np.zeros(shape = ((len(rt_list)), (len(RT))))

        #Calculates the signed distances, and append it in the distance matrix
        for i in range(len(rt_list)):
            rock = rt_list[i]

            for j in range(len(RT)):

                if math.isnan(RT[j]):
                    dist_matrix[i][j] = float('nan')

                elif RT[j] == rock:
                    dsmin = 1.0e21

                    for k in range(len(RT)):

                        if RT[j] != RT[k] and not math.isnan(RT[k]):
                            if (dist(X[j], Y[j], Z[j], X[k], Y[k], Z[k])) < dsmin:
                                dsmin = (dist(X[j], Y[j], Z[j], X[k], Y[k], Z[k]))

                        dist_matrix[i][j] = -dsmin

                else:
                    dsmin = 1.0e21

                    for k in range(len(RT)):

                        if RT[k] == rock:
                            if (dist(X[j], Y[j], Z[j], X[k], Y[k], Z[k])) < dsmin:
                                dsmin = (dist(X[j], Y[j], Z[j], X[k], Y[k], Z[k]))

                        dist_matrix[i][j] = dsmin

        #Creates the signed distances properties
        for i in range(len(dist_matrix)):
            list = dist_matrix[i].tolist()
            sgems.set_property(grid, 'Signed_Distances_RT_' + str(rt_list[i]), list)

        return True

    def finalize(self):
        return True

    def name(self):
        return "signed_distances"

###############################################################
def get_plugins():
    return ["signed_distances"]
