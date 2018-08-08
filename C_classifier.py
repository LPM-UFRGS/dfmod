#-------------------------------------------------------------------------------
#
# Signed Distances Function Interpolator
# **************************************
#
# This SGeMS plugin interpolates (OK) the signed distance function calculated
# for each data and rock type, and creates a geologic model based on the minimum
# estimated distance.
#
# AUTHORS: Flávio A. N. Amarante and Roberto Mentzingen Rolo
#
#-------------------------------------------------------------------------------

#!/bin/python
import sgems
import math
import numpy as np
import random
import copy
from scipy.stats import norm
import sys

#Creates a randon path given the size of the grid
def random_path(prop):
    nodes_not_nan = []
    for i in range(len(prop)):
        if not math.isnan(prop[i]):
            nodes_not_nan.append(i)
    random.shuffle(nodes_not_nan)
    return nodes_not_nan

#Calculates the proportions of variables on a list
def proportion(var, RT):
    rock_types =[]
    target_prop = []
    for k in range(len(RT)):
        target_prop.append(0)
        rock_types.append(int(RT[k].split('RT_')[-1]))
    rock_types.sort()
    var_not_nan = []
    for i in var:
        if not math.isnan(i):
            var_not_nan.append(i)
    for i in range(len(rock_types)):
        target_prop[i] = float(var.count(rock_types[i]))/len(var_not_nan)
    return target_prop

def create_variable(grid, name, list):
    lst_props_grid = sgems.get_property_list(grid)
    prop_final_data_name = name

    if (prop_final_data_name in lst_props_grid):
        flag = 0
        i = 1
        while (flag == 0):
            test_name = prop_final_data_name + '-' + str(i)
            if (test_name not in lst_props_grid):
                flag = 1
                prop_final_data_name = test_name
            i = i + 1

    sgems.set_property(grid, prop_final_data_name, list)

def classifier(grid, sim_val, interpol_val, C, var_n):
    #sim_val = sgems.get_property(grid, sim_val)
    #interpol_val = sgems.get_property(grid, interpol_val)

    trans_sim = sim_val
    for n,i in enumerate(sim_val):
        if not math.isnan(i):
            trans_sim[n] = (2*C*norm.cdf(i)-C)

    class_blocks = trans_sim
    for n,i in enumerate(sim_val):
        if not math.isnan(i):
            if i > interpol_val[n]:
                class_blocks[n] = 1
            if i < interpol_val[n]:
                class_blocks[n] = 0
    for n,i in enumerate(interpol_val):
        if i < -C:
            class_blocks[n] = 1
        if i > C:
            class_blocks[n] = 0

    FN = "Class_C_("+str(var_n)+")_C_"+ str(C)
    create_variable(grid, FN, class_blocks)

#Transform i,j,k in n
def ijk_in_n(grid, i, j, k):
    dims = sgems.get_dims(grid)
    n = k*dims[0]*dims[1]+j*dims[0]+i
    return n

#Crestes a list with indices of the neighbors valid blocks
def neighb(grid, indice):
        ijk = sgems.get_ijk(grid, indice)
        neighborhood = []
        for i in range(ijk[0]-1,ijk[0]+2):
            for j in range(ijk[1]-1,ijk[1]+2):
                for k in range(ijk[2]-1,ijk[2]+2):
                    ijk_blk = [i,j,k]
                    neighborhood.append(ijk_blk)
        dims = sgems.get_dims(grid)
        neighborhood_cp = copy.copy(neighborhood)
        for i in neighborhood_cp:
            if dims[2] == 1:
                if i[0] < 0 or i[1] < 0:
                    neighborhood.remove(i)
                elif i[0] > (dims[0] - 1) or i[1] > (dims[1] - 1):
                    neighborhood.remove(i)
                elif i[2] != 0:
                    neighborhood.remove(i)
                elif i == sgems.get_ijk(grid, indice):
                    neighborhood.remove(i)
            else:
                if i[0] < 0 or i[1] < 0 or i[2] < 0:
                    neighborhood.remove(i)
                elif i[0] > (dims[0] - 1) or i[1] > (dims[1] - 1) or i[2] > (dims[2] - 1):
                    neighborhood.remove(i)
                elif i == sgems.get_ijk(grid, indice):
                    neighborhood.remove(i)
        neighborhood_n = []
        for i in neighborhood:
            neighborhood_n.append(ijk_in_n(grid,i[0],i[1],i[2]))
        return neighborhood_n


# Shows every parameter of the plugin in the command pannel
def read_params(a, j=''):
    for i in a:
        if (type(a[i]) != type({'a': 1})):
            print j + "['" + str(i) + "']=" + str(a[i])
        else:
            read_params(a[i], j + "['" + str(i) + "']")

class C_classifier:
    def __init__(self):
        pass

    def initialize(self, params):
        self.params = params
        return True

    def execute(self):

        #Execute the funtion read_params
        #read_params(self.params)
        #print self.params'''

        #temp = sys.stdout
        #sys.stdout = sys.stderr
        #sys.stderr = temp

        # Get the grid and rock type propery
        grid = self.params['propertyselectornoregion']['grid']
        prop = self.params['propertyselectornoregion']['property']

        # Get the X, Y and Z coordinates and RT property
        X = sgems.get_property(grid, '_X_')
        Y = sgems.get_property(grid, '_Y_')
        Z = sgems.get_property(grid, '_Z_')
        RT_data = sgems.get_property(grid, prop)

        # Getting properties
        grid_krig = self.params['gridinter']['value']
        grid_var = self.params['gridselectorbasic']['value']
        props = (self.params['orderedpropertyselector']['value']).split(';')
        n_var = int(self.params['indicator_regionalization_input']['number_of_indicator_group'])
        n_prop = int(self.params['orderedpropertyselector']['count'])
        min_cond = self.params['spinBox_2']['value']
        max_cond = self.params['spinBox']['value']
        C_val = float(self.params['c_par']['value'])

        # Error messages
        if len(grid_var) == 0 or len(grid_krig) == 0:
            print 'Select the variables'
            return False

        if n_var != n_prop:
            print 'Number of variables and number of variograms models are diferent.'
            return False

        # Creating an empty list to store the interpolated distances
        SG_OK_list = []

        # Loop in every variable
        for i in xrange(0, n_var):

            #....... Interpolation .........#

            #Getting variables Interpolation
            prop_HD = props[i]
            prop_name = "Interpolated_" + str(prop_HD)
            prop_name_var = "Interpolated_" + str(prop_HD) + ' krig_var'
            var_str = ''
            indicator_group = "Indicator_group_" + str(i + 1)
            elipsoide = self.params['ellipsoidinput']['value']
            n_struct = int(self.params['indicator_regionalization_input'][indicator_group]['Covariance_input']['structures_count'])

            # Error message
            if n_struct == 0:
                print 'Variogram have no structures'
                #continue
                return False

            # Loop in every variogram structure - Interpolation
            for j in xrange(0, n_struct):

                # Getting variogram parameters
                Structure = "Structure_" + str(j + 1)
                cov_type = self.params['indicator_regionalization_input'][indicator_group]['Covariance_input'][Structure]['Two_point_model']['type']
                cont = self.params['indicator_regionalization_input'][indicator_group]['Covariance_input'][Structure]['Two_point_model']['contribution']

                if cov_type == 'Nugget Covariance':
                    # Writing variogram parameters on a variable in nugget effect case
                    var_str = var_str + '<{} type="{}">  <Two_point_model  contribution="{}"  type="{}"   >    </Two_point_model>    </Structure_1> '.format(Structure, 'Covariance', cont, cov_type, Structure)
                    print var_str
                else:
                    range1 = self.params['indicator_regionalization_input'][indicator_group]['Covariance_input'][Structure]['Two_point_model']['ranges']['range1']
                    range2 = self.params['indicator_regionalization_input'][indicator_group]['Covariance_input'][Structure]['Two_point_model']['ranges']['range2']
                    range3 = self.params['indicator_regionalization_input'][indicator_group]['Covariance_input'][Structure]['Two_point_model']['ranges']['range3']
                    rake = self.params['indicator_regionalization_input'][indicator_group]['Covariance_input'][Structure]['Two_point_model']['angles']['rake']
                    dip = self.params['indicator_regionalization_input'][indicator_group]['Covariance_input'][Structure]['Two_point_model']['angles']['dip']
                    azimuth = self.params['indicator_regionalization_input'][indicator_group]['Covariance_input'][Structure]['Two_point_model']['angles']['azimuth']

                    # Writing variogram parameters on a variable in other cases
                    var_str = var_str + '<{} type="{}">  <Two_point_model  contribution="{}"  type="{}"   >      <ranges range1="{}"  range2="{}"  range3="{}"   />      <angles azimuth="{}"  dip="{}"  rake="{}"   />    </Two_point_model>    </{}> '.format(Structure, 'Covariance', cont, cov_type, range1, range2, range3, azimuth, dip, rake, Structure)

            # Calling ordinary kriging for each variable, using the variograms parameters above
            sgems.execute('RunGeostatAlgorithm  kriging::/GeostatParamUtils/XML::<parameters>  <algorithm name="kriging" />     <Variogram  structures_count="{}" >    {}  </Variogram>    <ouput_kriging_variance  value="0"  />     <output_n_samples_  value="0"  />     <output_average_distance  value="0"  />     <output_sum_weights  value="0"  />     <output_sum_positive_weights  value="0"  />     <output_lagrangian  value="0"  />     <Nb_processors  value="-2"  />    <Grid_Name value="{}" region=""  />     <Property_Name  value="{}" />     <Hard_Data  grid="{}"   property="{}"   region=""  />     <Kriging_Type  type="Ordinary Kriging (OK)" >    <parameters />  </Kriging_Type>    <do_block_kriging  value="1"  />     <npoints_x  value="5" />     <npoints_y  value="5" />     <npoints_z  value="5" />     <Min_Conditioning_Data  value="{}" />     <Max_Conditioning_Data  value="{}" />     <Search_Ellipsoid  value="{}" />    <AdvancedSearch  use_advanced_search="0"></AdvancedSearch>  </parameters>'.format(n_struct, var_str, grid_krig, prop_name, grid_var, prop_HD, min_cond, max_cond, elipsoide))

            SG_OK_list = sgems.get_property(grid_krig, prop_name)

        	# Uncertainty Zone
            zone_name = "UZ_C_"+str(C_val)
            zone = np.array(SG_OK_list)
            m1 = zone < C_val
            m2 = zone > -C_val
            mask = np.logical_and(m1,m2)
            uzz = mask.astype('int').tolist()
            sgems.set_region(grid_krig,zone_name,uzz)

        	#......... Simulation .........#

        	#Creating an empty list to store the simulated distances
            SG_sim_list = []

            # Getting variables Simulation
            grid_sim_region = zone_name
            prop_name_sim = "Sim_C_"+ str(C_val)
            var_str_s = ''
            nb_realz = self.params['spinBox_3']['value']
            max_cond_hard_data = self.params['Max_Conditioning_Data']['value']
            max_cond_prev_sim_data = self.params['Max_Conditioning_Simul_Data']['value']
            number_pro = self.params['spinBox_6']['value']
            sim_group = "Indicator_group_" + str(i + 1)
            elipsoide_sim = self.params['Search_Ellipsoid_Sim']['value']
            n_struct_s = int(self.params['sim_regionalization_input'][sim_group]['Covariance_input']['structures_count'])
            HD = ""
            p_HD = ""

            #Error message
            if n_struct_s == 0:
                print 'Variogram have no structures'
                return False

      		# Loop in every variogram structure
            for j in xrange(0, n_struct):

                # Getting variogram parameters
                Structure_s = "Structure_" + str(j + 1)
                cov_type_s = self.params['sim_regionalization_input'][sim_group]['Covariance_input'][Structure_s]['Two_point_model']['type']
                cont_s = self.params['sim_regionalization_input'][sim_group]['Covariance_input'][Structure_s]['Two_point_model']['contribution']

                if cov_type_s == 'Nugget Covariance':
                    #Writing variogram parameters on a variable in nugget effect case
                    var_str_s = var_str + '<{} type="{}">  <Two_point_model  contribution="{}"  type="{}"   >    </Two_point_model>    </Structure_1> '.format(Structure, 'Covariance', cont, cov_type, Structure)

                else:
                    srange1 = self.params['sim_regionalization_input'][sim_group]['Covariance_input'][Structure_s]['Two_point_model']['ranges']['range1']
                    srange2 = self.params['sim_regionalization_input'][sim_group]['Covariance_input'][Structure_s]['Two_point_model']['ranges']['range2']
                    srange3 = self.params['sim_regionalization_input'][sim_group]['Covariance_input'][Structure_s]['Two_point_model']['ranges']['range3']
                    srake = self.params['sim_regionalization_input'][sim_group]['Covariance_input'][Structure_s]['Two_point_model']['angles']['rake']
                    sdip = self.params['sim_regionalization_input'][sim_group]['Covariance_input'][Structure_s]['Two_point_model']['angles']['dip']
                    sazimuth = self.params['sim_regionalization_input'][sim_group]['Covariance_input'][Structure_s]['Two_point_model']['angles']['azimuth']

                    # Writing variogram parameters on a variable in other cases
                    var_str_s = var_str + '<{} type="{}">  <Two_point_model  contribution="{}"  type="{}"   >      <ranges range1="{}"  range2="{}"  range3="{}"   />      <angles azimuth="{}"  dip="{}"  rake="{}"   />    </Two_point_model>    </{}> '.format(Structure_s, 'Covariance', cont_s, cov_type_s, srange1, srange2, srange3, sazimuth, sdip, srake, Structure_s)

            # Calling sgsim, using the variograms parameters above
            sgems.execute('RunGeostatAlgorithm  sgsim::/GeostatParamUtils/XML::<parameters>  <algorithm name="sgsim" />     <Simulation_seed type="clock"/>  <Path type="random"><Seed type="clock"/></Path>  <Nb_processors  value="{}"  />    <Grid_Name value="{}" region="{}"  />    <Assign_Hard_Data  value="0"  />     <Property_Name value="{}" reuse_property="0" />  <Nb_Realizations  value="{}" />     <Use_LVM  value="0"  />     <Max_Conditioning_Data  value="{}" />     <Max_Conditioning_Simul_Data  value="{}" />     <Search_Ellipsoid  value="{}" />    <AdvancedSearch  use_advanced_search="0"></AdvancedSearch>     <Hard_Data  grid="{}"   property="{}"   region=""  />     <Use_Target_Histogram value="0" /> <Export_nscore value="1" /> <Use_break_tie_index value="0" /> <nonParamCdf ref_in_distribution ="0" ref_on_file ="0" ref_on_grid ="0" break_ties_indices ="0" grid ="{}" region ="" property ="{}" is_weight = "0"><LTI_type function ="No extrapolation" extreme ="0" omega ="" /> <UTI_type function ="No extrapolation" extreme ="0" omega ="" /> </nonParamCdf>      <Covariance_input  structures_count="{}" > {} </Covariance_input>  </parameters>'.format(number_pro, grid_krig, grid_sim_region, prop_name_sim, nb_realz, max_cond_hard_data, max_cond_prev_sim_data, elipsoide_sim, HD, p_HD, grid_var, prop_HD, n_struct_s, var_str_s))

            Simulations = sgems.get_properties_in_group(grid_krig,prop_name_sim)

            for prop in Simulations:
                ss = sgems.get_property(grid_krig, prop)
                SG_sim_list.append(ss)

            #classifying categories
            for i in SG_sim_list:
                classifier(grid_krig,i,SG_OK_list,C_val,prop_HD)

			#Deleting interpolation and simulations
			if self.params['Keep_in']['value']=='0':
                sgems.execute('DeleteObjectProperties  {}::{}'.format(grid_krig, prop_name))
            if self.params['keep_sim']['value']=='0':
			    sgems.execute('DeleteObjectProperties  {}::{}'.format(grid_krig, prop_name_sim))

        return True

    def finalize(self):

        return True

    def name(self):

        return "C_classifier"

################################################################################
def get_plugins():
    return ["C_classifier"]
