#!/bin/python
import sgems
import math
import numpy as np
import random

def random_path(n):
    path = range(n)
    random.shuffle(path)
    return path

def proportion(var):
    items = []
    for i in var:
        if i not in items:
            items.append(i)
    items.sort()
    target_prop = []
    for i in items:
        target_prop.append(float(var.count(i))/len(var))
    return target_prop

# Shows every parameter of the plugin in the command pannel
def read_params(a, j=''):
    for i in a:
        if (type(a[i]) != type({'a': 1})):
            print j + "['" + str(i) + "']=" + str(a[i])
        else:
            read_params(a[i], j + "['" + str(i) + "']")

class interpolator:
    def __init__(self):
        pass

    def initialize(self, params):
        self.params = params
        return True

    def execute(self):

        # Execute the funtion read_params
        read_params(self.params)
        print self.params

        # Getting properties
        grid_krig = self.params['gridselectorbasic_2']['value']
        grid_var = self.params['gridselectorbasic']['value']
        props = (self.params['orderedpropertyselector']['value']).split(';')
        n_var = int(self.params['indicator_regionalization_input']['number_of_indicator_group'])
        n_prop = int(self.params['orderedpropertyselector']['count'])
        min_cond = self.params['spinBox_2']['value']
        max_cond = self.params['spinBox']['value']

        # Error messages
        if n_var != n_prop:
            print 'Number of variables and number of variograms models are diferent.'
            return False

        if len(grid_var) == 0 or len(grid_krig) == 0:
            print 'Select the variables'
            return False

        SG_OK_list = []

        # Loop in every variable
        for i in xrange(0, n_var):

            # Getting variables
            prop_HD = props[i]

            prop_name = "Interpolated_" + str(prop_HD)

            var_str = ''

            indicator_group = "Indicator_group_" + str(i + 1)

            elipsoide = self.params['ellipsoidinput']['value']

            n_struct = int(
                self.params['indicator_regionalization_input'][indicator_group]['Covariance_input']['structures_count'])

            # Error message
            if n_struct == 0:
                print 'Variogram have no structures'
                return False

            # Loop in every variogram structure
            for j in xrange(0, n_struct):
                # Getting variogram parameters
                Structure = "Structure_" + str(j + 1)

                range1 = self.params['indicator_regionalization_input'][indicator_group]['Covariance_input'][Structure]['Two_point_model']['ranges']['range1']
                range2 = self.params['indicator_regionalization_input'][indicator_group]['Covariance_input'][Structure]['Two_point_model']['ranges']['range2']
                range3 = self.params['indicator_regionalization_input'][indicator_group]['Covariance_input'][Structure]['Two_point_model']['ranges']['range3']

                cont = self.params['indicator_regionalization_input'][indicator_group]['Covariance_input'][Structure]['Two_point_model']['contribution']

                cov_type = self.params['indicator_regionalization_input'][indicator_group]['Covariance_input'][Structure]['Two_point_model']['type']

                rake = self.params['indicator_regionalization_input'][indicator_group]['Covariance_input'][Structure]['Two_point_model']['angles']['rake']
                dip = self.params['indicator_regionalization_input'][indicator_group]['Covariance_input'][Structure]['Two_point_model']['angles']['dip']
                azimuth = self.params['indicator_regionalization_input'][indicator_group]['Covariance_input'][Structure]['Two_point_model']['angles']['azimuth']

                # Writing variogram parameters on a variable
                var_str = var_str + '<{} type="{}">  <Two_point_model  contribution="{}"  type="{}"   >      <ranges range1="{}"  range2="{}"  range3="{}"   />      <angles azimuth="{}"  dip="{}"  rake="{}"   />    </Two_point_model>    </{}> '.format(Structure, 'Covariance', cont, cov_type, range1, range2, range3, azimuth, dip, rake, Structure)

            # Calling ordinary kriging for each variable, using the variograms parameters above
            sgems.execute('RunGeostatAlgorithm  kriging::/GeostatParamUtils/XML::<parameters>  <algorithm name="kriging" />     <Variogram  structures_count="{}" >    {}  </Variogram>    <ouput_kriging_variance  value="0"  />     <output_n_samples_  value="0"  />     <output_average_distance  value="0"  />     <output_sum_weights  value="0"  />     <output_sum_positive_weights  value="0"  />     <output_lagrangian  value="0"  />     <Nb_processors  value="-2"  />    <Grid_Name value="{}" region=""  />     <Property_Name  value="{}" />     <Hard_Data  grid="{}"   property="{}"   region=""  />     <Kriging_Type  type="Ordinary Kriging (OK)" >    <parameters />  </Kriging_Type>    <do_block_kriging  value="1"  />     <npoints_x  value="5" />     <npoints_y  value="5" />     <npoints_z  value="5" />     <Min_Conditioning_Data  value="{}" />     <Max_Conditioning_Data  value="{}" />     <Search_Ellipsoid  value="{}" />    <AdvancedSearch  use_advanced_search="0"></AdvancedSearch>  </parameters>'.format(n_struct, var_str, grid_krig, prop_name, grid_var, prop_HD, min_cond, max_cond, elipsoide))

            SG_OK_list.append(sgems.get_property(grid_krig, prop_name))

        RT = (self.params['orderedpropertyselector']['value']).split(';')

        #Determinig geomodel
        GeoModel = []

        for i in range(len(SG_OK_list[0])):
            sgmin = 10e21
            for j in range(len(SG_OK_list)):
                if SG_OK_list[j][i] < sgmin:
                    sgmin = SG_OK_list[j][i]
                    k = j
            GeoModel.append(int(RT[k][-1]))

        sgems.set_property(grid_krig, 'Geologic_Model', GeoModel)

        #Operating softmax transformation
        if self.params['softmax_check']['value']=='1':
            gamma =float( self.params['Gamma']['value'])

        Prob_list = SG_OK_list[:]
        if self.params['softmax_check']['value']=='1':
            for i in range(len(SG_OK_list[0])):
                soma = 0
                for j in range(len(SG_OK_list)):
                    soma = soma + math.exp(-SG_OK_list[j][i]/gamma)
                for j in range(len(SG_OK_list)):
                    Prob_list[j][i] =math.exp(-SG_OK_list[j][i]/gamma)/soma

            for i in range(len(Prob_list)):
                sgems.set_property(grid_krig,'Probability_RT'+str(RT[i][-1]),Prob_list[i])

        '''#Operating servo-system
        if self.params['servo_check']['value'] == '1':
            var_rt = sgems.get_property(self.params['targe_prop']['grid'],self.params['targe_prop']['property'])
            RT_prop = proportion(var_rt)

            ran_path = random_path(len(Prob_list[0][0]))

            #Reading lambda1
            lambda1 = float(self.params['Lambda']['value'])

            # Determine proportions in geomodel
            proportions_geomodel = []
            for i in range(len(Prob_list)):
                soma = 0
                for j in ran_path:
                    if GeoModel[j] == RT[i]:
                        soma = soma + 1
                proportion_v = 0
                proportion_vector =[]
                for j in ran_path:
                    if GeoModel[j] == RT[i]:
                        proportion_v = (proportion_v + 1)/soma
                        proportions_vector.append(proportion_v)
                proportions_geomodel.append(proportions_vector)

            for i in range(len(Prob_list)):
                for j in ran_path:
                    Prob_list[i][j] = Prob_list[i][j] + lambda1*(RT_prop[i]-proprotions_geomodel[i][j])

            print Prob_list'''





        return True

    def finalize(self):

        return True

    def name(self):

        return "interpolator"

###############################################################
def get_plugins():
    return ["interpolator"]
