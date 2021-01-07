#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Python routines for reading the input file TFP_parameters.in  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

def read_params_file():
   return open("CC_parameters.in","r").readlines()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
def read_NX():
   list_of_parameters = read_params_file()
   line_of_interest = [x for x in list_of_parameters if 'along X' in x]
   return int( line_of_interest[0].partition(" ")[0]) 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
def read_NZ():
   list_of_parameters = read_params_file()
   line_of_interest = [x for x in list_of_parameters if 'along Z' in x]
   return int( line_of_interest[0].partition(" ")[0]) 
   
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
def read_prandtl():
   list_of_parameters = read_params_file()
   line_of_interest = [x for x in list_of_parameters if 'Prandtl' in x]
   return float( line_of_interest[0].partition(" ")[0].replace('d','e')) 
   
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
def read_rayleigh():
   list_of_parameters = read_params_file()
   line_of_interest = [x for x in list_of_parameters if 'Rayleig' in x]
   return float( line_of_interest[0].partition(" ")[0].replace('d','e')) 
   
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
def read_Lbox_X():
   list_of_parameters = read_params_file()
   line_of_interest = [x for x in list_of_parameters if 'x/z aspect ratio' in x]
   return float( line_of_interest[0].partition(" ")[0].replace('d','e')) 
   
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
def read_H_over_ell():
   list_of_parameters = read_params_file()
   line_of_interest = [x for x in list_of_parameters if 'H_over_ell' in x]
   return float( line_of_interest[0].partition(" ")[0].replace('d','e')) 


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
def read_radiativeHeat_bool():
   list_of_parameters = read_params_file()
   line_of_interest = [x for x in list_of_parameters if 'Radiative' in x]
   string = line_of_interest[0]
   if (string[1:5].lower()=='true'):
      return True
   else:
      return False
