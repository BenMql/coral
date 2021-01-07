#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Python routines for reading the input file TFP_parameters.in  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

def read_params_file():
   return open("CLR_parameters.in","r").readlines()

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
def read_rayleigh():
   list_of_parameters = read_params_file()
   line_of_interest = [x for x in list_of_parameters if 'Rayleig' in x]
   return float( line_of_interest[0].partition(" ")[0].replace('d','e')) 
   
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
def read_Ekman():
   list_of_parameters = read_params_file()
   line_of_interest = [x for x in list_of_parameters if 'Ekman' in x]
   return float( line_of_interest[0].partition(" ")[0].replace('d','e')) 
   
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
def read_Rossby():
   list_of_parameters = read_params_file()
   line_of_interest = [x for x in list_of_parameters if 'Rossby' in x]
   return float( line_of_interest[0].partition(" ")[0].replace('d','e')) 
   
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
def read_Lbox_X():
   list_of_parameters = read_params_file()
   line_of_interest = [x for x in list_of_parameters if 'x/lambda' in x]
   return float( line_of_interest[0].partition(" ")[0].replace('d','e')) 
   
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
def read_depth():
   list_of_parameters = read_params_file()
   line_of_interest = [x for x in list_of_parameters if 'depth/lambda' in x]
   return float( line_of_interest[0].partition(" ")[0].replace('d','e')) 
