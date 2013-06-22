import os

def run_reconstruction_periodic_box(data, rand, zflag):
  if(zflag == 0):
    os.system("./rpb -d "+data+" -r "+rand)
  elif(zflag == 1):
    os.system("./rpb -z -d "+data+" -r "+rand)

def run_reconstruction_nonperiodic_box(data, rand, counts):
  os.system("./rnpb -d "+data+" -r "+rand+" -c "+counts)

def run_pk_periodic_box(data):
  os.system("./ppb -d "+data)

def run_pk_nonperiodic_box(data, rand):
  os.system("./pnpb -d "+data+" -r "+rand)

def run_xi_periodic_box(data, rand):
  os.system("./cpb -d "+data+" -r "+rand)

def run_xi_nonperiodic_box(data, rand):
  os.system("./cnpb -d "+data+" -r "+rand)

def run_xi_rec_periodic_box(data, rand, srand):
  os.system("./crpb -d "+data+" -r "+rand+" -s "+srand)

def run_xi_rec_nonperiodic_box(data, rand, srand):
  os.system("./crnpb -d "+data+" -r "+rand+" -s "+srand)

################################RECONSTRUCTION#################################

#===============================Periodic Box===================================
#USE:
#run_reconstruction_periodic_box([data_file_name],[random_file_name],zflag)
#zflag = 0 => real space, = 1 => redshift space

#INPUTS:
#-galaxies
#   [data_file_name].txt
#   Format: x, y, z
#-randoms
#   [random_file_name].txt
#   Format: x, y, z

#OUTPUTS:
#-reconstructed galaxies
#   [data_file_name]_Lbox[size of box]_Ngrid[size of grid].txt
#   Format: x, y ,z
#-reconstructed randoms
#   [data_file_name]_Lbox[size of box]_Ngrid[size of grid]_rand.txt
#   Format: x, y ,z

#RUN:
#real space
run_reconstruction_periodic_box("per_box","per_box_rand",0)
#redshift space
run_reconstruction_periodic_box("per_box_z","per_box_rand",1)
#==============================================================================

#============================Non-Periodic Box==================================
#USE:
#run_reconstruction_nonperiodic_box([data_file_name],[random_file_name],
#                                   [random_cell_counts_file_name])

#INPUT:
#-galaxies
#   [data_file_name].txt
#   Format: x, y ,z, weight
#-randoms
#   [random_file_name].txt
#   Format: x, y ,z, weight
#-random cell counts
#   [random_cell_counts_file_name].txt
#   Format: counts, ra, dec
#   ***This file is produced by randomly placing [rand_per_cell] randoms in 
#   each cell and counting how many fall within the survey volume. The variable
#   [rand_per_cell] needs to be updated in the code. The number of cells used
#   to generate this file must equal Ngrid and the Box-Size must match Lbox.
#   Also, the mins[] array may need to be adjusted depending on how the 'counts'
#   file is generated. mins[] is basically used to shift the particles such that
#   the volume they cover spatially matches the counts file.***

#OUTPUT:
#-reconstructed galaxies
#   [data_file_name]_Lbox[size of box]_Ngrid[size of grid].txt
#   Format: x, y ,z, weight
#-reconstructed randoms
#   [random_file_name]_Lbox[size of box]_Ngrid[size of grid].txt
#   Format: x, y ,z, weight

#RUN:
run_reconstruction_nonperiodic_box("nonper_box","nonper_box_rand",\
                                   "nonper_box_random_cell_counts")
#==============================================================================
###############################################################################

################################POWER SPECTRUM#################################
#===============================Periodic Box===================================
#USE:
#run_pk_periodic_box([data_file_name])

#INPUT:
#-galaxies
#   [data_file_name].txt
#   Format: x, y, z

#OUTPUT:
#-galaxy power spectrum
#   [data_file_name].pk
#   Format: k, P(k)

#RUN:
run_pk_periodic_box("per_box")
#==============================================================================

#============================Non-Periodic Box==================================
#USE:
#run_pk_nonperiodic_box([data_file_name],[random_file_name])

#INPUT:
#INPUTS MUST CONTAIN n(z)!!
#-galaxies
#   [data_file_name].txt 
#   Format: x, y, z, weight, n(z)
#-randoms
#   [random_file_name].txt
#   Format: x, y, z, weight, n(z)

#OUTPUT:
#-galaxy power spectrum
#   [data_file_name].pk
#   Format: k, P(k)

#RUN:
run_pk_nonperiodic_box("nonper_box","nonper_box_rand")
#==============================================================================
###############################################################################

##############################CORRELATION FUNCTION#############################
#=======================Periodic Box, Before Reconstruction====================
#USE:
#run_xi_periodic_box([data_file_name],[random_file_name])

#INPUT:
#-galaxies
#   [data_file_name].txt
#   Format: x, y, z
#-randoms
#   [random_file_name].txt
#   Format: x, y, z

#OUTPUT:
#-galaxy correlation function
#   [data_file_name].xi
#   Format: r, xi0(r), xi2(r), DD(r), DR(r), RR(r) [DD, DR & RR are normalized]

#RUN:
run_xi_periodic_box("per_box","per_box_rand")
#==============================================================================

#=====================Non-Periodic Box, Before Reconstruction==================
#USE:
#run_xi_nonperiodic_box([data_file_name],[random_file_name])

#INPUT:
#-galaxies
#   [data_file_name].txt
#   Format: x, y, z, weight
#-randoms
#   [random_file_name].txt
#   Format: x, y, z, weight

#OUTPUT:
#-galaxy correlation function
#   [data_file_name].xi
#   Format: r, xi0(r), xi2(r), DD(r), DR(r), RR(r) [DD, DR & RR are normalized]

#RUN:
run_xi_nonperiodic_box("nonper_box","nonper_box_rand")
#==============================================================================

#=======================Periodic Box, After Reconstruction=====================
#USE:
#run_xi_rec_periodic_box([reconstructed_data_file_name],[random_file_name],
#                        [reconstructed_random_file_name])

#INPUT:
#-reconstructed galaxies (this can be the output from 
#                         reconstruction_periodic_box.cpp)
#   [reconstructed_data_file_name].txt
#   Format: x, y, z
#-randoms
#   [random_file_name].txt
#   Format: x, y, z
#-reconsturcted randoms (this can be the output from 
#                        reconstruction_periodic_box.cpp)
#   [reconstructed_random_file_name]
#   Format: x, y, z

#OUTPUT:
#-reconstructed galaxy correlation function
#   [reconstructed_data_file_name].xi
#   Format: r, xi0(r), xi2(r), DD(r), DS(r), SS(r) [DD, DS & SS are normalized]

#RUN:
run_xi_rec_periodic_box("per_box_rec_Lbox1500_Ngrid512",\
                        "per_box_rand",\
                        "per_box_rec_Lbox1500_Ngrid512_rand")
#==============================================================================

#=====================Non-Periodic Box, After Reconstruction===================
#USE:
#run_xi_rec_nonperiodic_box([reconstructed_data_file_name],[random_file_name],
#                           [reconstructed_random_file_name])

#INPUT:
#-reconstructed galaxies (this can be the output from 
#                         reconstruction_non-periodic_box.cpp)
#   [reconstructed_data_file_name].txt
#   Format: x, y, z, weight
#-randoms
#   [random_file_name].txt
#   Format: x, y, z, weight
#-reconsturcted randoms (this can be the output from 
#                        reconstruction_non-periodic_box.cpp)
#   [reconstructed_random_file_name]
#   Format: x, y, z, weight

#OUTPUT:
#-reconstructed galaxy correlation function
#   [reconstructed_data_file_name].xi
#   Format: r, xi0(r), xi2(r), DD(r), DS(r), SS(r) [DD, DS & SS are normalized]

#RUN:
run_xi_rec_nonperiodic_box("nonper_box_rec_Lbox3400_Ngrid256",\
                           "nonper_box_rand",\
                           "nonper_box_rec_Lbox3400_Ngrid256_rand")
#==============================================================================
###############################################################################
