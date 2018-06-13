"""Jordan-Wigner transformation or when the doc string is longer than the actual code"""
#author tjdim, Cambridge MA 2018
from scipy.sparse import identity
import _physical_constants as pc
import numpy as np



def jordan_wigner_trafo(nsites):
    """Performs Jordan-Wigner transformation to built on-site creation and annihilation operators

        Args
        ----

            nsites (int): Number of sites in the real-space lattice
 
 
        Returns
        -------

            a_op_dagger (dict) : keys (int): index of site
                                 values (2^(nsites)-dim array): on-site creation operator a_j_dagger 
            
            a_op (dict) : keys (int): index of site j
                          values (2^(nsites)-dim array): on-site annihilation operators a_j


        References
        ----------

        .. [1] R.Babbush et al., NJP, 18 033032, 2016. 

        Examples
        --------

            >>> from qulearn.transformations import _jordan_wigner 
 
    """
    IDENTITY_2D = identity(2).toarray()    
    #initialize operators 
    creation_operator = {}
    annihilation_operator = {}
    #built on-site operator for each site_j
    if nsites > 6:
        # Exponential wall of many-body problem
        print "Warning: Too many sites! \n Try instead with fewer sites. "
    else:
        for site_j in range(1, nsites + 1):        
            #initialize  
            a_op_dagger = 1.0
            #loop over all sites to built on-site operator a_dagger_j
            for site_s in range(1, site_j):  
                a_op_dagger = np.kron(pc.SIGMA_Z, a_op_dagger)          
            #built kronecker product with sigma_plus
            a_op_dagger = np.kron(pc.SIGMA_PLUS, a_op_dagger)
            #for destruction operator: a_op[site_j] = np.kron(pc.SIGMA_MINUS, kron_product)
            #fill remaining tensor product with identity matrix to match dimension
            for site_s in range(site_j, nsites):
                a_op_dagger = np.kron(IDENTITY_2D, a_op_dagger)  
            creation_operator[site_j] = a_op_dagger
            #a^dagger is the transpose of the operator a
            annihilation_operator[site_j] = np.transpose(a_op_dagger)  
        return creation_operator, annihilation_operator

 
ns = input('Enter number of sites: ')

try:
    create, destroy = jordan_wigner_trafo(ns)
except:
    while ns > 6:
        ns = input('Enter number of sites: ')
        create, destroy = jordan_wigner_trafo(ns)
#

print np.matrix(create)
print('------------')
print np.matrix(destroy[1])

