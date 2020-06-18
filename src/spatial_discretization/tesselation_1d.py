import numpy as np
from .utils import utils_1d as u1d

class Tesselation1D:

    @staticmethod
    def are_valid_opts(opts):
        return \
            opts["dimension"] == 1 and \
            opts["type"]      == "dg" and \
            "flux_type" in opts and \
            opts["basis"]     == {"type": "nodal_lgl", "order": 1} and \
            "grid" in opts


    def __init__(self, opts):
        if not Tesselation1D.are_valid_opts(opts):
            raise ValueError("Tesselation parameters are not supported")

        self.n_order = opts["basis"]["order"]
        #compute basic Legendre Gauss Lobato grid
        self.jgl = u1d.jacobi_gauss_lobatto(0,0,self.n_order)

        #build reference element matrices
        self.vander_matrix = u1d.vandermonde(self.n_order, self.jgl)
        self.diff_matrix   = u1d.differentiation_matrix(self.n_order, self.jgl, self.vander_matrix)

        #create surface integral terms
        self.lift = u1d.surface_integral_dg(self.n_order, self.vander_matrix)

        #build coordinates of all nodes
        box = opts["grid"]["box"]
        self.k_elem = int((box[1]-box[0])/opts["grid"]["steps"])
        [_,vx,_,etov] = u1d.mesh_generator(box[0],box[1],self.k_elem)
        self.vx = vx
        self.etov = etov
        self.nodes_coord = u1d.nodes_coordinates(self.n_order,self.etov,self.vx)

        #compute normals at element faces
        self.normals = u1d.normals(self.k_elem)
        
        #calculate geometric factors
        [self.rx,self.jacobian] = u1d.geometric_factors(self.nodes_coord,self.diff_matrix)
        [self.etoe, self.etof] = u1d.connect(self.etov)
        [self.vmap_m, self.vmap_p, self.vmap_b, self.map_b] =  u1d.build_maps( \
                                                                self.n_order, \
                                                                self.nodes_coord, \
                                                                self.etoe, \
                                                                self.etof)

        fmask_1 = np.where(np.abs(self.jgl+1)<1e-10)[0][0]
        fmask_2 = np.where(np.abs(self.jgl-1)<1e-10)[0][0]
        fmask = [fmask_1,fmask_2]
#        fx = self.nodes_coord[fmask]
        self.f_scale = 1/(self.jacobian[fmask])

        #initialize fields
        self.vars = {}
        for field in opts["fields"]:
            if field not in self.vars.keys:
                self.vars[field] = np.zeros((self.n_order+1),self.k_elem)

    def field(self, field_type):
        """
            field returns a one dimensional numpy array containing the 
            fields of type field_type
        """
        field = self.vars[field_type] 
        return 
    
    def curl(self, field_type):
        """
            curl returns a one dimensional numpy array containing the 
            discrete curl of fields of type field_type
        """
        field = self.vars[field_type]
        curl  = self.rx*np.matmul(self.diff_matrix,field)
        return curl

    def flux(self, field_type):
        """
            flux returns a one dimensional numpy array containing the 
            numerical flux of fields of type field_type
        """
        print("TBD") #TODO

    def get_smallest_distance(self):
        """
            get_smallest_distance returns ... 
        """
        print("TBD") #TODO
        
