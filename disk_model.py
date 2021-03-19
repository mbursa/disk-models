# Python interface to disk model libraries.
#
# Defines a class to load a shared object libraray (.so) 
# with a disk model and to call its functions


import sys
import ctypes
import importlib
import logging


class DiskModel:
    """
    Disk model that links an external library/module.
    """

    def __init__(self, model_lib, bh_mass, bh_spin, options=''):
        """
        Initializes an external disk model.

        Args:
            model_lib: a shared library or a python module that provides disk model methods
            bh_mass: black hole mass parameter [M_sun]
            bh_spin: black hole spin parameter [0..1]
            options: other options tht are passed to the model library (a sequence of key=value pairs delimited by comma)
        """
        self.lib = None

        try:
            if (model_lib and model_lib.endswith('.so')):
                # load shared library and setup c_types interface for its functions
                # example: http://gestaltrevision.be/wiki/python/cfun
                self.lib = ctypes.cdll.LoadLibrary(model_lib)

                self.lib.diskmodel_init.argtypes = [ctypes.c_double, ctypes.c_double, ctypes.c_char_p]
                self.lib.diskmodel_init.restype = None

                self.lib.diskmodel_done.argtypes = []
                self.lib.diskmodel_done.restype = None

                self.lib.diskmodel_name.argtypes = []
                self.lib.diskmodel_name.restype = ctypes.c_char_p

                self.lib.diskmodel_r_min.argtypes = []
                self.lib.diskmodel_r_min.restype = ctypes.c_double

                self.lib.diskmodel_mdot.argtypes = []
                self.lib.diskmodel_mdot.restype = ctypes.c_double

                self.lib.diskmodel_lumi.argtypes = []
                self.lib.diskmodel_lumi.restype = ctypes.c_double

                self.lib.diskmodel_flux.argtypes = [ctypes.c_double]
                self.lib.diskmodel_flux.restype = ctypes.c_double

                self.lib.diskmodel_sigma.argtypes = [ctypes.c_double]
                self.lib.diskmodel_sigma.restype = ctypes.c_double

                self.lib.diskmodel_l.argtypes = [ctypes.c_double]
                self.lib.diskmodel_l.restype = ctypes.c_double

                self.lib.diskmodel_vr.argtypes = [ctypes.c_double]
                self.lib.diskmodel_vr.restype = ctypes.c_double

                self.lib.diskmodel_h.argtypes = [ctypes.c_double]
                self.lib.diskmodel_h.restype = ctypes.c_double

                self.lib.diskmodel_dhdr.argtypes = [ctypes.c_double]
                self.lib.diskmodel_dhdr.restype = ctypes.c_double

                self.lib.diskmodel_dump.argtypes = []
                self.lib.diskmodel_dump.restype = None

            elif (model_lib and model_lib.endswith('.py')):
                raise Exception('"'+model_lib+'"'+' is not an .so file')
            #end if
        except AttributeError as e:
            print('Missing required function in disk module.\n', e.args)
            raise

        # init model
        if (self.lib):
            self.lib.diskmodel_init(bh_mass, bh_spin, str.encode(options))
            self.r_min = self.lib.diskmodel_r_min()
            self.name  = self.lib.diskmodel_name().decode()
            self.mdot  = self.lib.diskmodel_mdot()
            self.lumi  = self.lib.diskmodel_lumi()
            logging.info("# Disk model: '%s' M=%.1f a=%.3f mdot=%.3e (%.5e) lum=%.3e", self.name, bh_mass, bh_spin, self.mdot, self.mdot*bh_mass*2.225475942e+18, self.lumi)
        #end if

    #end of def


    def __del__(self):
        if (self.lib and self.lib.diskmodel_done): self.lib.diskmodel_done()
    #end of def

    def close(self):
        if (self.lib and self.lib.diskmodel_done): self.lib.diskmodel_done()
        self.lib = None
    #end of def

    def flux(self, R): return self.lib.diskmodel_flux(R) if self.lib else 0.0

    def sigma(self, R): return self.lib.diskmodel_sigma(R) if self.lib else 0.0

    def l(self, R): return self.lib.diskmodel_l(R) if self.lib else 0.0

    def vr(self, R): return self.lib.diskmodel_vr(R) if self.lib else 0.0

    def h(self, R): return self.lib.diskmodel_h(R) if self.lib else 0.0

    def dhdr(self, R): return self.lib.diskmodel_dhdr(R) if self.lib else 0.0

    def eval(self, R, p): return self.lib.diskmodel_eval(R,p) if self.lib else None

    def dump(self): return self.lib.diskmodel_dump() if self.lib else None

#end of class





# function to test the class
def test():
    if (len(sys.argv) < 4):
        print("Usage: %s sd|nt spin mdot")
        sys.exit()
    #end if

    model = sys.argv[1]    
    spin  = float(sys.argv[2])
    mdot  = float(sys.argv[3])

    m = DiskModel('models/nt/disk-'+model+'.so', 10.0, spin, 'mdot=%f'%(mdot))

    print('# model dump:', m.name)
    print('# rmin:', m.r_min)
    print('# mdot:', m.mdot)
    print('# lumi:', m.lumi)
    print('# [radius  flux  sigma  ell  vr  h  dhdr]')
    r = m.r_min
    while (r < 1e4):
        print(r, m.flux(r), m.sigma(r), m.l(r), m.vr(r), m.h(r), m.dhdr(r))
        r *= 1.05
    #end while
#end def


# if executed as main file then
# run the main function and exit with returned error code
if __name__ == "__main__": sys.exit(test())


