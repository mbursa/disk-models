## Disk model interface

Each disk model implements a set of functions that defines the common interface. The following list of functions specifies the quantities that each model provides along with a specification of the function parameters and resulting values.

**Note 1**: Whereever a function has radius _R_ as a parameter, it is assumed to be the cylindrical radius, i.e. the distance to the point measured parallel to the equatorial plane.

**Note 2**: If possible, the function takes arguments and return results in dimensionless units. E.g. radius _R_ in units of gravitational radius `r_g = GM/c^2`, velocity in units of the speed of light, etc.

#### `int diskmodel_init(char *modellib, double M, double a, char *params)`
Model initialization.  
Loads the shared library, links it into the program and calls its initialization routine passing the paramters of the model.

**Parameters**:  
* **modellib**: filesystem path to the `.so` library (string) 
* **M**: mass of the central black hole [Msun] 
* **a**: spin of the central black hole [0..1] 
* **params**: model-specific parametres that are passed to the library initialization function (string in key1=value1,key2=value2,... format)

**Return value**:  
Returns the result of the library's initialization function (typically 0 on success) or -1 of the library could not loaded or be initialized. 

#### `void diskmodel_done()`
Model finitialization.  
Closes the model, frees memory and unlinks the libraray. 

#### `char* diskmodel_name()`
Model name.

**Return value**:  
Returns a pointer to a string with the model's name. Note: the pointer references a library's internal string buffer and it is  memory-managed by the library.

<!--
#### `void diskmodel_help()`
Model name.
Returns a pointer to a string with the model's name. 
-->

#### `double diskmodel_r_min()`
Minimal radius of the disk (disk inner edge).
Gives minimal value for radius for which the functions provide valid results. E.g. for Novikov-Thorne disk, this corresponds to the radius of the marginally stable orbit (ISCO).

**Return value**:  
Radius of the disk inner edge [GM/c2] 

#### `double diskmodel_flux(double R)`
Local flux from one side of the disk.
Gives the local radiative flux from the surface of the disk at equatorial radius _R_, i.e. flux measured by an observer that is at rest with respect to the orbiting fluid.

**Parameters**: 
* **R**: radius of emission [r_g]

**Return value**:  
Total outgoing flux from unit area on one side of the disk [erg cm-2 s-1]. 

#### `double diskmodel_lumi()`
Total disk luminosity.
Luminosity is obtained by integrating local flux over the surface area of the disk (both sides) going into the whole sky (4pi solid angle). The integration makes a proper transformation of the flux from local to coordinate frame, but it ignores other relativistic effects, e.g. light bending.
```math
L = 2 * 2\pi \int F(r) (-U_t) r dr
```
**Return value**:  
Total disk luminosity from both surfaces [erg s-1] 

#### `double diskmodel_mdot()`
Mass accretion rate.
Returns mass accretion rate in Eddington units of (Mdot_Edd*M). See [SIM5](https://github.com/mbursa/sim5) for the definition of Eddington mass accretion rate _Mdot_Edd_.

**Return value**:  
Mass accretion rate in Eddington units. 

#### `double diskmodel_sigma(double R)`
Column mass density.
Returns a column mass density of the fluid, i.e. the fluid density vertically integrated from midplane to the disk surface, at a given radius.

**Parameters**:  
* **R**: radius (measured in equatorial plane) [r_g]

**Return value**:  
Midplane column density in [g/cm2]. 


#### `double diskmodel_l(double R)`
Specific angular momentum.
Returns specific angular momentum of the fluid at a given radius.

**Parameters**:  
* **R**: radius (measured in equatorial plane) [r_g]
* 
**Return value**:  
Specific angular momentum [in geometrical units]. 

#### `double diskmodel_vr(double R)`
Radial velocity.
Returns bulk radial velocity of the fluid at given radius as measured by an observer in the co-rotating frame.

**Parameters**:  
* **R**: radius (measured in equatorial plane) [r_g]

**Return value**:  
Radial velocity in [speed_of_light]. Positive value means outflow, negative value means inflow.

#### `double diskmodel_h(double R)`
Surface height.
Returns the scale-height of the surface of the disk above midplane at given radius.

**Parameters**:  
* **R**: radius (measured in equatorial plane) [r_g]

**Return value**:  
Scale-height [r_g]. Note: This should represent the location of the photosphere, where the effective optical depth tau ~= 1.

#### `double diskmodel_dhdr(double R)`
Derivative of surface height.
Returns surface profile as a derivative $dH/dR$ of its height above midplane at given radius.

**Parameters**:  
* **R**: radius (measured in equatorial plane) [r_g]

**Return value**:  
Derivative of surface height. 

#### `double diskmodel_eval(double r, int quantity)`
Other quantity evaluation.
Returns the value of a given quantity. The disk model may provide more quantities than the standard set. Additional quantities may be accessed using this function by providing the quantity identifier.

**Parameters**:  
* **R**: radius (measured in equatorial plane) [r_g] 
* **quantity**: quantity identified code

**Return value**:  
The value of the requested quantity. 

#### `void diskmodel_params(FILE* output)`
Prints model parameters.
Writes down the parameters of the model to the given stream.

 **Parameters**:  
* **stream**: stream to write to 

#### `void diskmodel_dump(FILE* filename)`
Prints the disk structure as a function of radius.
The function prints the profile of all quantities as a function of radius from r_ms to some outer radius (~2000 rg). It prints to a file identified by its path (overwrites existing) and if that is empty it prints to STDOUT.

**Parameters**:  
* **filename**: Path to a file that should be written. If NULL then it prints to STDOUT.



