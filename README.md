
# Accretion Disk Radial Structure Models

A collection of radial structure models of various accretion disk solutions, i.e. the radial dependence of their structural and thermodynamic quantities.

## Models

The following models are included:

* [Novikov-Thorne thin disk model](models/nt) (in models/nt)
* [Sadowski slim disk model](models/slimdisk) (in models/slimdisk)

See the model's README files for information about each model.

## Installation

To download and compile the model collection use

    git clone --recursive https://github.com/mbursa/disk-models.git
    cd disk-models
    make

The recursive option makes sure all submodules are included too (see bellow). 

Individual models can be compiled separately by running `make` in the directory of a particular model or by running
    
    make -C models/<model_name>
    
### Dependencies
Some of the models may have external dependencies, which are handled using [git submodules](https://git-scm.com/book/en/v2/Git-Tools-Submodules) that are located under `libs` subfolder. Submodules are downloaded and initialized automatically when `git clone --recurse-submodules` (for version >=2.13 of git) or `git clone --recursive`  (for git versions >1.6.5 and <2.13) is used to clone the repository. 

For already cloned repos, or older git versions (<1.6), use the following to initialize submodules:
```
git clone https://github.com/mbursa/disk-models.git
cd disk-models
git submodule update --init --recursive
``` 

Most models use [SIM5](https://github.com/mbursa/sim5) library.

## Usage

The compilation produces a linux shared library (`.so` extension) for each disk model. These libraries can be linked to programs at compile-time or they can be loaded by a program dynamically at run-time. 

The models can be used quite independently. However, the idea has been that the models are used along with [SIM5](https://github.com/mbursa/sim5) library to compute and compare the observed spectra and other characteristics.

### Dynamic linking

<!--
```
gcc -Ldisk-models/nt main.c -l:nt -o your_program`
```
-->
tbd

### Loading at runtime
Run-time loading in a program can be done using [`dlopen()`](http://man7.org/linux/man-pages/man3/dlopen.3.html) function. 
```C
void* lib_handle;
int (*diskmodel_init)(double M, double a, char* params);
// open library
void* lib_handle = dlopen(modellib, RTLD_NOW);
// assign library function to a local variable
*(void**)(&diskmodel_init) = dlsym(handle, "diskmodel_init");
// execute library function 
diskmodel_init(10.0, 0.0, "");
// close the library
dlclose(handle);
```

A complete example of how to dynamicly load and use a disk model library in a program is provided in [examples/runtime-loading](examples/runtime-loading).

### Loading with SIM5 library

[SIM5](https://github.com/mbursa/sim5) library contains a ready-to-use [interface for loading the disk models](https://github.com/mbursa/sim5/blob/public/doc/sim5lib-doc.md#sim5disk). It essentially uses  the approach described above and provides a C interface to access the functions of a disk model library that takes care of implementing the dynamic linking. A brief example how it can be used:
```C
#include "sim5lib.h"
// link the disk model libray disk-nt.so
diskmodel_init("./disk-nt.so", 10.0, 0.0, "mdot=0.1,alpha=0.1");
// call a funtion from the library (effective flux at raduis r = 10 r_g)
float F = diskmodel_flux(10.0);
// close the library
diskmodel_done();

```

### Interface
<!--
A complete description of the interface is given in a separate document.
-->
Each disk model implements the following set of functions that together define a common interface:

### `int diskmodel_init(char *modellib, double M, double a, char *params)`
Model initialization.  
Loads the shared library, links it into the program and calls its initialization routine.
**Parameters**  
* **modellib**: filesystem path to the `.so` library (string) 
* **M**: mass of BH [Msun] 
* **a**: spin of BH [0..1] 
* **params**: parametres that are passed to the library initialization function (string in key1=value1,key2=value2,... format)
**Return value**
Return the result of the library's initialization function or -1 of the library could not loaded or be initialized. 

### `void diskmodel_done()`
External disk model finitialization.  
Frees memory and unlinks the libraray. 

### `char* diskmodel_name()`
Model name.
Returns a pointer to a string with the model's name. 

### `void diskmodel_help()`

### `double diskmodel_r_min()`
Minimal radius of the disk (disk inner edge).
Gives minimal value for radius for which the functions provide valid results. E.g. for NT disk, this corresponds to the radius of the marginally stable orbit.
**Return value**
Radius of disk inner edge [GM/c2] 

### `double diskmodel_flux(double r)`
Local flux from one side of the disk.
rovides radial radiation flux dependence measured in local frame, i.e. flux measured by an observer that is at rest with respect to the fluid.
**Parameters**
* **R**: radius of emission [GM/c2]
**Return value**
Total outgoing flux from unit area on one side of the disk [erg cm-2 s-1]. 

### `double diskmodel_lumi()`
Total disk luminosity.
Luminosity is obtained by integrating local flux over the surface area of the disk (both sides) going into the whole sky (4pi solid angle). The integration makes a proper transformation of the flux from local to coordinate frame, but it ignores other relativistic effects, e.g. light bending.
```math
L = 2 * 2\pi \int F(r) (-U_t) r dr
```
**Return value**
Total disk luminosity of both surfaces [erg s-1] 

### `double diskmodel_mdot()`
Mass accretion rate.
Returns mass accretion rate in Eddington units of (Mdot_Edd*M).
**Return value**
Mass accretion rate in Eddington units. 

### `double diskmodel_sigma(double r)`
Column density.
Returns midplane column density of the fluid, i.e. the fluid density integrated from midplane to the disk surface, at a given radius.
**Parameters**
* **R**: radius (measured in equatorial plane) [rg]

**Return value**
Midplane column density in [g/cm2]. 


### `double diskmodel_l(double r)`
Specific angular momentum.
 double diskmodel_ell(double R)
Returns specific angular momentum of the fluid at given radius.
**Parameters**
* **R**: radius (measured in equatorial plane) [rg]
**Return value**
Specific angular momentum in [g.u.]. 

### `double diskmodel_vr(double r)`
Radial velocity.
 double diskmodel_vr(double R)
Returns bulk radial velocity of the fluid at given radius as measured by aan observer in the co-rotating frame.
**Parameters**
* **R**: radius (measured in equatorial plane) [rg]
**Return value**
Radial velocity in [speed_of_light]. 

### `double diskmodel_h(double r)`
Surface height.
 double diskmodel_h(double R)
Returns the scale-height of the surface of the disk above midplane at given radius.
**Parameters**
* **R**: radius (measured in equatorial plane) [rg]
**Return value**
Scale-height [rg]. 

### `double diskmodel_dhdr(double r)`
Derivative of surface height.
 double diskmodel_dhdr(double R)
Returns surface profile as derivative $`dH/dR`$ of its height above midplane at given radius.
**Parameters**
* **R**: radius (measured in equatorial plane) [rg]
**Return value**
Derivative of surface height. 

### `double diskmodel_eval(double r, int quantity)`
Other quantity evaluation.
 double diskmodel_eval(double R, int quantity)
Returns the value of a given quantity. The disk model may provide more quantities than the standard set. Additional quantities may be accessed using this function by providing the quantity identifier.
**Parameters**
* **R**: radius (measured in equatorial plane) [rg] 
* **quantity**: quantity identified code
**Return value**
The value of the requested quantity. 

### `void diskmodel_params(FILE* output)`
Prints model parameters.
 void diskmodel_params(FILE *stream)/! if(stream) _diskmodel_params(stream)
Writes down the parameters of the model to the given stream. **Parameters**
* **stream**: stream to write to 

### `void diskmodel_dump()`
Prints the disk structure as a function of radius.
 void diskmodel_dump(char *filename)
The function prints the profile of all quantities as a function of radius from r_ms to some outer radius (~2000 rg). It prints to a file identified by its path (overwrites existing) and if that is empty it prints to STDOUT.
**Parameters**
* **filename**: Path to a file that should be written. If NULL then it prints to STDOUT.

## Citing

If you have used this software in your work, please acknowledge it by citing its ASCL record:  
[ascl:XXXX.XXX](XXXX.XXX). See their [citing guidelines](https://ascl.net/wordpress/about-ascl/citing-ascl-code-entries/).


## License

Disk models are released under the MIT License.

