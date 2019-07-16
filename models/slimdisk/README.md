Relativistic slim disk model
============================

Relativistic slim disk model by Alexander Sadowski.

This model gives the structure of the polytropic relativistic slim disk model as described in Sadowski et al. (2011). 

The model interpolates through a set of tabulated numerical solutions. These solution are computed for a reference mass M=10 Msun, and although the model can rescale the disk structure to any mass, one should not go very much off the stellar mass range. Let's say masses 5-20 M_sun should give a reasonably good results.

## Installing

The model reads and interpolates tables with numerical slim disk solutions. Those tables are included in `data` subfolder in a compressed archive and have to be untared before using the model

```bash
cd models/slimdisk/data
tar -xvf slimdisk-polytropic.xz
```

## Options

The initialization function `diskmodel_init()` takes the following model-specific options through the _params_ parameter:

| Option  | Description |
|---------|-------------|
| path    | disk path to the directory with slim disk solutions (default current directory) |
| alpha   | value of alpha-viscosity paramater (range 0.005-0.1; default 0.1) |
| mdot    | mass accretion rate in units of [Mdot_Edd](https://github.com/mbursa/sim5/blob/public/doc/sim5lib-doc.md#sim5consth---constants-and-unit-conversions) (range 0.01-500; default 0.1) |
| lumi    | total disk luminosity in units of [L_Edd](https://github.com/mbursa/sim5/blob/public/doc/sim5lib-doc.md#sim5consth---constants-and-unit-conversions); takes precedence over _mdot_ option |

## References

* Sadowski 2009 ([ADS](https://adsabs.harvard.edu/abs/2009ApJS..183..171S)) - original formulation of the slim disk model
* Sadowski+2011 ([ADS](https://adsabs.harvard.edu/abs/2011A%26A...527A..17S)) - polytropic slim disk model

