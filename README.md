# k4-project-template


This repository can be a starting point and template for projects using the Key4HEP software stack.


## Dependencies

* ROOT

* PODIO

* Gaudi

* k4FWCore

## Installation


``` bash
source /cvmfs/sw.hsf.org/key4hep/setup.sh
mkdir build install
cd build;
cmake .. -DCMAKE_INSTALL_PREFIX=../install
make install
```

Alternatively you can source the nightlies instead of the releases:

``` bash
source /cvmfs/sw-nightlies.hsf.org/key4hep/setup.sh
```

Note that if you source the releases and use the current version of this
repository this is not guaranteed to work as there could be changes since this
repository was built for the release. What you can do in this case is to
checkout a previous tag, for example:

``` bash
git checkout v0.3.0
```

This is because the releases are only built with tagged versions of the
packages. With the nightlies this repository should always work; if it doesn't
please [open an
issue](https://github.com/key4hep/k4-project-template/issues/new/choose).

## Execute Examples

Make sure that `../install/lib` and `../install/python` are in `LD_LIBRARY_PATH`
and `PYTHONPATH` respectively, by doing:
``` bash
export LD_LIBRARY_PATH=$PWD/../install/lib:$LD_LIBRARY_PATH
export PYTHONPATH=$PWD/../install/python:$PYTHONPATH
```
and then run the examples like this:

``` bash
k4run ../k4ProjectTemplate/options/createHelloWorld.py
k4run ../k4ProjectTemplate/options/createExampleEventData.py
```


## References:
These could perhaps be usefule for newcomers.
1. [lhcb-98-064 COMP](https://cds.cern.ch/record/691746/files/lhcb-98-064.pdf)
2. [Hello World in the Gaudi Framework](https://lhcb.github.io/DevelopKit/02a-gaudi-helloworld)
