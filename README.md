## RFsPy: Receiver Function Processing and Plotting codes. 

Pascal Audet and Andrew Schaeffer

August 2019

### Installation

### Dependencies

You will need **Python 3.7+**.
Also, the following packages are required:

- [`obspy`](https://github.com/obspy/obspy/wiki)

Note that `numpy` and `matplotlib` are installed as dependencies of `obspy`. See below for full installation details. 

### Conda environment

You can create a custom [conda environment](https://conda.io/docs/user-guide/tasks/manage-environments.html)
where `RFsPy` can be installed along with its dependencies.

Clone the repository:
```bash
git clone https://gitlab.com/uottawa-geophysics/SeismoPy/RFsPy.git
cd RFsPy
```

Create the environment from the `rfpy_env.yml` file:
```bash
conda env create -f rfpy_env.yml
```
Activate the newly created environment:
```bash
conda activate rfpy
```

### 1) Installing using pip

Once the previous steps are performed, you can install `RFsPy` using pip:
```bash
pip install .
```

Installation after local or remote (gitlab/github) code updates:

```bash
pip install --upgrade .
```

### 2) Building and Installing

Alternatively, you can build and install the project (from the root of the source tree, e.g., inside the cloned `RFsPy` directory):

```bash
python setup.py build 
python setup.py install
```


Please note, if you are actively working on the code, or making frequent edits, it is advisable
to perform the pip installation with the -e flag. This enables an editable installation, where
symbolic links are used rather than straight copies. This means that any changes made in the
local folders will be reflected in the packages available on the system.


## Package Contents

<!-- ### StDb Module:

*  Classes.py -> contains the Class definitions
*  convert.py -> subroutines for converting to/from csv types
*  io.py -> input and output routines for loading and writing the station databases.

### Scripts: 
Python scripts making use of the module for manipulation and creation of databases

* ls_stdb.py -> script for easily viewing the contents of a database
* gen_stdb.py -> script to create a databases from a text file (specific csv format)
* query_stdb.py -> script to create a new database based on a query to a network client
* dump_stdb.py -> script to export a database into csv format
* edit_stdb.py -> interactive entry-by-entry editing of a database file (not recommended if editing large numbers of entries)
* append_stdb.py -> add new entries to an existing database
* merge_stdb.py -> combine multiple station databases together
* update_stdb.py -> convert an old format station database to the new format (v1 to v2)
 -->

### Examples



## Components

Note that presently (1.0.1) this code uses an outdated version of the station database that is hard-coded within the io_utils component of the rf module. This should be upgraded to use the Standalone StDb module in the future.

