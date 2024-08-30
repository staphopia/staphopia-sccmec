# Archived - Please use https://github.com/rpetit3/sccmec



*staphopia-sccmec* reads Staphopia results or assembly to assign a SCCmec type.

# staphopia-sccmec
[Staphopia](https://staphopia.emory.edu) includes a primer based SCCmec Typing scheme.
Unfortunately there was a disconnect between running the ananlysis and making the call.
The primers are aligned by BLAST in the [Staphopia Anlaysis Pipeline](https://github.com/staphopia/staphopia-ap) 
but the logic making the type and subtype calls was in the [Staphopia Website](https://github.com/staphopia/staphopia-web).

*staphopia-sccmec* has taken this process and merged into a single command. 

You can either provide results already produced by Staphopia or give an assembly
and a path to SCCmec data ([included in this repo](https://github.com/staphopia/staphopia-sccmec/tree/master/share/staphopia-sccmec/data)).


# Installation
### Bioconda
```
conda create -n staphopia-sccmec -c conda-forge -c bioconda staphopia-sccmec
```

### From Source
```
git clone git@github.com:staphopia/staphopia-sccmec.git
cd staphopia-sccmec
./bin/staphopia-sccmec.py
```

### Usage
```
usage: staphopia-sccmec [-h] [--assembly ASSEMBLY|ASSEMBLY_DIR|STAPHOPIA_DIR] [--staphopia STAPHOPIA_DIR] [--sccmec SCCMEC_DATA] [--ext STR] [--evalue EVALUE] [--hamming] [--json]
                        [--debug] [--depends] [--test] [--citation] [--version]

Determine SCCmec Type/SubType

optional arguments:
  -h, --help            show this help message and exit

Options:

  --assembly ASSEMBLY|ASSEMBLY_DIR|STAPHOPIA_DIR
                        Input assembly (FASTA format), directory of assemblies to predict SCCmec. (Cannot be used with --staphopia)
  --staphopia STAPHOPIA_DIR
                        Input directory of samples processed by Staphopia. (Cannot be used with --assembly)
  --sccmec SCCMEC_DATA  Directory where SCCmec reference data is stored (Default: /home/yaopeng/staphopia-sccmec/share/staphopia-sccmec/data).
  --ext STR             Extension used by assemblies. (Default: fna)
  --evalue EVALUE       evalue required by blast
  --hamming             Report the hamming distance of each type.
  --json                Report the output as JSON (Default: tab-delimited)
  --debug               Print debug related text.
  --depends             Verify dependencies are installed/found.
  --test                Run with example test data.
  --citation            Print citation information for using Staphopia SCCmec
  --version             show program's version number and exit
```


# Example Usage
*staphopia-sccmec* has two ways of being used, with Staphopia results or with an assembly.

### Typing Staphopia Results
Part of the Staphopia Anlaysis Pipeline is to BLAST the SCCmec primers against your assembly.
As mentioned previously, the logic to actually use those results and make a type call was
included on the database side.

*staphopia-sccmec* allows you to type using your existing Staphopia results. You just need to provide 
the path to a set of samples processed by staphopia and use the `--staphopia` parameter.


Below is an example of typing the Staphopia results for SRX085180.
```
staphopia-sccmec --staphopia /data/storage/semaphore/staphopia-v1/staphopia/ | head
sample  I       II      III     IV      V       VI      VII     VIII    IX      meca    Ia      IIa     IIb     IIIa    IVa     IVb     IVc     IVd     IVg     IVh
S.200218.00785  False   False   False   False   False   False   False   False   False   False   False   False   False   False   False   False   False   False   False   False
S.200218.00787  False   False   False   False   False   False   False   False   False   False   False   False   False   False   False   False   False   False   False   False
S.200218.00789  False   False   False   False   False   False   False   False   False   False   False   False   False   False   False   False   False   False   False   False
S.200218.00791  False   False   False   False   False   False   False   False   False   False   False   False   False   False   False   False   False   False   False   False
S.200218.00793  False   False   False   False   False   False   False   False   False   False   False   False   False   False   False   False   False   False   False   False
S.200218.00795  False   True    False   True    False   False   False   False   False   True    False   True    False   False   False   False   False   False   False   False
S.200218.00797  False   False   False   False   False   False   False   False   False   False   False   False   False   False   False   False   False   False   False   False
S.200218.00824  False   False   False   False   False   False   False   False   False   False   False   False   False   False   False   False   False   False   False   False
S.200218.00827  False   False   False   False   False   False   False   False   False   False   False   False   False   False   False   False   False   False   False   False
```

### Typing an Assembly
As an alternative to Staphopia results, you can provide an assembly, or a directory of assemblies. 
When an assembly is provided, a BLAST database is created of the assembly and the primers are blasted
against it.

Below is an example of typing the GCF_001580515 assembly, available in the `test/` directory.
#### Single Assembly
```
staphopia-sccmec --assembly test/GCF_001580515.1.fna
sample  I       II      III     IV      V       VI      VII     VIII    IX      meca    Ia      IIa     IIb     IIIa    IVa     IVb     IVc     IVd     IVg     IVh
GCF_001580515        False   False   False   True    False   False   False   False   False   True    False   False   False   False   True    False   False   False   False   False
```

#### Directory of Assemblies
```
ls ~/bactopia-dev/test-assemblies/
ERX140798.fna   ERX1711341.fna  ERX204841.fna   ERX3565118.fna  ERX385151.fna   ERX514151.fna  ERX956727.fna   SRX1885573.fna  SRX477083.fna   SRX6900463.fna
ERX1666543.fna  ERX1830501.fna  ERX2543896.fna  ERX3837876.fna  ERX3969884.fna  ERX770670.fna  SRX1885362.fna  SRX3883084.fna  SRX5659541.fna  SRX7775692.fna

staphopia-sccmec --assembly ~/bactopia-dev/test-assemblies/
sample  I       II      III     IV      V       VI      VII     VIII    IX      meca    Ia      IIa     IIb     IIIa    IVa     IVb     IVc     IVd     IVg     IVh
ERX204841       False   False   False   True    False   False   False   False   False   True    False   False   False   False   True    False   False   False   False   False
ERX2543896      False   False   False   False   False   False   False   False   False   True    False   False   False   False   False   False   False   False   False   False
ERX3565118      False   False   False   False   False   False   False   False   False   False   False   False   False   False   False   False   False   False   False   False
SRX3883084      False   False   True    False   True    False   True    False   False   True    False   False   False   True    False   False   False   False   False   False
SRX6900463      False   False   False   True    False   False   False   False   False   True    False   False   False   False   True    False   False   False   False   False
ERX956727       False   False   False   False   False   False   False   False   False   False   False   False   False   False   False   False   False   False   False   False
SRX1885362      False   False   False   False   False   False   False   False   False   False   False   False   False   False   False   False   False   False   False   False
SRX7775692      False   False   False   True    False   False   False   False   False   True    False   False   False   False   True    False   False   False   False   False
ERX514151       False   False   False   True    False   False   False   False   False   True    False   False   False   False   False   False   False   False   False   True
SRX477083       False   False   False   True    False   False   False   False   False   True    False   False   False   False   True    False   False   False   False   False
SRX5659541      False   True    False   False   False   False   False   False   False   True    False   True    False   False   False   False   False   False   False   False
ERX3969884      False   False   False   False   False   False   False   False   False   False   False   False   False   False   False   False   False   False   False   False
ERX1830501      False   False   False   False   False   False   False   False   False   False   False   False   False   False   False   False   False   False   False   False
ERX1711341      False   False   False   False   False   False   False   False   False   False   False   False   False   False   False   False   False   False   False   False
ERX1666543      False   False   False   False   False   False   False   False   False   False   False   False   False   False   False   False   False   False   False   False
ERX385151       False   False   False   False   False   False   False   False   False   False   False   False   False   False   False   False   False   False   False   False
ERX3837876      False   False   False   True    False   False   False   False   False   True    False   False   False   False   False   False   True    False   False   False
ERX140798       False   False   True    False   False   False   False   False   False   True    False   False   False   True    False   False   False   False   False   False
SRX1885573      False   False   False   False   False   False   False   False   False   False   False   False   False   False   False   False   False   False   False   False
ERX770670       False   False   False   True    False   False   False   False   False   True    False   False   False   False   False   False   False   False   False   True
```

### Alternate Outputs
#### JSON Output
You can also switch from a tab-delimited output to JSON, using the `--json` option.
```
staphopia-sccmec --assembly test/GCF_001580515.fna --json
[
    {
        "sample": "GCF_001580515",
        "I": false,
        "II": false,
        "III": false,
        "IV": true,
        "V": false,
        "VI": false,
        "VII": false,
        "VIII": false,
        "IX": false,
        "meca": true,
        "Ia": false,
        "IIa": false,
        "IIb": false,
        "IIIa": false,
        "IVa": true,
        "IVb": false,
        "IVc": false,
        "IVd": false,
        "IVg": false,
        "IVh": false
    }
]
```

#### Hamming Distance
By default, `staphopia-sccmec` reports `True` for exact primer matches and `False` for at least 1 base pair difference. As an alternative you can print the [Hamming Distance](https://en.wikipedia.org/wiki/Hamming_distance) instead. The Hamming Distance outputs the number of mismatches, with 0 being a perfect match.

```
staphopia-sccmec --assembly test/GCF_001580515.fna --hamming --json
[
    {
        "sample": "GCF_001580515",
        "I": 8,
        "II": 6,
        "III": 14,
        "IV": 0,
        "V": 16,
        "VI": 14,
        "VII": 16,
        "VIII": 20,
        "IX": 8,
        "meca": 0,
        "Ia": 19,
        "IIa": 20,
        "IIb": 14,
        "IIIa": 14,
        "IVa": 0,
        "IVb": 14,
        "IVc": 12,
        "IVd": 17,
        "IVg": 12,
        "IVh": 14
    }
]
```

# License
[MIT License](https://raw.githubusercontent.com/staphopia/staphopia-sccmec/master/LICENSE)

# Citation
Petit III RA, Read TD, [*Staphylococcus aureus viewed from the perspective of 40,000+ genomes.*](http://dx.doi.org/10.7717/peerj.5261) PeerJ 6, e5261 (2018).

# Author 

* Robert A. Petit III
* Twitter: [@rpetit3](https://twitter.com/rpetit3)

