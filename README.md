*staphopia-sccmec* reads Staphopia results or assembly to assign a SCCmec type.

# staphopia-sccmec
[Staphopia](https://staphopia.emory.edu) includes a primer based SCCmec Typing scheme.
Unfortunately there was a disconnect between running the ananlysis and making the call.
The primers are aligned by BLAST in the [Staphopia Anlaysis Pipeline](https://github.com/staphopia/staphopia-ap) 
but the logic making the type and subtype calls was in the [Staphopia Website](https://github.com/staphopia/staphopia-web).

*staphopia-sccmec* has taken this process and merged into a single command. 

You can either provide results already produced by Staphopia or give an assembly
and a path to SCCmec data ([included in this repo](https://github.com/staphopia/staphopia-sccmec/tree/master/data)).


# Installation
### Bioconda
*staphopia-sccmec* is currently not available from Bioconda, but it is being worked on.

### From Source
```
git clone git@github.com:staphopia/staphopia-sccmec.git
cd staphopia-sccmec.git
./staphopia-sccmec.py
```

### Usage
```
usage: staphopia-sccmec.py [-h] [--prefix STR] [--staphopia] [--hamming]
                           [--cpus INT] [--keep_files] [--debug] [--quiet]
                           [--depends] [--version]
                           ASSEMBLY|PRIMERS_JSON SCCMEC_DATA|SUBTYPES_JSON
                           OUTPUT_DIR

Determine SCCmec Type/SubType

optional arguments:
  -h, --help            show this help message and exit

Options:

  ASSEMBLY|PRIMERS_JSON
                        Input assembly (FASTA format) to predict SCCmec. Or,
                        "primers.json" from Staphopia, this requires "--
                        staphopia"
  SCCMEC_DATA|SUBTYPES_JSON
                        Directory where SCCmec reference data is stored. Or,
                        "subtypes.json" from Staphopia, this requires "--
                        staphopia"
  OUTPUT_DIR            Directory to output results to
  --prefix STR          Prefix (e.g. sample name) to use for outputs (Default:
                        assembly file without extension)
  --staphopia           Inputs are results from Staphopia. The first arguement
                        (ASSEMBLY) should be "primers.json" and the second
                        arguement (SCCMEC_DATA) should be "subtypes.json".
                        These files are found in "analyses/sccmec" folder of
                        Staphopia results.
  --hamming             Report the hamming distance of each type.
  --cpus INT            Number of processors to use.
  --keep_files          Keep all output files (Default: remove blastdb.
  --debug               Print debug related text.
  --quiet               Only critical errors will be printed.
  --depends             Verify dependencies are installed/found.
  --version             show program's version number and exit

```


# Example Usage
*staphopia-sccmec* has two ways of being used, with Staphopia results or with an assembly.

### Typing Staphopia Results
Part of the Staphopia Anlaysis Pipeline is to BLAST the SCCmec primers against your assembly.
As mentioned previously, the logic to actually use those results and make a type call was
included on the database side.

*staphopia-sccmec* allows you to type using your existing Staphopia results. You just need to provide 
the `primers.json` and `subtypes.json` outputs from Staphopia and the `--staphopia` parameter. These 
files are located in the `analyses/sccmec` folder.


Below is an example of typing the Staphopia results for SRX085180.
```
./staphopia-sccmec.py  SRX085180/analyses/sccmec/primers.json SRX085180/analyses/sccmec/subtypes.json ./ --staphopia --prefix SRX085180
2020-02-21 23:53:33:root:INFO - Processing Staphopia outputs
2020-02-21 23:53:33:root:INFO - Writing predicted SCCmec type based on primers to .//SRX085180/sccmec-primer-type.json
{
    "cassette": {
        "sample": "SRX085180",
        "I": false,
        "II": false,
        "III": false,
        "IV": true,
        "V": false,
        "VI": false,
        "VII": false,
        "VIII": false,
        "IX": false,
        "meca": true
    },
    "subtype": {
        "sample": "SRX085180",
        "Ia": false,
        "IIa": false,
        "IIb": false,
        "IIIa": false,
        "IVa": false,
        "IVb": false,
        "IVc": false,
        "IVd": true,
        "IVg": false,
        "IVh": false
    }
}
```

Based on the results it appears that SRX085180 is SCCmec Type IVd.

### Typing an Assembly
As an alternative to Staphopia results, you can provide an assembly. When an assembly is 
provided, a BLAST database is create of the assembly and the primers are blasted against
it.

Below is an example of typing the SRX085180 based on its assembly.
```
./staphopia-sccmec.py  /data/storage/servers/merlin-home/rpetit/cgc-staphopia/staphopia/SRX085/SRX085180/SRX085180/analyses/assembly/SRX085180.contigs.fasta.gz data/ ./
2020-02-21 23:57:49:root:INFO - Make BLAST database for /data/storage/servers/merlin-home/rpetit/cgc-staphopia/staphopia/SRX085/SRX085180/SRX085180/analyses/assembly/SRX085180.contigs.fasta.gz
2020-02-21 23:57:49:root:INFO - BLAST SCCmec primers against /data/storage/servers/merlin-home/rpetit/cgc-staphopia/staphopia/SRX085/SRX085180/SRX085180/analyses/assembly/SRX085180.contigs.fasta.gz
2020-02-21 23:57:50:root:INFO - Writing predicted SCCmec type based on primers to .//SRX085180.contigs.fasta/sccmec-primer-type.json
{
    "cassette": {
        "sample": "SRX085180.contigs.fasta",
        "I": false,
        "II": false,
        "III": false,
        "IV": true,
        "V": false,
        "VI": false,
        "VII": false,
        "VIII": false,
        "IX": false,
        "meca": true
    },
    "subtype": {
        "sample": "SRX085180.contigs.fasta",
        "Ia": false,
        "IIa": false,
        "IIb": false,
        "IIIa": false,
        "IVa": false,
        "IVb": false,
        "IVc": false,
        "IVd": true,
        "IVg": false,
        "IVh": false
    }
}
```

Similar to the Staphopia results, SRX085180 based on its assembly appears to be 
SCCmec Type IVd.

