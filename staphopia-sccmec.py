#! /usr/bin/env python3
"""

"""
import json
import logging
import os
import sys
from collections import OrderedDict


PROGRAM = os.path.basename(sys.argv[0])
VERSION = "0.0.2"
DESCRIPTION = "A standalone version of Staphopia's SCCmec typing method."
STDOUT = 11
STDERR = 12
logging.addLevelName(STDOUT, "STDOUT")
logging.addLevelName(STDERR, "STDERR")


def validate_requirements():
    """Validate the required programs are available, if not exit (1)."""
    from shutil import which
    programs = {
        'makeblastdb': which('makeblastdb'), 'blastn': which('blastn'),
        'tblastn': which('tblastn')
    }

    missing = False
    for prog, path in programs.items():
        if path:
            logging.debug(f'{prog}: command found!')
        else:
            logging.error(f'{prog}: command not found in PATH.')
            missing = True

    if missing:
        logging.error("Program requirement missing, exiting")
        sys.exit(1)


def validate_datasets(data_dir, input_assembly):
    """Validate required reference datasets are available, if not exit (1)."""
    datasets = {
        'primers': f'{data_dir}/primers.fasta',
        'proteins': f'{data_dir}/proteins.fasta',
        'subtypes': f'{data_dir}/subtypes.fasta',
        'mapping': f'{data_dir}/sccmec_cassettes.bwt',
        'input_assembly': input_assembly
    }

    missing = False
    for dataset, path in datasets.items():
        if os.path.exists(path):
            logging.debug(f'{dataset} ({path}) found!')
        else:
            logging.error(f'{dataset} ({path}) not found, please verify path.')
            missing = True

    if missing:
        logging.error("Dataset requirement missing, exiting")
        sys.exit(1)


def set_log_level(error, debug):
    """Set the output log level."""
    return logging.ERROR if error else logging.DEBUG if debug else logging.INFO


def get_log_level():
    """Return logging level name."""
    return logging.getLevelName(logging.getLogger().getEffectiveLevel())


def execute(cmd, directory=os.getcwd(), capture=False, capture_stderr=False,
            stdout_file=None, stderr_file=None, silent=True):
    """A simple wrapper around executor."""
    from executor import ExternalCommand
    command = ExternalCommand(
        cmd, directory=directory, capture=capture,
        capture_stderr=capture_stderr, stdout_file=stdout_file,
        stderr_file=stderr_file, silent=silent
    )

    command.start()
    if get_log_level() == 'DEBUG':
        logging.log(STDOUT, command.decoded_stdout)
        logging.log(STDERR, command.decoded_stderr)

    if capture:
        return command.decoded_stdout


def makeblastdb(input, outdir, prefix):
    """ Make a BLAST database for a given input FASTA. """
    blastdb = f'{outdir}/blastdb'
    execute(f'mkdir -p {blastdb}')
    execute((f'makeblastdb -in {input} -input_type fasta -dbtype "nucl" '
             f'-title "{prefix}" -out {blastdb}/{prefix}'))
    return f'{blastdb}/{prefix}'


def blastn(query, blastdb, output):
    """ BLAST SCCmec related primers against an input assembly. """
    execute((
        f'blastn -db {blastdb} -max_target_seqs 1 -dust no -word_size 7 '
        f'-perc_identity 100 -outfmt 15 -query {query} > {output}'
    ))
    return read_blast_json(output)


def tblastn(query, blastdb, output, cpu=1):
    """ BLAST SCCmec related proteins against an input assembly. """
    execute((
        f'tblastn -db {blastdb} -outfmt 15 -query {query} -evalue 0.0001 '
        f'-outfmt 15 -num_threads {cpu} -max_target_seqs 1 > {output}'
    ))
    return read_blast_json(output)


def mapping(fastq, ):
    """ Map FASTQ reads against reference SCCmec cassettes. """
    """
    bwa-align.sh "!{fq}" !{staphopia_data}/sccmec/sccmec_cassettes !{read_length} !{cpu} "!{p}" "!{n}"
    samtools view -bS bwa.sam | samtools sort -o sccmec.bam -
    genomeCoverageBed -ibam sccmec.bam -d | gzip --best - > cassette-coverages.gz
    cp .command.err sccmec-mapping-stderr.log
    cp .command.out sccmec-mapping-stdout.log
    """


def read_blast_json(json_file):
    """ Process input BLAST JSON and return hits as a list. """
    hits = []
    json_data = None
    with open(json_file, 'rt') as blast_fh:
        json_data = json.load(blast_fh)

    for entry in json_data['BlastOutput2']:
        hit = entry['report']['results']['search']
        if len(hit['hits']):
            # Only storing the top hit
            hsp = hit['hits'][0]['hsps'][0]

            # Includes mismatches and gaps
            mismatch = hsp['align_len'] - hsp['identity']

            # Hamming distance
            hd = mismatch
            if hit['query_len'] > hsp['align_len']:
                # Include those bases that weren't aligned
                hd = hit['query_len'] - hsp['align_len'] + mismatch

            hits.append({
                'title': hit['query_title'],
                'length': hit['query_len'],
                'bitscore': int(hsp['bit_score']),
                'evalue': hsp['evalue'],
                'identity': hsp['identity'],
                'mismatch': mismatch,
                'gaps': hsp['gaps'],
                'hamming_distance': hd,
                'query_from': hsp['query_from'],
                'query_to': hsp['query_to'],
                'hit_from': hsp['hit_from'],
                'hit_to': hsp['hit_to'],
                'align_len': hsp['align_len'],
                'qseq': hsp['qseq'],
                'hseq': hsp['hseq'],
                'midline': hsp['midline']
            })
    return hits


def max_primer_hamming_distance():
    return ({
        "IS1272-F": 20, "IS7": 24, "IS5": 23, "IS2": 25, "mecI-R": 19,
        "mecI-F": 21, "mI6": 23, "mI4": 20, "mI3": 25, "mecR1-R": 24,
        "mcR5": 21, "mcR3": 20, "mcR2": 20, "mA7": 23, "mA6": 18,
        "mA2": 21, "mA1": 21, "ccrCf-B": 27, "ccrCr-B": 28, "ccrCr-A": 22,
        "ccrCf-A": 20, "ccrB4": 20, "ccrB": 22, "ccrA4": 20, "ccrA3": 23,
        "ccrA2": 24, "ccrA1": 24
    })


def predict_type_by_primers(prefix, blast_results, hamming_distance=False):
    dist = max_primer_hamming_distance()
    primers = OrderedDict()
    for key in dist:
        primers[key] = False

    for hit in blast_results:
        if '|' in hit['title']:
            name, primer = hit['title'].split('|')
        else:
            name = hit['title']

        # Differentiate the mec class primers
        if name in ['mecR1', 'mecI', 'mecA', 'IS1272', 'IS431']:
            name = primer
        elif name == 'ccrCf':
            if int(hit['length']) == 27:
                name = 'ccrCf-B'
            else:
                name = 'ccrCf-A'
        elif name == 'ccrCr':
            if int(hit['length']) == 28:
                name = 'ccrCr-B'
            else:
                name = 'ccrCr-A'

        dist[name] = int(hit['hamming_distance'])
        if int(hit['hamming_distance']) == 0:
            # Require perfect matches
            primers[name] = True

    # Determine ccrC
    dist['ccrC'] = min((dist['ccrCr-B'] + dist['ccrCf-B']),
                       (dist['ccrCr-A'] + dist['ccrCf-A']))
    if ((primers['ccrCr-B'] and primers['ccrCf-B']) or
            (primers['ccrCr-A'] and primers['ccrCf-A'])):
        primers['ccrC'] = True
    else:
        primers['ccrC'] = False

    # Determine mec class
    mec_class = {'meca': False, 'A': False, 'B': False, 'C': False,
                 'AB': False, 'ABC': False}
    mec_dist = {'meca': 0, 'A': 0, 'B': 0, 'C': 0, 'AB': 0, 'ABC': 0}

    if hamming_distance:
        mec_dist['meca'] = dist['mA1'] + dist['mA2']
        # Class A
        mec_dist['A'] = mec_dist['meca'] + min(
            (dist['mI4'] + dist['mI3'] + dist['mcR2'] + dist['mcR5']),
            (dist['mI4'] + dist['mcR3'])
        )

        # Class B
        mec_dist['B'] = mec_dist['meca'] + dist['IS5'] + dist['mA6']

        # Class C
        mec_dist['C'] = mec_dist['meca'] + dist['IS2']

        # Class A,B
        mec_dist['AB'] = (mec_dist['meca'] + dist['mecI-R'] +
                          dist['mecI-F'] + dist['IS1272-F'] +
                          dist['mecR1-R'])

        # Class A,B,C
        mec_dist['ABC'] = (mec_dist['meca'] + dist['mecI-R'] +
                           dist['mecI-F'] + dist['IS1272-F'] +
                           dist['mecR1-R'])

    # True/False Classes
    if primers['mA1'] and primers['mA2']:
        mec_class['meca'] = True
        # Class A
        if (primers['mI4'] and primers['mI3'] and primers['mcR2'] and
                primers['mcR5']) or (primers['mI4'] and primers['mcR3']):
                mec_class['A'] = True

        # Class B
        if primers['IS5'] and primers['mA6']:
            mec_class['B'] = True

        # Class C
        if primers['IS2']:
            mec_class['C'] = True

        # Class A,B
        if (primers['mecI-R'] and primers['mecI-F'] and
                primers['IS1272-F'] and primers['mecR1-R']):
            mec_class['A'] = True
            mec_class['B'] = True
            mec_class['AB'] = True

        # Class A,B,C
        if (primers['mI6'] and primers['IS7'] and
                primers['IS2'] and primers['mA7']):
            mec_class['A'] = True
            mec_class['B'] = True
            mec_class['C'] = True
            mec_class['ABC'] = True

    mec = OrderedDict([
        ('sample', prefix),
        ('I', False), ('II', False), ('III', False), ('IV', False),
        ('V', False), ('VI', False), ('VII', False), ('VIII', False),
        ('IX', False), ('meca', mec_class['meca'])
    ])

    if hamming_distance:
        mec['meca'] = mec_dist['meca']
        mec['I'] = dist['ccrA1'] + dist['ccrB'] + mec_dist['B']
        mec['II'] = dist['ccrA2'] + dist['ccrB'] + mec_dist['A']
        mec['III'] = dist['ccrA3'] + dist['ccrB'] + mec_dist['A']
        mec['IV'] = dist['ccrA2'] + dist['ccrB'] + mec_dist['B']
        mec['V'] = dist['ccrC'] + mec_dist['C']
        mec['VI'] = dist['ccrA4'] + dist['ccrB4'] + mec_dist['B']
        mec['VII'] = dist['ccrC'] + mec_dist['C']
        mec['VIII'] = dist['ccrA4'] + dist['ccrB4'] + mec_dist['A']
        mec['IX'] = dist['ccrA1'] + dist['ccrB'] + mec_dist['C']
    else:
        if primers['ccrA1'] and primers['ccrB'] and mec_class['B']:
            mec['I'] = True

        if primers['ccrA2'] and primers['ccrB'] and mec_class['A']:
            mec['II'] = True

        if primers['ccrA3'] and primers['ccrB'] and mec_class['A']:
            mec['III'] = True

        if primers['ccrA2'] and primers['ccrB'] and mec_class['B']:
            mec['IV'] = True

        if primers['ccrC'] and mec_class['C']:
            mec['V'] = True

        if primers['ccrA4'] and primers['ccrB4'] and mec_class['B']:
            mec['VI'] = True

        if primers['ccrC'] and mec_class['C']:
            mec['VII'] = True

        if primers['ccrA4'] and primers['ccrB4'] and mec_class['A']:
            mec['VIII'] = True

        if primers['ccrA1'] and primers['ccrB'] and mec_class['C']:
            mec['IX'] = True

    return mec


def max_subtype_hamming_distance():
    return ({
        "ivh-r": 20, "ivh-f": 20, "ivg-r": 20, "ivg-f": 20, "ivd-r": 20,
        "ivd-f": 20, "ivce-r": 20, "ivce-f": 20, "ivbf-r": 20, "ivbf-f": 20,
        "iva-r": 20, "iva-f": 20, "iiia-3ab": 20, "iiia-3a1": 20,
        "iib-2b4": 22, "iib-2b3": 21, "iia-kdpB2": 22, "iia-kdpB1": 22,
        "ia-1a4": 21, "ia-1a3": 23
    })


def predict_subtype_by_primers(prefix, blast_results, hamming_distance=False):
    dist = max_subtype_hamming_distance()
    primers = {}
    for key in dist:
        primers[key] = False

    for hit in blast_results:
        name = hit['title']
        if int(hit['hamming_distance']) == 0:
            primers[name] = True
        dist[name] = int(hit['hamming_distance'])

    subtypes = OrderedDict([
        ('sample', prefix),
        ('Ia', False), ('IIa', False), ('IIb', False), ('IIIa', False),
        ('IVa', False), ('IVb', False), ('IVc', False), ('IVd', False),
        ('IVg', False), ('IVh', False),
    ])

    if hamming_distance:
        subtypes['Ia'] = dist['ia-1a3'] + dist['ia-1a4']
        subtypes['IIa'] = dist['iia-kdpB1'] + dist['iia-kdpB2']
        subtypes['IIb'] = dist['iib-2b3'] + dist['iib-2b4']
        subtypes['IIIa'] = dist['iiia-3a1'] + dist['iiia-3ab']
        subtypes['IVa'] = dist['iva-f'] + dist['iva-r']
        subtypes['IVb'] = dist['ivbf-f'] + dist['ivbf-r']
        subtypes['IVc'] = dist['ivce-f'] + dist['ivce-r']
        subtypes['IVd'] = dist['ivd-f'] + dist['ivd-r']
        subtypes['IVg'] = dist['ivg-f'] + dist['ivg-r']
        subtypes['IVh'] = dist['ivh-f'] + dist['ivh-r']
    else:
        if primers['ia-1a3'] and primers['ia-1a4']:
            subtypes['Ia'] = True

        if primers['iia-kdpB1'] and primers['iia-kdpB2']:
            subtypes['IIa'] = True

        if primers['iib-2b3'] and primers['iib-2b4']:
            subtypes['IIb'] = True

        if primers['iiia-3a1'] and primers['iiia-3ab']:
            subtypes['IIIa'] = True

        if primers['iva-f'] and primers['iva-r']:
            subtypes['IVa'] = True

        if primers['ivbf-f'] and primers['ivbf-r']:
            subtypes['IVb'] = True

        if primers['ivce-f'] and primers['ivce-r']:
            subtypes['IVc'] = True

        if primers['ivd-f'] and primers['ivd-r']:
            subtypes['IVd'] = True

        if primers['ivg-f'] and primers['ivg-r']:
            subtypes['IVg'] = True

        if primers['ivh-f'] and primers['ivh-r']:
            subtypes['IVh'] = True

    return subtypes


def print_predictions(cassettes, subtypes, output, quiet=False):
    """ Output the predictions in JSON format. """
    prediction = {
        'cassette': cassettes,
        'subtype': subtypes
    }
    with open(output, 'w') as json_fh:
        json.dump(prediction, json_fh, indent=4)

    if not quiet:
        print(json.dumps(prediction, indent=4))


if __name__ == '__main__':
    import argparse as ap
    import textwrap
    parser = ap.ArgumentParser(
        prog=PROGRAM,
        conflict_handler='resolve',
        description=(
            f'{PROGRAM} (v{VERSION}) - {DESCRIPTION}'
        ),
        formatter_class=ap.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent(f'''
            example usage:
              {PROGRAM} ST2019.fasta sccmec_dir --prefix sample-2019
              {PROGRAM} test/GCF_001580515.1.fna data/ ./
        ''')
    )
    parser = ap.ArgumentParser(prog='staphopia-sccmec.py',
                               conflict_handler='resolve',
                               description='Determine SCCmec Type/SubType')
    group1 = parser.add_argument_group('Options', '')
    group1.add_argument('assembly', metavar="ASSEMBLY|PRIMERS_JSON", type=str,
                        help=('Input assembly (FASTA format) to predict SCCmec. Or, '
                              '"primers.json" from Staphopia, this requires '
                              '"--staphopia"'))
    group1.add_argument('sccmec_data', metavar="SCCMEC_DATA|SUBTYPES_JSON", type=str,
                        help=('Directory where SCCmec reference data is stored. Or, '
                              '"subtypes.json" from Staphopia, this requires '
                              '"--staphopia"'))
    group1.add_argument('outdir', metavar="OUTPUT_DIR", type=str,
                        help='Directory to output results to')
    group1.add_argument('--prefix', metavar="STR", type=str,
                        help=('Prefix (e.g. sample name) to use for outputs '
                              '(Default: assembly file without extension)'))
    group1.add_argument('--staphopia', action='store_true', 
                        help=('Inputs are results from Staphopia. The first '
                              'arguement (ASSEMBLY) should be "primers.json" '
                              'and the second arguement (SCCMEC_DATA) should '
                              'be "subtypes.json". These files are found in '
                              '"analyses/sccmec" folder of Staphopia results.'))
    group1.add_argument('--hamming', action='store_true',
                        help='Report the hamming distance of each type.')
    group1.add_argument('--cpus', metavar='INT', type=int, default=1,
                        help='Number of processors to use.')
    group1.add_argument('--keep_files', action='store_true',
                        help='Keep all output files (Default: remove blastdb.')
    group1.add_argument('--debug', action='store_true',
                        help='Print debug related text.')
    group1.add_argument('--quiet', action='store_true',
                        help='Only critical errors will be printed.')
    group1.add_argument('--depends', action='store_true',
                        help='Verify dependencies are installed/found.')
    group1.add_argument('--version', action='version',
                        version=f'{PROGRAM} {VERSION}')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()

    # Setup logs
    FORMAT = '%(asctime)s:%(name)s:%(levelname)s - %(message)s'
    logging.basicConfig(format=FORMAT, datefmt='%Y-%m-%d %H:%M:%S',)
    logging.getLogger().setLevel(set_log_level(args.quiet, args.debug))

    # Check dependencies
    if args.depends:
        validate_requirements()
        validate_datasets(args.sccmec_data, args.assembly)
        sys.exit(0)
    else:
        validate_requirements()
        if not args.staphopia:
            validate_datasets(args.sccmec_data, args.assembly)

    # Get prefix for file names
    prefix = os.path.splitext(os.path.basename(args.assembly))[0]
    if args.prefix:
        prefix = args.prefix
    outdir = f'{args.outdir}/{prefix}'
    execute(f'mkdir -p {outdir}')

    primer_hits = None
    subtype_hits = None
    if not args.staphopia:
        # Make temporary BLAST database
        logging.info(f'Make BLAST database for {args.assembly}')
        blastdb = makeblastdb(args.assembly, outdir, prefix)

        # BLAST SCCmec Primers (including subtypes)
        logging.info(f'BLAST SCCmec primers against {args.assembly}')
        primer_hits = blastn(f'{args.sccmec_data}/primers.fasta', blastdb,
                             f'{outdir}/primers.json')
        subtype_hits = blastn(f'{args.sccmec_data}/subtypes.fasta', blastdb,
                              f'{outdir}/subtypes.json')
    else:
        # Read Staphopia outputs
        logging.info(f'Make BLAST database for {args.assembly}')
        primer_hits = read_blast_json(args.assembly)
        subtype_hits = read_blast_json(args.sccmec_data)

    primer_prediction = predict_type_by_primers(prefix, primer_hits,
                                                hamming_distance=args.hamming)
    subtype_prediction = predict_subtype_by_primers(prefix, subtype_hits,
                                                    hamming_distance=args.hamming)

    primer_out = f'{outdir}/sccmec-primer-type.json'
    logging.info(f'Writing predicted SCCmec type based on primers to {primer_out}')
    print_predictions(primer_prediction, subtype_prediction, primer_out,
                      quiet=args.quiet)

    # Clean up
    if not args.keep_files:
        execute(f'rm -rf {outdir}/blastdb')
