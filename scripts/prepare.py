"""
from Bactopia!
"""

import logging
import re
import sys
import textwrap
from pathlib import Path

import rich
import rich.console
import rich.traceback
import rich_click as click
from rich.logging import RichHandler

# Set up Rich
stderr = rich.console.Console(stderr=True)
rich.traceback.install(console=stderr, width=200, word_wrap=True, extra_lines=1)
click.rich_click.USE_RICH_MARKUP = True
click.rich_click.OPTION_GROUPS = {
    "prepare": [
        {"name": "Required Options", "options": ["--path"]},
        {
            "name": "Matching Options",
            "options": [
                "--fastq-ext",
                "--fastq-separator",
                "--pe1-pattern",
                "--pe2-pattern",
                "--ont",
                "--recursive",
                "--prefix",
            ],
        },
        {
            "name": "Additional Options",
            "options": [
                "--examples",
                "--verbose",
                "--silent",
                "--version",
                "--help",
            ],
        },
    ]
}


def search_path(path, pattern, recursive=False):
    if recursive:
        return Path(path).rglob(pattern)
    else:
        return Path(path).glob(pattern)


def get_path(fastq, abspath, prefix):
    fastq_path = str(fastq.absolute())
    if prefix:
        return fastq_path.replace(str(abspath), prefix).replace("///", "//")
    return fastq_path


def print_examples():
    print(
        textwrap.dedent(
            """
        # Example '*.fastq.gz' FASTQ files:
        prepare.py --path fastqs/
        sample          runtype         r1                            r2
        sample01        paired-end      fastqs/sample01_R1.fastq.gz   fastqs/sample01_R2.fastq.gz
        sample02        single-end      fastqs/sample02.fastq.gz
        sample03        paired-end      fastqs/sample03_R1.fastq.gz   fastqs/sample03_R2.fastq.gz

        # Example '*_001.fastq.gz' FASTQ files:
        prepare.py --path fastqs/ --fastq-ext '_001.fastq.gz'
        sample          runtype         r1                                r2
        sample01        paired-end      fastqs/sample01_R1_001.fastq.gz   fastqs/sample01_R2_001.fastq.gz
        sample02        paired-end      fastqs/sample02_R1_001.fastq.gz   fastqs/sample02_R2_001.fastq.gz
        sample03        paired-end      fastqs/sample03_R1_001.fastq.gz   fastqs/sample03_R2_001.fastq.gz

        # Example '*.fq.gz' FASTQ files:
        prepare.py --path fastqs --fastq-ext '.fq.gz'
        sample         runtype          r1                      r2   
        sample01       single-end      fastqs/sample01.fq.gz
        sample02       single-end      fastqs/sample02.fq.gz
        sample03       single-end      fastqs/sample03.fq.gz

        # Example Nanopore FASTQ files:
        prepare.py --path fastqs/ --ont
        sample          runtype       r1                         r2
        sample01        ont           fastqs/sample01.fastq.gz
        sample02        ont           fastqs/sample02.fastq.gz
        sample03        ont           fastqs/sample03.fastq.gz

        # Example changing the separator, do this if you have underscores in your sample names and you want SE reads:
        prepare.py --path ext/ --fastq-separator '.'
        sample          runtype         r1                          r2      
        sample_01       single-end      fastqs/sample_01.fastq.gz
        sample_02       single-end      fastqs/sample_02.fastq.gz
        sample_03       single-end      fastqs/sample_03.fastq.gz
    """
        )
    )
    sys.exit(0)


@click.command(epilog='Created by rpetit3 https://github.com/bactopia/bactopia-py then butchered by Ammar Aziz.')
@click.option(
    "--path", "-p", required=True, help="Directory where FASTQ files are stored"
)
@click.option(
    "--fastq-ext",
    "-f",
    default=".fastq.gz",
    show_default=True,
    help="Extension of the FASTQs",
)
@click.option(
    "--pe1-pattern",
    "-1",
    default="[Aa]|[Rr]1|1",
    show_default=True,
    help="Designates difference first set of paired-end reads",
)
@click.option(
    "--pe2-pattern",
    "-2",
    default="[Bb]|[Rr]2|2",
    show_default=True,
    help="Designates difference second set of paired-end reads",
)
@click.option(
    "--fastq-separator",
    "-s",
    default="_",
    show_default=True,
    help="Split FASTQ name on the last occurrence of the separator",
)
@click.option(
    "--recursive", 
    "-r", 
    is_flag=True, 
    help="Directories will be traversed recursively"
)
@click.option(
    "--ont",
    "-o",
    is_flag=True,
    help="Single-end reads should be treated as Oxford Nanopore reads",
)
@click.option(
    "--prefix", 
    "-x",
    default=None, 
    help="Prefix to add to the path")
@click.option(
    "--examples", 
    "-e",
    is_flag=True, 
    help="Print example usage")
@click.option(
    "--verbose",
    "-v",
    is_flag=True, 
    help="Increase the verbosity of output")
@click.option("--silent", 
    is_flag=True, 
    help="Only critical errors will be printed")
def prepare(
    path,
    fastq_ext,
    pe1_pattern,
    pe2_pattern,
    fastq_separator,
    recursive,
    ont,
    prefix,
    examples,
    verbose,
    silent,
):
    """
    \b
    Create a 'file of filenames' (FOFN) of samples for pipeline input
    """
    # Setup logs
    logging.basicConfig(
        # format="%(asctime)s:%(name)s:%(levelname)s - %(message)s",
        format="%(levelname)s - %(message)s",
        datefmt="!!!",
        # datefmt="%Y-%m-%d %H:%M:%S",
        handlers=[
            RichHandler(rich_tracebacks=True, 
                        console=rich.console.Console(stderr=True), 
                        show_level=False, 
                        markup=True)
        ], 
    )
    logging.getLogger().setLevel(
        logging.ERROR if silent else logging.DEBUG if verbose else logging.INFO
    )

    abspath = Path(path).absolute()
    SAMPLES = {}

    # Match FASTQS
    for fastq in search_path(abspath, f"*{fastq_ext}", recursive=recursive):
        fastq_name = fastq.name.replace(fastq_ext, "")
        # Split the fastq file name on separator
        # Example MY_FASTQ_R1.rsplit('_', 1) becomes ['MY_FASTQ', 'R1'] (PE)
        # Example MY_FASTQ.rsplit('_', 1) becomes ['MY_FASTQ'] (SE)
        split_vals = fastq_name.rsplit(fastq_separator, 1)
        sample_name = split_vals[0]
        if sample_name not in SAMPLES:
            SAMPLES[sample_name] = {
                "pe": {"r1": [], "r2": []},
                "se": [],
                # "assembly": [],
            }

        if len(split_vals) == 1:
            # single-end
            SAMPLES[sample_name]["se"].append(get_path(fastq, abspath, prefix))
        else:
            # paired-end
            pe1 = re.compile(pe1_pattern)
            pe2 = re.compile(pe2_pattern)
            if pe1.match(split_vals[1]):
                SAMPLES[sample_name]["pe"]["r1"].append(
                    get_path(fastq, abspath, prefix)
                )
            elif pe2.match(split_vals[1]):
                SAMPLES[sample_name]["pe"]["r2"].append(
                    get_path(fastq, abspath, prefix)
                )
            else:
                logging.error(
                    f'Could not determine read set for "{fastq_name}".'
                )
                logging.error(
                    f"Found {split_vals[1]} expected (R1: {pe1_pattern} or R2: {pe2_pattern}). \nIf input is SE and your sample name contains underscores, add --fastq-separator '-'"
                )
                logging.error(
                    "Please use --pe1-pattern and --pe2-pattern to correct and try again."
                )
                sys.exit(1)

    FOFN = []
    for sample, vals in sorted(SAMPLES.items()):
        r1_reads = vals["pe"]["r1"]
        r2_reads = vals["pe"]["r2"]
        se_reads = vals["se"]
        errors = []
        is_single_end = False
        multiple_read_sets = False
        pe_count = len(r1_reads) + len(r2_reads)

        if len(r1_reads) != len(r2_reads):
            # PE reads must be a pair
            errors.append(
                f'[bold yellow]"{sample}" must have equal paired-end read sets (R1 has {len(r1_reads)} and R2 has {len(r2_reads)}), please check.[/bold yellow] [spring_green3]!!![/spring_green3]'
            )
        elif pe_count > 2:
            # PE reads must be a pair
            # if merge:
            #     multiple_read_sets = True
            # else:
            errors.append(
                f'[bold yellow]"{sample}" cannot have more than two paired-end FASTQ, please check loggic and input. [/bold yellow][spring_green3]!!![/spring_green3]'
            )
        if ont:
            if not pe_count and len(se_reads):
                is_single_end = True
            ### needed??
            
            elif pe_count and len(se_reads): #and not hybrid and not short_polish:
                errors.append(
                    f'[bold yellow]"{sample}" cannot have paired and single-end FASTQs with --ont flag, please check. I think you hit a section from bactopia-prepare.py! [spring_green3]!!![/spring_green3]'
                )
        else:
            if len(se_reads) > 1:
                # Can't have multiple SE reads
                # if merge:
                #     multiple_read_sets = True
                # else:
                errors.append(
                    f'[bold yellow]"{sample}" has more than two single-end FASTQs, please check. [/bold yellow][spring_green3]!!![/spring_green3]'
                )
            elif pe_count and len(se_reads):
                # Can't have SE and PE reads unless long reads
                errors.append(
                    f'[bold yellow]"{sample}" has paired and single-end FASTQs, please check. [/bold yellow][spring_green3]!!![/spring_green3]'
                )

        if errors:
            logging.error("\n".join(errors))
        else:
            runtype = ""
            r1 = ""
            r2 = ""

            if pe_count:
                if multiple_read_sets:
                    ### needed???
                    if ont:
                        pass
                        # if hybrid:
                        #     runtype = "hybrid-merge-ont"
                        # elif short_polish:
                        #     runtype = "short_polish-merge-ont"
                    else:
                        runtype = "merge-pe"
                    r1 = ",".join(sorted(r1_reads))
                    r2 = ",".join(sorted(r2_reads))
                else:
                    runtype = "paired-end"
                    r1 = r1_reads[0]
                    r2 = r2_reads[0]

            if se_reads:
                if ont and is_single_end:
                    runtype = "ont"
                    r1 = se_reads[0]
                else:
                    # if multiple_read_sets:
                    #     runtype = "merge-se"
                    #     r1 = ",".join(se_reads)
                    # else:
                    runtype = "single-end"
                    r1 = se_reads[0]

            FOFN.append(
                [sample, runtype, r1, r2]
            )

    if FOFN:
        print("#sample,runtype,r1,r2")
        for line in FOFN:
            print(",".join(line))
    else:
        logging.error(
            f"Unable to find any samples in {path}. Please try adjusting the following parameters to fit your needs."
        )
        logging.error("Values Used:")
        logging.error(f"    --fastq-ext => {fastq_ext}")
        logging.error(f"    --fastq-separator => {fastq_separator}")
        logging.error(f"    --pe1-pattern => {pe1_pattern}")
        logging.error(f"    --pe2-pattern => {pe2_pattern}")
        logging.error("")
        logging.error(
            "You can also use '--examples' to see a few examples of using prepare.py"
        )
        sys.exit(1)


def main():
    if len(sys.argv) == 1:
        prepare.main(["--help"])
    elif "--examples" in sys.argv:
        print_examples()
    else:
        prepare()


if __name__ == "__main__":
    main()
