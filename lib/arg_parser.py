import argparse


VERSION = "0.1.0"


def get_args():
    parser = argparse.ArgumentParser(
        description="""
            DESCRIPTION

            EXAMPLE:
                    """,
        formatter_class=argparse.RawTextHelpFormatter,
    )

    parser.add_argument(
        "--fasta",
        "-f",
        type=str,
        help="""Input fasta file to search. Use '-' to read the file from stdin.""",
        required=False,
        default="",
    )
    parser.add_argument(
        "--sequence",
        "-s",
        type=str,
        help="""Nucleotide Sequence""",
        required=False,
        default="",
    )
    parser.add_argument(
        "--sequence_id",
        "-sid",
        type=str,
        help="""Nucleotide Sequence Identifier""",
        required=False,
        default="",
    )
    parser.add_argument(
        "--classifier",
        required=False,
        default="lib/models/G4Boost_classifier.json",
        help="""Use specified classifier (G4Boost_classifier.sav)""",
    )
    parser.add_argument(
        "--regressor",
        required=False,
        default="lib/models/G4Boost_regressor.json",
        help="""Use specified classifier (G4Boost_regressor.sav)""",
    )
    parser.add_argument(
        "--maxloop",
        "-N",
        type=int,
        required=False,
        default=12,
        help="""Maximum length of the loop. Default is to report up to 12nt.""",
    )
    parser.add_argument(
        "--minloop",
        "-n",
        type=int,
        required=False,
        default=1,
        help="""Minimum length of the loop. Default is to report up to 1nt.""",
    )
    parser.add_argument(
        "--maxG",
        "-G",
        type=int,
        required=False,
        default=7,
        help="""Maximum number of consecutive G bases within a G-stem. Default is to report up to 7 Gs.""",
    )
    parser.add_argument(
        "--minG",
        "-g",
        type=int,
        required=False,
        default=1,
        help="Maximum number of consecutive G bases within a G-stem. Default is to report up to 1 Gs."
        "",
    )
    parser.add_argument(
        "--loops",
        "-l",
        type=int,
        required=False,
        default=11,
        help="""Maximum number of flexible loops separating the G-stems. Default is to report up to 11 Gs.""",
    )
    parser.add_argument(
        "--noreverse",
        action="store_true",
        help="""Do not search the reverse complement of the input fasta.""",
    )
    parser.add_argument(
        "--quiet",
        "-q",
        action="store_true",
        help="""Do not print progress report (i.e. sequence names as they are scanned).""",
    )
    parser.add_argument(
        "--version", "-v", action="version", version="%(prog)s " + VERSION
    )
    parser.add_argument(
        "--output",
        "-o",
        type=str,
        required=False,
        default="G4s",
        help="""Output file name of G4Boost.""",
    )

    args = parser.parse_args()

    if (len(args.fasta) == 0) and (len(args.sequence) == 0):
        parser.error("Exactly one of --fasta or --sequence must be provided.")

    if (len(args.fasta) != 0) and (len(args.sequence) != 0):
        parser.error(
            "FASTA file and Sequence both are provided. Inference will default to the FASTA file")

    if (len(args.sequence) != 0) and (len(args.sequence_id) == 0):
        parser.error(
            "Sequence ID must be provided when sequence argument is provided.")

    return args
