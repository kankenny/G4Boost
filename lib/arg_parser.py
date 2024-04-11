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
        required=True,
    )
    parser.add_argument(
        "--classifier",
        required=False,
        default="G4Boost_classifier.json",
        help="""Use specified classifier (G4Boost_classifier.sav)""",
    )
    parser.add_argument(
        "--regressor",
        required=False,
        default="G4Boost_regressor.json",
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
        help="""Maximum number of consecutive G bases within a G-stem. Default is to report up to 1 Gs.""",
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
    parser.add_argument("--version", "-v", action="version",
                        version="%(prog)s " + VERSION)

    args = parser.parse_args()

    return args
