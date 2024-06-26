import pandas as pd
import xgboost as xgb
from lib.arg_parser import get_args
from lib.preprocess import process_sequences


def predict_stability(features, args):
    print("Starting stability prediction!")

    regressor = xgb.XGBRegressor()
    classifier = xgb.XGBClassifier()

    regressor.load_model(args.regressor)
    classifier.load_model(args.classifier)

    selected = [
        "seq_length",
        "length",
        "loops",
        "G-quartet",
        "maxlbase",
        "minlbase",
        "G",
        "C",
        "GG",
        "CC",
    ]
    features = pd.DataFrame.from_dict(features)
    X_test = features[selected]

    g4_pred = classifier.predict(X_test)
    g4_pred_proba = classifier.predict_proba(X_test)[:, 1]
    mfe_pred = regressor.predict(X_test)

    features["g4_pred"] = g4_pred
    features["g4_prob"] = g4_pred_proba
    features["mfe_pred"] = mfe_pred
    features["loops"] = [loop - 1 for loop in features["loops"]]

    features.to_csv(f"output/{args.output}_preds.csv", sep="\t", index=False)

    print("G4Boost completed screening!")


def main():
    args = get_args()
    features = process_sequences(args)
    predict_stability(features, args)


def test():
    import os
    import sys

    fastas = [
        "APC Exons NR_104673.1.exons .fa",
        "APC Whole Genome Sequence.fasta",
        "BRCA1 Exons.fa",
        "BRCA1 Whole genome sequence.fasta",
        "BRCA2_NM_000059.4.exons.fa",
        "BRCA2_sequence.fasta",
        "EGFR_NM_005228.5.exons.fa",
        "EGFR_sequence.fasta",
        "FTO_NM_001080432.3.exons.fa",
        "FTO_sequence.fasta",
        "GAPDH_NM_002046.7.exons.fa",
        "GAPDH_sequence.fasta",
        "PIK3CA Exons XM_006713658.5.exons.fa",
        "PIK3CA whole genome sequence.fasta",
        "PTEN Exons  NM_001126049.2.exons.fa",
        "PTEN Whole genome sequence.fasta",
        "TNF Exons NM_130463.4.exons.fa",
        "TNF Whole Genome Sequence.fasta",
        "TP53_NM_000546.6.exons.fa",
        "TP53_sequence.fasta"
    ]

    for fasta in fastas:
        file_path = os.path.join("example", fasta)
        sys.argv = ["G4Boost.py", "--fasta", file_path, "--output", fasta]
        print(sys.argv)
        main()


if __name__ == "__main__":
    import sys

    try:
        if sys.argv[1] == "test":
            test()
        else:
            main()
    except IndexError:
        main()
