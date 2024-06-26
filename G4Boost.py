import pandas as pd
import xgboost as xgb
from lib.arg_parser import get_args
from lib.preprocess import process_sequences
from tests.integration import inference_on_sample_fasta_files


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


if __name__ == "__main__":
    import sys

    try:
        if sys.argv[1] == "test":
            inference_on_sample_fasta_files()
        else:
            main()
    except IndexError:
        main()
