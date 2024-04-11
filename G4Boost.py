import re
import sys
import pandas as pd
import xgboost as xgb
from collections import defaultdict

from lib.arg_parser import get_args
from lib.util import (
    chrom_name,
    revcomp,
    update_dataFrame,
    sort_table
)

" -----------------------------[  Functions ]-------------------------------- "

args = get_args()


if args.fasta != "-":
    ref_seq_fh = open(args.fasta)
    output = args.fasta + ".gff"
else:
    ref_seq_fh = sys.stdin
    output = "G4Boost_quadruplexes.gff"

ref_seq = []
line = ref_seq_fh.readline()
chrom = chrom_name(line)
if chrom != "noID":
    line = ref_seq_fh.readline()
else:
    chrom = line.strip()
gquad_list = []
eof = False

gb = range(args.minG, args.maxG + 1)[::-1]
gs = range(3, args.loops + 1)[::-1]
longest = (args.maxG + args.maxloop) * args.loops + args.maxG
features = defaultdict(list)


while True:
    if not args.quiet:
        print("Processing %s\n" % (chrom))
    while line.startswith(">") is False:
        ref_seq.append(line.strip())
        line = ref_seq_fh.readline()
        if line == "":
            eof = True
            break
    ref_seq = "".join(ref_seq)
    ref_seq = ref_seq.upper().replace("U", "T")
    rev_ref_seq = revcomp(ref_seq)
    seqlen = len(ref_seq)
    for g in gb:
        for s in gs:
            gstem_base = ""
            for i in range(g):
                gstem_base += "G"
            reg = ""
            for i in range(s):
                reg += "([gG]{%d}\w{%d,%d})" % (g, args.minloop, args.maxloop)
            reg += "([gG]{%d})" % (g)
            for m in re.finditer(reg, ref_seq):
                seq = m.group(0)
                start = m.start()
                end = m.end()
                if len(ref_seq) > longest:
                    ref = seq
                else:
                    ref = ref_seq
                quad_id = chrom + "_" + str(m.start()) + "_" + str(m.end())
                gquad_list.append(
                    [chrom, start, end, quad_id, len(seq), "+", seq])
                if seq not in features["g4motif"]:
                    features = update_dataFrame(features, reg, seq, ref)
                    features["seq"].append(chrom)
                temp = ""
                for i in range(start, end):
                    temp += "N"
                ref_seq = ref_seq[:start] + temp + ref_seq[end:]
            if args.noreverse is False:
                for m in re.finditer(reg, rev_ref_seq):
                    seq = m.group(0)
                    start = m.start()
                    end = m.end()
                    if len(rev_ref_seq) > longest:
                        ref = seq
                    else:
                        ref = rev_ref_seq
                    quad_id = chrom + "_" + str(m.start()) + "_" + str(m.end())
                    gquad_list.append(
                        [
                            chrom,
                            seqlen - end,
                            seqlen - start,
                            quad_id,
                            len(seq),
                            "-",
                            seq,
                        ]
                    )
                    if seq not in features["g4motif"]:
                        features = update_dataFrame(
                            features, reg, seq, ref)
                        features["seq"].append(chrom)
                    temp = ""
                    for i in range(start, end):
                        temp += "N"
                    rev_ref_seq = rev_ref_seq[:start] + \
                        temp + rev_ref_seq[end:]
            gquad_sorted = sort_table(gquad_list, (1, 2, 3))
            gquad_list = []
            for xline in gquad_sorted:
                xline = "\t".join([str(x) for x in xline])
                with open(output, "a") as out:
                    out.write(xline + "\n")
    if eof:
        break
    chrom = chrom_name(line)
    ref_seq = []
    line = ref_seq_fh.readline()
    if line == "":
        break

print("Starting stability prediction!\n\n")
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

if args.fasta != "-":
    output = args.fasta + ".g4scores.csv"
else:
    output = "G4Boost_quadruplexes.g4.csv"
features.to_csv(output, sep="\t", index=False)

print("G4Boost completed screening!\n\n")
