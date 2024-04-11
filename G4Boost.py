#!/usr/bin/env python

import re
import sys
import operator
import pandas as pd
import xgboost as xgb
from lib.arg_parser import get_args

" ------------------------------[  Functions ]--------------------------------- "

args = get_args()


def sort_table(table, cols):
    for col in reversed(cols):
        table = sorted(table, key=operator.itemgetter(col))
    return table


def chrom_name(header):
    if not header.startswith(">"):
        #        raise Exception('FASTA header does not start with ">":\n%s' % header)
        return "noID"
    chr = re.sub("^>\s*", "", header)
    chr = re.sub("\s.*", "", chr)
    return chr


def revcomp(seq):
    complement = {"A": "T", "C": "G", "G": "C", "T": "A", "U": "A", "N": "N"}
    return "".join(complement.get(base, base) for base in reversed(seq))


def findall(seq, search):
    count = -1
    loc = 0
    newloc = 0
    while newloc > -1:
        newloc = seq[loc:].find(search)
        loc = loc + newloc + 1
        count += 1
    return count


def initialize_dataFrame():
    header = [
        "seq",
        "seq_length",
        "g4motif",
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
    data_dict = {}
    for h in header:
        data_dict[h] = []
    return data_dict


def topology(reg, seq):
    split_seq = re.split(reg, seq)
    if len(split_seq[-1]) == 0:
        gstem_base = split_seq[-2]
    else:
        gstem_base = split_seq[-1]
    g = len(gstem_base)
    loops = [len(sp_seq) - g for sp_seq in split_seq]
    loops = [lbase for lbase in loops if lbase > 0]
    maxlbase = max(loops)
    minlbase = min(loops)
    test = gstem_base
    for sp_seq in split_seq:
        if len(sp_seq) > g:
            test += sp_seq[g:].lower()
            test += gstem_base
    return [test, len(test), len(loops) + 1, g, maxlbase, minlbase]


def update_dataFrame(features, reg, seq, ref):
    [test, length, maxgstem, maxgbase, maxlbase, minlbase] = topology(reg, seq)
    features["g4motif"].append(test)
    features["length"].append(length)
    features["seq_length"].append(len(ref))
    features["loops"].append(maxgstem)
    features["G-quartet"].append(maxgbase)
    features["maxlbase"].append(maxlbase)
    features["minlbase"].append(minlbase)
    features["G"].append(int(findall(ref, "G") * 100 / len(ref)))
    features["GG"].append(int(findall(ref, "GG") * 100 / len(ref)))
    features["C"].append(int(findall(ref, "C") * 100 / len(ref)))
    features["CC"].append(int(findall(ref, "CC") * 100 / len(ref)))
    return features


def findmotifs(reg, seq, start):
    gquad_list = []
    for m in re.finditer(reg, seq):
        seq = m.group(0)
        quad_id = chrom + "_" + \
            str(m.start() + start) + "_" + str(m.end() + start)
        gquad_list.append(
            [
                chrom,
                m.start() + start,
                m.end() + start,
                quad_id,
                len(m.group(0)),
                "+",
                seq,
            ]
        )
    return gquad_list


# -----------------------------------------------------------------------------


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
features = initialize_dataFrame()

# if args.fasta != '-': output= args.fasta+'.gff'
# else: output = 'G4Boost_quadruplexes.gff'
# out=open(output, 'w')

while True:
    if not args.quiet:
        sys.stderr.write("Processing %s\n" % (chrom))
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
                        features = update_dataFrame(features, reg, seq, ref)
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

# ---------------

sys.stderr.write("Starting stability prediction!\n\n")
regressor = xgb.XGBRegressor()
classifier = xgb.XGBClassifier()
regressor.load_model(args.regressor)
classifier.load_model(args.classifier)

# classifier = pickle.load(open(args.classifier, 'rb'))
# regressor = pickle.load(open(args.regressor, 'rb'))
# selected='length seq_length G-quartet loops maxlbase minlbase G C GG CC'.split(' ')
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
# X_test=xgb.DMatrix(X_test)
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

sys.stderr.write("G4Boost completed screening!\n\n")
