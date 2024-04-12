import re
import sys
from collections import defaultdict
from lib.util import chrom_name, revcomp, update_dataFrame, sort_table


def process_sequences(args):
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
                    reg += "([gG]{%d}\w{%d,%d})" % (g,
                                                    args.minloop, args.maxloop)
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
                        quad_id = chrom + "_" + \
                            str(m.start()) + "_" + str(m.end())
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
    return features
