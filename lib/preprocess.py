import re
from collections import defaultdict
from lib.util import chrom_name, revcomp, update_dataFrame, sort_table


def process_sequences(args):
    ref_seq_file = open(args.fasta)

    ref_seq = []
    line = ref_seq_file.readline()

    chrom = chrom_name(line)
    if chrom != "noID":
        line = ref_seq_file.readline()
    else:
        chrom = line.strip()
    gquad_list = []
    eof = False

    gb = range(args.minG, args.maxG + 1)[::-1]
    gs = range(3, args.loops + 1)[::-1]
    features = defaultdict(list)
    longest = (args.maxG + args.maxloop) * args.loops + args.maxG

    while True:
        if not args.quiet:
            print(f"Processing {chrom}\n")
        while line.startswith(">") is False:
            ref_seq.append(line.strip())
            line = ref_seq_file.readline()
            if line == "":
                eof = True
                break
        ref_seq = "".join(ref_seq)
        ref_seq = ref_seq.upper().replace("U", "T")

        for g in gb:
            for s in gs:
                gstem_base = ""
                for i in range(g):
                    gstem_base += "G"
                reg = ""

                for i in range(s):
                    reg += "([gG]{%d}\w{%d,%d})" % (g,
                                                    args.minloop,
                                                    args.maxloop
                                                    )
                reg += "([gG]{%d})" % (g)

                ref_seq = _process_forward_seq(
                    chrom, features, gquad_list, reg, ref_seq, longest)
                if args.noreverse is False:
                    ref_seq = _process_reverse_seq(
                        chrom, features, gquad_list, reg, ref_seq, longest)

                _write_to_gff(gquad_list, args.gff_output)
                gquad_list = []

        if eof:
            break
        chrom = chrom_name(line)
        ref_seq = []
        line = ref_seq_file.readline()
        if line == "":
            break

    return features


def _process_forward_seq(chrom, features, gquad_list, reg, ref_seq, longest):
    for m in re.finditer(reg, ref_seq):
        seq = m.group(0)
        start = m.start()
        end = m.end()
        if len(ref_seq) > longest:
            ref = seq
        else:
            ref = ref_seq
        quad_id = chrom + "_" + str(m.start()) + "_" + str(m.end())
        gquad_list.append([chrom,
                           start,
                           end,
                           quad_id,
                           len(seq),
                           "+",
                           seq
                           ])
        if seq not in features["g4motif"]:
            features = update_dataFrame(features, reg, seq, ref)
            features["seq"].append(chrom)
        temp = ""
        for i in range(start, end):
            temp += "N"
        ref_seq = ref_seq[:start] + temp + ref_seq[end:]
    return ref_seq


def _process_reverse_seq(chrom, features, gquad_list, reg, ref_seq, longest):
    rev_ref_seq = revcomp(ref_seq)
    seq_len = len(ref_seq)
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
        gquad_list.append([chrom,
                           seq_len - end,
                           seq_len - start,
                           quad_id,
                           len(seq),
                           "-",
                           seq,
                           ])
        if seq not in features["g4motif"]:
            features = update_dataFrame(
                features, reg, seq, ref)
            features["seq"].append(chrom)
        temp = ""
        for i in range(start, end):
            temp += "N"
        rev_ref_seq = rev_ref_seq[:start] + \
            temp + rev_ref_seq[end:]
    return ref_seq


def _write_to_gff(gquad_list, dest):
    gquad_sorted = sort_table(gquad_list, (1, 2, 3))
    for xline in gquad_sorted:
        xline = "\t".join([str(x) for x in xline])
        with open(dest, "a") as out:
            out.write(xline + "\n")
