import re
from Bio import SeqIO
from collections import defaultdict
from lib.util import revcomp, update_dataFrame, sort_table, transcribe_sequence


def process_sequences(args):
    fasta_sequences = SeqIO.parse(open(args.fasta), 'fasta')

    features = defaultdict(list)

    for fasta in fasta_sequences:
        seq, chrom_name = fasta.seq, fasta.id

        if not args.quiet:
            print(f"Processing {chrom_name}\n")

        _process_sequence(args, seq, chrom_name, features)

    return features


def _process_sequence(args, seq, chrom_name, features):
    gb = range(args.minG, args.maxG + 1)[::-1]
    gs = range(3, args.loops + 1)[::-1]
    longest = (args.maxG + args.maxloop) * args.loops + args.maxG

    seq = transcribe_sequence(seq)

    gquad_list = []

    for g in gb:
        for s in gs:
            regex = _create_regex(g, s, args.minloop, args.maxloop)

            seq = _process_forward_seq(chrom_name,
                                       features,
                                       gquad_list,
                                       regex,
                                       seq,
                                       longest
                                       )

            if args.noreverse is False:
                seq = _process_reverse_seq(chrom_name,
                                           features,
                                           gquad_list,
                                           regex,
                                           seq,
                                           longest
                                           )

            _write_to_gff(gquad_list, args.gff_output)
            gquad_list = []


def _process_forward_seq(chrom, features, gquad_list, regex, ref_seq, longest):
    for m in re.finditer(regex, ref_seq):
        seq = m.group(0)
        start = m.start()
        end = m.end()
        if len(ref_seq) > longest:
            ref = seq
        else:
            ref = ref_seq
        quad_id = chrom + "_" + str(start) + "_" + str(end)
        gquad_list.append([chrom,
                           start,
                           end,
                           quad_id,
                           len(seq),
                           "+",
                           seq
                           ])
        if seq not in features["g4motif"]:
            features = update_dataFrame(features, regex, seq, ref)
            features["seq"].append(chrom)
        temp = ""
        for i in range(start, end):
            temp += "N"
        ref_seq = ref_seq[:start] + temp + ref_seq[end:]
    return ref_seq


def _process_reverse_seq(chrom, features, gquad_list, regex, ref_seq, longest):
    rev_ref_seq = revcomp(ref_seq)
    seq_len = len(ref_seq)
    for m in re.finditer(regex, rev_ref_seq):
        seq = m.group(0)
        start = m.start()
        end = m.end()
        if len(rev_ref_seq) > longest:
            ref = seq
        else:
            ref = rev_ref_seq
        quad_id = chrom + "_" + str(start) + "_" + str(end)
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
                features, regex, seq, ref)
            features["seq"].append(chrom)
        temp = ""
        for i in range(start, end):
            temp += "N"
        rev_ref_seq = rev_ref_seq[:start] + temp + rev_ref_seq[end:]
    return ref_seq


def _write_to_gff(gquad_list, dest):
    gquad_sorted = sort_table(gquad_list, (1, 2, 3))
    for xline in gquad_sorted:
        xline = "\t".join([str(x) for x in xline])
        with open(dest, "a") as out:
            out.write(xline + "\n")


def _create_regex(g, s, min_loop, max_loop):
    regex = ""

    for i in range(s):
        regex += "([gG]{%d}\w{%d,%d})" % (g, min_loop,
                                          max_loop
                                          )
    regex += "([gG]{%d})" % (g)

    return regex
