import re
import operator


def sort_table(table, cols):
    for col in reversed(cols):
        table = sorted(table, key=operator.itemgetter(col))
    return table


def chrom_name(header):
    if not header.startswith(">"):
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


def findmotifs(reg, seq, start, chrom):
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
