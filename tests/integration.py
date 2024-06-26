def inference_on_sample_fasta_files():
    import os
    import sys
    from G4Boost import main

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
