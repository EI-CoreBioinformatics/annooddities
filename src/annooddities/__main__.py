#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to find annotation oddities from GFF3 and genome FASTA files
"""

import os
import sys
import argparse
import csv
from collections import Counter, defaultdict
import logging
import shutil
import re
from pyfaidx import Fasta
import json
import yaml
import toml

script = os.path.basename(sys.argv[0])
executed_command = " ".join(sys.argv)


CANONICAL_INTRON_MOTIFS = {"GT-AG", "GC-AG", "AT-AC"}

METRIC_ODDITIES = [
    "{exon_num} == 1",
    "{exon_num} > 1",
    "{five_utr_length} > 10000",
    "{five_utr_num} > 5",
    "{three_utr_length} > 10000",
    "{three_utr_num} > 4",
    "not {is_complete}",
    "not {has_start_codon}",
    "not {has_stop_codon}",
    "{is_fragment}",
    "{has_inframe_stop}",
    "{max_exon_length} > 10000",
    "{max_intron_length} > 500000",
    "{min_exon_length} <= 5",
    "0 < {min_intron_length} <= 5",
    "{selected_cds_fraction} <= 0.3",
    "{canonical_intron_proportion} != 1",
    # TODO: lacking data to implement
    # "{non_verified_introns_num} >= 1",
    "{only_non_canonical_splicing}",
    # TODO: lacking data to implement
    # "{proportion_verified_introns} <= 0.5",
    "{suspicious_splicing}",
]

MIKADO_STATS_ORDER = [
    "Number of genes",
    "Number of Transcripts",
    "Transcripts per gene",
    "Number of monoexonic genes",
    "Monoexonic transcripts",
    "Transcript mean size cDNA (bp)",
    "Transcript median size cDNA (bp)",
    "Min cDNA",
    "Max cDNA",
    "Total exons",
    "Exons per transcript",
    "Exon mean size (bp)",
    "CDS mean size (bp)",
    "Transcript mean size CDS (bp)",
    "Transcript median size CDS (bp)",
    "Min CDS",
    "Max CDS",
    "Intron mean size (bp)",
    "5'UTR mean size (bp)",
    "3'UTR mean size (bp)",
]


class AnnoOddities:
    def __init__(self, args):
        self.args = args
        self.genome = os.path.abspath(self.args.genome_fasta)
        self.gff = os.path.abspath(self.args.gff3_file)
        self.force = self.args.force
        # check if input and output files/directories exist
        if not os.path.isfile(self.genome):
            logging.error(f"Reference genome file {self.genome} not found")
            sys.exit(1)
        if not os.path.isfile(self.gff):
            logging.error(f"GFF3 file {self.gff} not found")
            sys.exit(1)

        # use canonical intron motifs from args if provided, else use default
        global CANONICAL_INTRON_MOTIFS
        if self.args.canonical_intron_motifs:
            CANONICAL_INTRON_MOTIFS = set(self.args.canonical_intron_motifs.split(","))
            logging.info(
                f"Using user-defined canonical intron motifs: {CANONICAL_INTRON_MOTIFS}"
            )
        else:
            logging.info(
                f"Using default canonical intron motifs: {CANONICAL_INTRON_MOTIFS}"
            )

        # write all stats to all_stats.tsv in current directory
        self.all_stats_tsv = os.path.join(
            self.args.output_prefix + ".AnnoOddities.all_stats.tsv"
        )
        self.oddity_summary = os.path.join(
            self.args.output_prefix + ".AnnoOddities.oddity_summary.txt"
        )
        self.genome_name = os.path.basename(self.genome)
        self.gff_name = os.path.basename(self.gff)

        # genome dna index
        self.genome_db = str()

        self.sanitised_gff = str()
        self.agat_statistics_yaml = str()
        self.agat_statistics_json = str()
        self.agat_statistics_toml = str()
        self.gffread_table = str()
        self.mikado_tab_stats = str()
        self.mikado_stats = str()
        self.mikado_stats_yaml = str()
        self.mikado_stats_json = str()
        self.mikado_stats_toml = str()
        self.mikado_basic_stats = str()

        # for extraction of oddities from gff3, store tid, gid and oddities matched
        self.oddity_tid_info = defaultdict(dict)
        self.oddity_gid_info = defaultdict(list)

    @staticmethod
    def _format_number(x):
        try:
            int_x = int(x.replace(",", ""))
            return str(int_x)
        except ValueError:
            return "{:.02f}".format(float(x.replace(",", "")))

    @staticmethod
    def parse_mikado_stats(mikado_stats):
        """Parse Mikado stats to create basic summary stats"""

        output_stats_yaml = mikado_stats.replace(".tsv", ".yaml")
        output_stats_json = mikado_stats.replace(".tsv", ".json")
        output_stats_toml = mikado_stats.replace(".tsv", ".toml")

        output_stats = mikado_stats.replace(
            ".mikado_stats.tsv", ".mikado_summary_stats.tsv"
        )

        stats_info = defaultdict()
        data_list = []

        with open(mikado_stats) as f:
            reader = csv.DictReader(f, delimiter="\t")
            data_list = list(reader)
        for row in data_list:
            # Generate basic mikado stats info
            try:
                stat = row["Stat"].strip()
            except KeyError:
                continue
            if stat.lower() in {
                "number of genes",
                "transcripts per gene",
                "number of monoexonic genes",
                "monoexonic transcripts",
                "exons per transcript",
            }:
                stats_info[stat] = AnnoOddities._format_number(row["Total"])
                if stat.lower() == "transcripts per gene":
                    stats_info["Number of Transcripts"] = stats_info[stat]
                    stats_info[stat] = AnnoOddities._format_number(row["Average"])
                elif stat.lower() == "exons per transcript":
                    stats_info["Total exons"] = stats_info[stat]
                    stats_info[stat] = AnnoOddities._format_number(row["Average"])

            elif stat.lower() in {
                "exon lengths",
                "cds exon lengths",
                "intron lengths",
                "5'utr length",
                "3'utr length",
            }:
                stats_info["{} mean size (bp)".format(stat.split(" ")[0])] = (
                    AnnoOddities._format_number(row["Average"])
                )

            elif stat.lower() in {"cdna lengths", "cds lengths"}:
                stats_info[
                    "Transcript mean size {} (bp)".format(stat.split(" ")[0])
                ] = AnnoOddities._format_number(row["Average"])
                stats_info[
                    "Transcript median size {} (bp)".format(stat.split(" ")[0])
                ] = AnnoOddities._format_number(row["Median"])
                stats_info["Min {}".format(stat.split(" ")[0])] = (
                    AnnoOddities._format_number(row["Min"])
                )
                stats_info["Max {}".format(stat.split(" ")[0])] = (
                    AnnoOddities._format_number(row["Max"])
                )
        with open(output_stats, "w") as out_stats:
            for stat in MIKADO_STATS_ORDER:
                print(
                    stat,
                    stats_info[stat.replace("cDNA", "CDNA")],
                    sep="\t",
                    file=out_stats,
                )
        data = defaultdict(dict)
        for row in data_list:
            stat_name = row["Stat"]
            data["Stat"][stat_name] = {
                "Total": row["Total"],
                "Average": row["Average"],
                "Mode": row["Mode"],
                "Min": row["Min"],
                "1%": row["1%"],
                "5%": row["5%"],
                "10%": row["10%"],
                "25%": row["25%"],
                "Median": row["Median"],
                "75%": row["75%"],
                "90%": row["90%"],
                "95%": row["95%"],
                "99%": row["99%"],
                "Max": row["Max"],
            }
        with open(output_stats_yaml, "w") as yaml_file:
            yaml.dump(dict(data), yaml_file)

        with open(output_stats_json, "w") as json_file:
            json.dump(data, json_file, indent=4)

        with open(output_stats_toml, "w") as toml_file:
            toml.dump(dict(data), toml_file)

        return output_stats_yaml, output_stats_json, output_stats_toml, output_stats

    @staticmethod
    def parse_mikado_tab_stats(mikado_tab_stats):
        """Parse Mikado tab-stats file and return dict of transcript metrics"""

        transcripts = defaultdict(dict)

        with open(mikado_tab_stats) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                # Generate per transcript info
                tid = row["TID"]
                transcripts[tid] = {
                    # moved to gffread table parsing
                    # "exon_num": int(row["Exon number"]),
                    # "total_exon_length": int(row["cDNA length"]),
                    # patch fix as +1 was not added when calculating intron length in Mikado. Add Exon number to account for this
                    # "total_intron_length": int(row["Intronic length"]),
                    # "total_cds_length": int(row["CDS length"]),
                    # "cds_exon_num": int(row["# coding exons"]),
                    # convert from percentage to float for easier calculations
                    # "cds_cdna_ratio": float(row["cDNA/CDS ratio"]) / 100.0,
                    "five_utr_length": int(row["5'UTR"]),
                    "five_utr_num": int(row["5'UTR exons"]),
                    "three_utr_length": int(row["3'UTR"]),
                    "three_utr_num": int(row["3'UTR exons"]),
                }

        return transcripts

    @staticmethod
    def _get_coords(tid, coord_type, coords_str):
        """
        Convert exons string to list of (start, end) tuples.
        Exons are given as start-end,start-end
        param coords_str: str, exons or cds string from gffread table
        """
        exons = []
        try:
            int(coords_str.split(",")[0].split("-")[0])
        except ValueError:
            logging.debug(
                f"Warning: Unable to parse {coord_type} coordinates for transcript {tid}. The coordinates are: '{coords_str}'"
            )
            return []
        for exon in coords_str.split(","):
            start, end = map(int, exon.split("-"))
            if start > end:
                start, end = end, start
            exons.append((start, end))
        return exons

    @staticmethod
    def _get_introns(exons_list):
        """
        Convert list of (start, end) tuples to list of intron (start, end) tuples.
        param exons_list: list of (start, end) tuples
        """
        introns = []
        for i in range(1, len(exons_list)):
            intron_start = exons_list[i - 1][1] + 1
            intron_end = exons_list[i][0] - 1
            introns.append((intron_start, intron_end))
        return introns

    @staticmethod
    def _get_cumulative_length(lengths_list):
        """
        Convert list of (start, end) tuples to cumulative length.
        param lengths_list: list of (start, end) tuples
        """
        return sum(end - start + 1 for start, end in lengths_list)

    @staticmethod
    def _get_lengths(lengths_list):
        """
        Convert list of (start, end) tuples to list of lengths.
        param lengths_list: list of (start, end) tuples
        """
        return [end - start + 1 for start, end in lengths_list]

    @staticmethod
    def _get_max_length(lengths_list):
        """
        Get maximum length from list of (start, end) tuples.
        param lengths_list: list of (start, end) tuples
        """
        return (
            max(end - start + 1 for start, end in lengths_list) if lengths_list else 0
        )

    @staticmethod
    def _get_min_length(lengths_list):
        """
        Get minimum length from list of (start, end) tuples.
        param lengths_list: list of (start, end) tuples
        """
        return (
            min(end - start + 1 for start, end in lengths_list) if lengths_list else 0
        )

    @staticmethod
    def _reverse_complement(seq):
        """Get reverse complement of a sequence"""
        complement = str.maketrans("ACGTacgt", "TGCAtgca")
        return seq.translate(complement)[::-1]

    @staticmethod
    def _get_junction_classification(genome_db, chrom, strand, introns):
        """Get intron motif from genome fasta"""
        # create defalultdict list to store junctions
        # GT-AG:++:count;percentage;
        # GC-AG:++:count;percentage;
        # AT-AC:++:count;percentage;
        # GT-AG:+-:count;percentage;
        # GC-AG:+-:count;percentage;
        # AT-AC:+-:count;percentage;
        # GT-AG;-+;count;percentage;
        # GC-AG;-+;count;percentage;
        # AT-AC;-+;count;percentage;
        # GT-AG:--:count;percentage;
        # GC-AG:--:count;percentage;
        # AT-AC:--:count;percentage;
        known_strand_junctions = defaultdict(Counter)
        unknown_strand_junctions = defaultdict(Counter)
        # logic to handle strand to be uniformly +, - or .
        strand = strand if strand in ["+", "-"] else "."
        for intron in introns:
            intron_start, intron_end = intron
            donor_start = intron_start - 1
            donor_end = intron_start + 1
            acceptor_start = intron_end - 2
            acceptor_end = intron_end

            # DEBUG
            # print(
            #     "Intron:",
            #     chrom,
            #     intron_start,
            #     intron_start + 1,
            #     intron_end - 1,
            #     intron_end,
            #     strand,
            #     donor_start,
            #     donor_end,
            #     acceptor_start,
            #     acceptor_end,
            #     genome_db[chrom][donor_start:donor_end].seq,
            #     genome_db[chrom][acceptor_start:acceptor_end].seq,
            #     -genome_db[chrom][donor_start:donor_end],
            #     -genome_db[chrom][acceptor_start:acceptor_end],
            # )

            # my notes
            # :exp:++:	# canonical intron expected strand
            # :exp:--:	# canonical intron expected strand
            # :sus:+-:	# suspicious intron, if guilty!
            # :sus:-+:	# suspicious intron, if guilty!
            # :und:.+:	# unknown strand, could be canonical intron
            # :und:.-:	# unknown strand, could be canonical intron

            # https://github.com/mdshw5/pyfaidx
            donor_seq = genome_db[chrom][donor_start:donor_end].seq.upper()
            acceptor_seq = genome_db[chrom][acceptor_start:acceptor_end].seq.upper()
            junction = donor_seq + "-" + acceptor_seq

            if strand == "+":
                # :++:
                # since I am reading +ve strand, i can use the current junction as is
                known_strand_junctions[junction][":exp:++:"] += 1
                # :+-:
                junction = AnnoOddities._reverse_complement(junction)
                known_strand_junctions[junction][":sus:+-:"] += 1
            elif strand == "-":
                # :--:
                junction = AnnoOddities._reverse_complement(junction)
                known_strand_junctions[junction][":exp:--:"] += 1
                junction = AnnoOddities._reverse_complement(junction)
                known_strand_junctions[junction][":sus:-+:"] += 1
            else:
                # :.+:
                unknown_strand_junctions[junction][":und:.+:"] += 1
                # :.-:
                # none exists as we do not know the strand
                # :.-:
                junction = AnnoOddities._reverse_complement(junction)
                unknown_strand_junctions[junction][":und:.-:"] += 1
                # :.+:
                # none exists as we do not know the strand

        # Check each transcript junctions dict for canonical, suspicious and non-canonical introns

        (
            canonical_count,
            non_canonical_count,
            suspicious_count,
            unknown_predicted_strand,
        ) = AnnoOddities.classify_junctions(
            known_strand_junctions, junction_type="known"
        )
        (
            unknown_canonical_count,
            unknown_non_canonical_count,
            unknown_suspicious_count,
            unknown_predicted_strand,
        ) = AnnoOddities.classify_junctions(
            unknown_strand_junctions, junction_type="unknown"
        )

        # count of intron matching canonical intron then canonical proportion is 1
        canonical_intron_proportion = 0.0
        only_non_canonical_splicing = True
        suspicious_splicing = False
        if len(introns) == 0:
            canonical_intron_proportion = 1.0
            only_non_canonical_splicing = False
            suspicious_splicing = False
        else:
            if strand in ["+", "-"]:
                canonical_intron_proportion = float(canonical_count / len(introns))
                only_non_canonical_splicing = canonical_count == 0 and len(introns) > 0
                suspicious_splicing = suspicious_count > 0
            else:
                canonical_intron_proportion = float(
                    unknown_canonical_count / len(introns)
                )
                only_non_canonical_splicing = (
                    unknown_canonical_count == 0 and len(introns) > 0
                )
                suspicious_splicing = unknown_suspicious_count > 0

        # DEBUG
        # print(f"Known Junctions dict: {known_strand_junctions}")
        # print(f"Unknown Junctions dict: {unknown_strand_junctions}")
        # print(
        #     f"Final counts - Canonical introns: {canonical_count}, Non-canonical introns: {non_canonical_count}, Suspicious introns: {suspicious_count}, Unknown strand Canonical introns: {unknown_canonical_count}, Unknown strand Non-canonical introns: {unknown_non_canonical_count}, Unknown strand Suspicious introns: {unknown_suspicious_count}, Unknown predicted strand: {unknown_predicted_strand}"
        # )
        # print(
        #     f"Canonical intron proportion: {canonical_intron_proportion}, Only non-canonical splicing: {only_non_canonical_splicing}, Suspicious splicing: {suspicious_splicing}"
        # )

        return (
            canonical_intron_proportion,
            only_non_canonical_splicing,
            suspicious_splicing,
            known_strand_junctions,
            unknown_strand_junctions,
            canonical_count,
            non_canonical_count,
            suspicious_count,
            unknown_canonical_count,
            unknown_non_canonical_count,
            unknown_suspicious_count,
            unknown_predicted_strand,
        )

    @staticmethod
    def classify_junctions(junctions, junction_type=None):
        """
        Classify junctions into canonical, non-canonical and suspicious.
        At this point, we have junctions and their count
        param junctions: dict of junctions
        """
        canonical_count = 0
        non_canonical_count = 0
        suspicious_count = 0
        unknown_positive_junction_found = False
        unknown_negative_junction_found = False
        unknown_positive_canonical_count = 0
        unknown_negative_canonical_count = 0
        unknown_predicted_strand = None
        for junction, counts in junctions.items():
            for strand_sign, count in counts.items():
                # DEBUG
                # print(
                #     f"Junction: {junction}, Strand sign: {strand_sign}, Count: {count}"
                # )
                # expected canonical intron are :exp:++: and :exp:--: only these count as canonical introns
                if (
                    strand_sign in [":exp:++:", ":exp:--:"]
                    and junction in CANONICAL_INTRON_MOTIFS
                ):
                    canonical_count += count
                # suspicious introns are all :sus:+-: and :sus:-+: entries but only if they match canonical intron motifs
                elif (
                    strand_sign in [":sus:+-:", ":sus:-+:"]
                    and junction in CANONICAL_INTRON_MOTIFS
                ):
                    # DEBUG
                    # print(
                    #     f"Suspicious intron found: {junction} with strand sign {strand_sign}"
                    # )
                    suspicious_count += count
                # non-canonical introns are the ones not matching canonical motifs but in expected strand
                else:
                    if (junction_type == "known") and (
                        strand_sign in [":exp:++:", ":exp:--:"]
                        and junction not in CANONICAL_INTRON_MOTIFS
                    ):
                        non_canonical_count += count
                    # elif (junction_type == "unknown") and (
                    #     strand_sign in [":und:.+:", ":und:.-:"]
                    #     and junction not in CANONICAL_INTRON_MOTIFS
                    # ):
                    #     non_canonical_count += count
                    # else:
                    #     continue
                    # for unknown junctions, if one juction is in ":und:.+: and also in CANONICAL_INTRON_MOTIFS and other is in :und:.-:" and also in CANONICAL_INTRON_MOTIFS then count as suspicious. This is because we cannot be sure of the strand, so if both junctions are canonical, then it is suspicious. And if all the junctions belongs to only ":und:.+: or :und:.-:" and also in CANONICAL_INTRON_MOTIFS then it is canonical. All other junctions are non-canonical.

                    # logic: if for a transcript with 10 introns, 4 on the positive strand and 6 on the negative strand, and all 10 junctions are canonical, then canonical intron count is based on the higher count of canonical junctions on either strand, i.e., 6 in this case or 0.6 proportion and suspicious. If for a transcript with 10 introns, 4 on the positive strand and 6 on the negative strand, and only 2 junctions are canonical (1 on each strand), then canonical intron count is 1 (higher of the two), proportion is 0.1 and suspicious.
                    elif junction_type == "unknown":
                        if (
                            strand_sign in [":und:.+:"]
                            and junction in CANONICAL_INTRON_MOTIFS
                        ):
                            unknown_positive_canonical_count += count
                            unknown_positive_junction_found = True
                        if (
                            strand_sign in [":und:.-:"]
                            and junction in CANONICAL_INTRON_MOTIFS
                        ):
                            unknown_negative_canonical_count += count
                            unknown_negative_junction_found = True

                        if all(
                            [
                                unknown_positive_junction_found,
                                unknown_negative_junction_found,
                            ]
                        ):
                            suspicious_count += count
                            unknown_predicted_strand = "?"
                        else:
                            if strand_sign in [":und:.+:"]:
                                unknown_predicted_strand = "+"
                            elif strand_sign in [":und:.-:"]:
                                unknown_predicted_strand = "-"

                        if (
                            strand_sign in [":und:.+:", ":und:.-:"]
                            and junction not in CANONICAL_INTRON_MOTIFS
                        ):
                            non_canonical_count += count

        # for unknown junction_type, the canonical_count is the highest of the unknown_positive_canonical_count and unknown_negative_canonical_count
        if junction_type == "unknown":
            canonical_count = max(
                unknown_positive_canonical_count, unknown_negative_canonical_count
            )

        return (
            canonical_count,
            non_canonical_count,
            suspicious_count,
            unknown_predicted_strand,
        )

    @staticmethod
    def parse_gffread_table(genome_db, gffread_file):
        """Parse gffread table file and return dict of transcript metrics"""
        transcripts = defaultdict(dict)
        # there is no header line in gffread --table output, so we define the fieldnames
        fieldnames = [
            "TID",
            "GID",
            "chr",
            "start",
            "end",
            "strand",
            "numexons",
            "exons",
            "cds",
            "covlen",
            "cdslen",
            "InFrameStop",
            "partialness",
        ]
        model_categories = {
            ".": "complete",
            "5_3": "fragment",
            "3": "3prime_partial",
            "5": "5prime_partial",
        }
        model_inframestop = {".": "no_inframe_stop", "true": "has_inframe_stop"}
        with open(gffread_file) as f:
            reader = csv.DictReader(f, delimiter="\t", fieldnames=fieldnames)
            for row in reader:
                tid = row["TID"]
                transcripts[tid] = {
                    "has_start_codon": model_categories.get(row["partialness"], "NA")
                    not in ["5prime_partial", "fragment"],
                    "has_stop_codon": model_categories.get(row["partialness"], "NA")
                    not in ["3prime_partial", "fragment"],
                    "is_complete": model_categories.get(row["partialness"], "NA")
                    == "complete",
                    "is_fragment": model_categories.get(row["partialness"], "NA")
                    == "fragment",
                    "has_inframe_stop": model_inframestop.get(row["InFrameStop"], "NA")
                    == "has_inframe_stop",
                }
                transcripts[tid]["gene_id"] = str(row["GID"])
                transcripts[tid]["chrom"] = str(row["chr"])
                transcripts[tid]["start"] = int(row["start"])
                transcripts[tid]["end"] = int(row["end"])
                # region chr:start..end
                transcripts[tid]["region"] = (
                    f"{row['chr']}:{row['start']}..{row['end']}"
                )
                transcripts[tid]["strand"] = str(row["strand"])
                transcripts[tid]["exon_num"] = int(row["numexons"])
                # lets work out the exon and intron lengths from the exons field
                transcripts[tid]["exons"] = AnnoOddities._get_coords(
                    tid, "exon", row["exons"]
                )
                transcripts[tid]["exon_lengths"] = AnnoOddities._get_lengths(
                    transcripts[tid]["exons"]
                )
                transcripts[tid]["total_exon_length"] = row["covlen"]
                transcripts[tid]["max_exon_length"] = AnnoOddities._get_max_length(
                    transcripts[tid]["exons"]
                )
                transcripts[tid]["min_exon_length"] = AnnoOddities._get_min_length(
                    transcripts[tid]["exons"]
                )
                # CDS length
                transcripts[tid]["cds_exons"] = AnnoOddities._get_coords(
                    tid, "cds", row["cds"]
                )
                # cds_exon_num
                transcripts[tid]["cds_exon_num"] = len(transcripts[tid]["cds_exons"])
                transcripts[tid]["cds_exon_lengths"] = AnnoOddities._get_lengths(
                    transcripts[tid]["cds_exons"]
                )
                transcripts[tid]["total_cds_length"] = (
                    0 if row["cdslen"] == "." else int(row["cdslen"])
                )
                transcripts[tid]["max_cds_length"] = AnnoOddities._get_max_length(
                    transcripts[tid]["cds_exons"]
                )
                transcripts[tid]["min_cds_length"] = AnnoOddities._get_min_length(
                    transcripts[tid]["cds_exons"]
                )

                # calculate CDS/cDNA ratio as float
                if int(row["cdslen"]) > 0 and int(row["covlen"]) > 0:
                    transcripts[tid]["cds_cdna_ratio"] = float(row["cdslen"]) / float(
                        row["covlen"]
                    )
                else:
                    transcripts[tid]["cds_cdna_ratio"] = 0.0

                # get intron coordinates and calculate lengths
                transcripts[tid]["introns"] = AnnoOddities._get_introns(
                    transcripts[tid]["exons"]
                )
                # intron_num is len of introns
                transcripts[tid]["intron_num"] = len(transcripts[tid]["introns"])
                transcripts[tid]["intron_lengths"] = AnnoOddities._get_lengths(
                    transcripts[tid]["introns"]
                )
                transcripts[tid]["total_intron_length"] = (
                    AnnoOddities._get_cumulative_length(transcripts[tid]["introns"])
                )
                transcripts[tid]["max_intron_length"] = AnnoOddities._get_max_length(
                    transcripts[tid]["introns"]
                )
                transcripts[tid]["min_intron_length"] = AnnoOddities._get_min_length(
                    transcripts[tid]["introns"]
                )

                (
                    transcripts[tid]["canonical_intron_proportion"],
                    transcripts[tid]["only_non_canonical_splicing"],
                    transcripts[tid]["suspicious_splicing"],
                    transcripts[tid]["known_strand_junctions"],
                    transcripts[tid]["unknown_strand_junctions"],
                    transcripts[tid]["canonical_intron_count"],
                    transcripts[tid]["non_canonical_intron_count"],
                    transcripts[tid]["suspicious_intron_count"],
                    transcripts[tid]["unknown_canonical_intron_count"],
                    transcripts[tid]["unknown_non_canonical_intron_count"],
                    transcripts[tid]["unknown_suspicious_intron_count"],
                    transcripts[tid]["unknown_predicted_strand"],
                ) = AnnoOddities._get_junction_classification(
                    genome_db=genome_db,
                    chrom=row["chr"],
                    strand=row["strand"],
                    introns=transcripts[tid]["introns"],
                )

        return transcripts

    @staticmethod
    def _join_splice_site_dict(splice_site_dict):
        """
        Join splice site dict to string
        {'GT-AG': {':exp:--:': 6}, 'CT-AC': {':sus:-+:': 6}}
        to
        'GT-AG:exp:--:6:CT-AC:sus:-+:6'
        """
        return (
            ";".join(
                f"{splice_site}{inner_key}{inner_count}"
                # 1. Outer iteration (for splice_site in splice_site_dict)
                for splice_site, counter_obj in splice_site_dict.items()
                # 2. Inner iteration (for inner_key in counter_obj)
                for inner_key, inner_count in counter_obj.items()
            )
            if splice_site_dict
            else "NA"
        )

    headers = [
        "transcript_id",
        "gene_id",
        "chromosome",
        "start",
        "end",
        "region",
        "strand",
        "exon_num",
        "exons",
        "exon_lengths",
        "max_exon_length",
        "min_exon_length",
        "total_exon_length",
        "cds_exon_num",
        "cds_exons",
        "cds_exon_lengths",
        "max_cds_length",
        "min_cds_length",
        "total_cds_length",
        "cds_cdna_ratio",
        "intron_num",
        "introns",
        "intron_lengths",
        "max_intron_length",
        "min_intron_length",
        "total_intron_length",
        "five_utr_num",
        "five_utr_length",
        "three_utr_num",
        "three_utr_length",
        "has_start_codon",
        "has_stop_codon",
        "is_complete",
        "is_fragment",
        "has_inframe_stop",
        "known_strand_junctions_dict",
        "known_strand_junctions_str",
        "unknown_strand_junctions_dict",
        "unknown_strand_junctions_str",
        "canonical_intron_count",
        "non_canonical_intron_count",
        "suspicious_intron_count",
        "unknown_canonical_intron_count",
        "unknown_non_canonical_intron_count",
        "unknown_suspicious_intron_count",
        "unknown_predicted_strand",
        "canonical_intron_proportion",
        "only_non_canonical_splicing",
        "suspicious_splicing",
        "matched_oddities",
    ]

    @staticmethod
    def calculate_oddities(mikado_data, gffread_data, oddities, all_stats_tsv):
        """Calculate oddities for all transcripts"""
        oddity_counts = Counter({oddity: 0 for oddity in oddities})
        oddity_tid_info = defaultdict(dict)
        oddity_gid_info = defaultdict(list)

        with open(all_stats_tsv, "w") as stats_out:
            stats_out.write("\t".join(AnnoOddities.headers) + "\n")
            # using mikado as the base, as it should have all transcripts we care about
            for tid, mikado_metrics in mikado_data.items():
                gffread_metrics = gffread_data.get(tid, {})

                combined_metrics = {**mikado_metrics, **gffread_metrics}
                # DEBUG
                # print(f"Processing transcript: {tid}")
                # print(f"Mikado metrics: {mikado_metrics}")
                # print(f"gffread metrics: {gffread_metrics}")
                # # Combine metrics
                # print(f"Combined metrics before derived for {tid}: {combined_metrics}")

                # Add derived metrics - for oddity calculations
                combined_metrics["selected_cds_fraction"] = combined_metrics[
                    "cds_cdna_ratio"
                ]
                # DEBUG
                # print(f"Processing transcript: {tid} with metrics: {combined_metrics}")

                # Check each oddity condition
                matched_oddities = []
                matched_oddities_tsv = []

                for oddity in oddities:
                    try:
                        if eval(oddity.format(**combined_metrics)):
                            oddity_counts[oddity] += 1
                            metric_name = AnnoOddities._clean_metric_name(oddity)
                            # remove exon_num_gt_1 and exon_num_eq_1 from matched oddities tsv as they are not real oddities but retain in the main GFF
                            if metric_name not in ["exon_num_gt_1", "exon_num_eq_1"]:
                                matched_oddities_tsv.append(metric_name)
                            matched_oddities.append(metric_name)
                            oddity_tid_info[tid]["gene_id"] = combined_metrics.get(
                                "gene_id", "NA"
                            )
                            oddity_tid_info[tid]["matched_oddities"] = ";".join(
                                matched_oddities
                            )
                            oddity_gid_info[
                                combined_metrics.get("gene_id", "NA")
                            ].append(metric_name)
                    except Exception as e:
                        logging.error(f"Error evaluating oddity '{oddity}': {e}")
                # if matched_oddities_tsv is empty, set to None
                if not matched_oddities_tsv:
                    matched_oddities_tsv.append("None")
                combined_metrics["matched_oddities"] = ";".join(matched_oddities_tsv)

                # DEBUG
                # print(f"Combined metrics after derived for {tid}: {combined_metrics}")

                # Write combined metrics to file
                stats_out.write(
                    "\t".join(
                        [
                            tid,
                            # gene id
                            combined_metrics.get("gene_id", "NA"),
                            # chromosome
                            combined_metrics.get("chrom", "NA"),
                            # start
                            str(combined_metrics.get("start", "NA")),
                            # end
                            str(combined_metrics.get("end", "NA")),
                            # region
                            combined_metrics.get("region", "NA"),
                            # strand
                            str(combined_metrics.get("strand", "NA")),
                            # exon number
                            str(combined_metrics.get("exon_num", "NA")),
                            # exon coordinates
                            str(combined_metrics.get("exons", "NA")),
                            # exon coordinates lengths
                            str(combined_metrics.get("exon_lengths", "NA")),
                            # max exon coordinate length
                            str(combined_metrics.get("max_exon_length", "NA")),
                            # min exon coordinate length
                            str(combined_metrics.get("min_exon_length", "NA")),
                            # total exon length
                            str(combined_metrics.get("total_exon_length", "NA")),
                            # cds exon number
                            str(combined_metrics.get("cds_exon_num", "NA")),
                            # cds coordinates
                            str(combined_metrics.get("cds_exons", "NA")),
                            # cds coordinates lengths
                            str(combined_metrics.get("cds_exon_lengths", "NA")),
                            # max cds coordinate length
                            str(combined_metrics.get("max_cds_length", "NA")),
                            # min cds coordinate length
                            str(combined_metrics.get("min_cds_length", "NA")),
                            # total cds length
                            str(combined_metrics.get("total_cds_length", "NA")),
                            # cds/cdna ratio
                            str(combined_metrics.get("cds_cdna_ratio", "NA")),
                            # intron number
                            str(combined_metrics.get("intron_num", "NA")),
                            # introns
                            str(combined_metrics.get("introns", "NA")),
                            # intron_lengths
                            str(combined_metrics.get("intron_lengths", "NA")),
                            # max_intron_length
                            str(combined_metrics.get("max_intron_length", "NA")),
                            # min_intron_length
                            str(combined_metrics.get("min_intron_length", "NA")),
                            # total intron length
                            str(combined_metrics.get("total_intron_length", "NA")),
                            # five_utr_num
                            str(combined_metrics.get("five_utr_num", "NA")),
                            # five_utr_length
                            str(combined_metrics.get("five_utr_length", "NA")),
                            # three_utr_num
                            str(combined_metrics.get("three_utr_num", "NA")),
                            # three_utr_length
                            str(combined_metrics.get("three_utr_length", "NA")),
                            # has_start_codon
                            str(combined_metrics.get("has_start_codon", "NA")),
                            # has_stop_codon
                            str(combined_metrics.get("has_stop_codon", "NA")),
                            # is_complete
                            str(combined_metrics.get("is_complete", "NA")),
                            # is_fragment
                            str(combined_metrics.get("is_fragment", "NA")),
                            # has_inframe_stop
                            str(combined_metrics.get("has_inframe_stop", "NA")),
                            # known_strand_junctions
                            # known_strand_junctions_dict
                            str(
                                {
                                    k: dict(v)
                                    for k, v in combined_metrics.get(
                                        "known_strand_junctions", {}
                                    ).items()
                                }
                            ),
                            # known_strand_junctions_str
                            AnnoOddities._join_splice_site_dict(
                                combined_metrics.get("known_strand_junctions", {})
                            ),
                            # unknown_strand_junctions
                            # unknown_strand_junctions_dict
                            str(
                                {
                                    k: dict(v)
                                    for k, v in combined_metrics.get(
                                        "unknown_strand_junctions", {}
                                    ).items()
                                }
                            ),
                            # unknown_strand_junctions_str
                            AnnoOddities._join_splice_site_dict(
                                combined_metrics.get("unknown_strand_junctions", {})
                            ),
                            # canonical_intron_count
                            str(combined_metrics.get("canonical_intron_count", "NA")),
                            # non_canonical_intron_count
                            str(
                                combined_metrics.get("non_canonical_intron_count", "NA")
                            ),
                            # suspicious_intron_count
                            str(
                                combined_metrics.get("suspicious_intron_count", "NA"),
                            ),
                            # unknown_canonical_intron_count
                            str(
                                combined_metrics.get(
                                    "unknown_canonical_intron_count", "NA"
                                )
                            ),
                            # unknown_non_canonical_intron_count
                            str(
                                combined_metrics.get(
                                    "unknown_non_canonical_intron_count", "NA"
                                ),
                            ),
                            # unknown_suspicious_intron_count
                            str(
                                combined_metrics.get(
                                    "unknown_suspicious_intron_count", "NA"
                                ),
                            ),
                            # unknown_predicted_strand
                            str(combined_metrics.get("unknown_predicted_strand", "NA")),
                            # canonical_intron_proportion
                            str(
                                combined_metrics.get(
                                    "canonical_intron_proportion", "NA"
                                )
                            ),
                            # only_non_canonical_splicing
                            str(
                                combined_metrics.get(
                                    "only_non_canonical_splicing", "NA"
                                )
                            ),
                            # suspicious_splicing
                            str(combined_metrics.get("suspicious_splicing", "NA")),
                            # matched_oddities
                            str(combined_metrics.get("matched_oddities", "NA")),
                        ]
                    )
                )
                stats_out.write("\n")

        return oddity_counts, oddity_tid_info, oddity_gid_info

    @staticmethod
    def _clean_metric_name(oddity):
        return (
            oddity.replace("{", "")
            .replace("}", "")
            .replace("==", "eq")
            .replace(">", "gt")
            .replace("<=", "lte")
            .replace("<", "lt")
            .replace(">=", "gte")
            .replace("!=", "neq")
            .replace(" ", "_")
        )

    def run_agat_standardise_gff(self):
        """Run AGAT to standardise GFF3 file"""
        logging.info("Standardising GFF3 file using AGAT")
        output_gff = self.args.output_prefix + ".agat_standardised.gff"
        output_log = output_gff.replace(".gff", ".log")

        if os.path.isfile(output_gff) and not self.force:
            logging.info(
                f"Skipping agat standardisation, '{output_gff}' already exists. Use --force to rerun."
            )
            self.sanitised_gff = output_gff
            return

        cmd = f"agat_convert_sp_gxf2gxf.pl --gff {self.gff} --output {output_gff} > {output_log} 2>&1"
        self.process_cmd(cmd)
        self.sanitised_gff = output_gff

    def run_agat_sp_statistics(self):
        """Run AGAT to generate statistics from GFF3 file"""
        logging.info("Generating GFF3 statistics using AGAT")
        output_stats = self.sanitised_gff.replace(
            ".agat_standardised.gff", ".agat_sp_statistics.txt"
        )
        output_stats_yaml = output_stats.replace(".txt", ".yaml")
        output_stats_json = output_stats.replace(".txt", ".json")
        output_stats_toml = output_stats.replace(".txt", ".toml")
        output_log = output_stats.replace(".txt", ".log")

        if (
            os.path.isfile(output_stats)
            and os.path.isfile(output_stats_yaml)
            and os.path.isfile(output_stats_json)
            and os.path.isfile(output_stats_toml)
            and not self.force
        ):
            logging.info(
                f"Skipping agat statistics, '{output_stats}', '{output_stats_yaml}', '{output_stats_json}' and '{output_stats_toml}' already exist. Use --force to rerun."
            )
            self.agat_statistics_yaml = output_stats_yaml
            self.agat_statistics_json = output_stats_json
            self.agat_statistics_toml = output_stats_toml
            return

        cmd = (
            "agat_sp_statistics.pl "
            f"-g {self.genome} "
            f"--gff {self.sanitised_gff} "
            f"--yaml "
            f"--output {output_stats.replace('.txt', '')} > {output_log} 2>&1 "
            f" && mv {output_stats.replace('.txt', '')} {output_stats}"
        )
        self.process_cmd(cmd)

        # convert yaml to json and toml
        with open(output_stats_yaml) as yaml_file:
            yaml_data = yaml.safe_load(yaml_file)

        with open(output_stats_json, "w") as json_file:
            json.dump(yaml_data, json_file, indent=4)
        with open(output_stats_toml, "w") as toml_file:
            toml.dump(yaml_data, toml_file)
        self.agat_statistics_json = output_stats_json
        self.agat_statistics_yaml = output_stats_yaml
        self.agat_statistics_toml = output_stats_toml

    def run_gffread_table(self):
        """Run gffread to generate table from GFF3 file"""
        logging.info("Generating GFF3 table using gffread")
        output_table = self.sanitised_gff.replace(
            ".agat_standardised.gff", ".gffread_table.tbl"
        )
        output_log = output_table.replace(".tbl", ".log")

        if os.path.isfile(output_table) and not self.force:
            logging.info(
                f"Skipping gffread table generation, '{output_table}' already exists. Use --force to rerun."
            )
            self.gffread_table = output_table
            return

        cmd = (
            "gffread "
            f"{self.sanitised_gff} "
            f"-g {self.genome} "
            "-P "
            "--table @id,@geneid,@chr,@start,@end,@strand,@numexons,@exons,@cds,@covlen,@cdslen,InFrameStop,partialness "
            f" -o {output_table} > {output_log} 2>&1"
        )
        self.process_cmd(cmd)
        self.gffread_table = output_table

    def run_mikado_tab_stats(self):
        """Run Mikado tab-stats to generate metrics from GFF3 file"""
        logging.info("Generating Mikado tab-stats from GFF3 file")
        output_tab_stats = self.sanitised_gff.replace(
            ".agat_standardised.gff", ".mikado_tab_stats.tsv"
        )
        output_stats = output_tab_stats.replace(
            ".mikado_tab_stats.tsv", ".mikado_stats.tsv"
        )
        output_log = output_tab_stats.replace(".tsv", ".log")
        output_basic_stats = output_stats.replace(
            ".mikado_stats.tsv", ".mikado_summary_stats.tsv"
        )
        output_stats_yaml = output_stats.replace(".tsv", ".yaml")
        output_stats_json = output_stats.replace(".tsv", ".json")
        output_stats_toml = output_stats.replace(".tsv", ".toml")

        if (
            os.path.isfile(output_tab_stats)
            and os.path.isfile(output_stats)
            and os.path.isfile(output_basic_stats)
            and os.path.isfile(output_stats_yaml)
            and os.path.isfile(output_stats_json)
            and os.path.isfile(output_stats_toml)
            and not self.force
        ):
            logging.info(
                f"Skipping Mikado tab-stats generation, '{output_tab_stats}', '{output_stats}', '{output_basic_stats}', '{output_stats_yaml}', '{output_stats_json}' and '{output_stats_toml}' already exist. Use --force to rerun."
            )
            self.mikado_tab_stats = output_tab_stats
            self.mikado_stats = output_stats
            self.mikado_basic_stats = output_basic_stats
            self.mikado_stats_yaml = output_stats_yaml
            self.mikado_stats_json = output_stats_json
            self.mikado_stats_toml = output_stats_toml
            return

        cmd = (
            "mikado util stats "
            f"{self.sanitised_gff} "
            f"--tab-stats {output_tab_stats} > {output_stats} 2> {output_log}"
        )
        self.process_cmd(cmd)

        output_stats_yaml, output_stats_json, output_stats_toml, output_basic_stats = (
            self.parse_mikado_stats(output_stats)
        )
        self.mikado_tab_stats = output_tab_stats
        self.mikado_stats = output_stats
        self.mikado_stats_yaml = output_stats_yaml
        self.mikado_stats_json = output_stats_json
        self.mikado_stats_toml = output_stats_toml
        self.mikado_basic_stats = output_basic_stats
        logging.info(f"Mikado tab-stats written to {self.mikado_tab_stats}")
        logging.info(f"Mikado stats written to {self.mikado_stats}")
        logging.info(f"Mikado stats YAML written to {self.mikado_stats_yaml}")
        logging.info(f"Mikado stats JSON written to {self.mikado_stats_json}")
        logging.info(f"Mikado stats TOML written to {self.mikado_stats_toml}")
        logging.info(f"Mikado summary stats written to {self.mikado_basic_stats}")

    def update_oddities(self):
        """Update oddities with user-specified thresholds"""
        METRIC_ODDITIES_UPDATED = []
        for oddity in METRIC_ODDITIES:
            if "five_utr_length" in oddity:
                updated_oddity = oddity.replace(
                    "10000", str(self.args.five_prime_utr_length)
                )
                METRIC_ODDITIES_UPDATED.append(updated_oddity)
            elif "three_utr_length" in oddity:
                updated_oddity = oddity.replace(
                    "10000", str(self.args.three_prime_utr_length)
                )
                METRIC_ODDITIES_UPDATED.append(updated_oddity)
            elif "five_utr_num" in oddity:
                updated_oddity = oddity.replace("5", str(self.args.five_utr_num))
                METRIC_ODDITIES_UPDATED.append(updated_oddity)
            elif "three_utr_num" in oddity:
                updated_oddity = oddity.replace("4", str(self.args.three_utr_num))
                METRIC_ODDITIES_UPDATED.append(updated_oddity)
            elif "min_exon_length" in oddity:
                updated_oddity = oddity.replace("5", str(self.args.min_exon_length))
                METRIC_ODDITIES_UPDATED.append(updated_oddity)
            elif "max_exon_length" in oddity:
                updated_oddity = oddity.replace("10000", str(self.args.max_exon_length))
                METRIC_ODDITIES_UPDATED.append(updated_oddity)
            elif "min_intron_length" in oddity:
                updated_oddity = oddity.replace("5", str(self.args.min_intron_length))
                METRIC_ODDITIES_UPDATED.append(updated_oddity)
            elif "max_intron_length" in oddity:
                updated_oddity = oddity.replace(
                    "500000", str(self.args.max_intron_length)
                )
                METRIC_ODDITIES_UPDATED.append(updated_oddity)
            elif "selected_cds_fraction" in oddity:
                updated_oddity = oddity.replace(
                    "0.3", str(self.args.selected_cds_fraction)
                )
                METRIC_ODDITIES_UPDATED.append(updated_oddity)
            else:
                METRIC_ODDITIES_UPDATED.append(oddity)

        logging.debug(
            f"Updated METRIC_ODDITIES with user thresholds: {METRIC_ODDITIES_UPDATED}"
        )

        # print the updated oddities for user
        logging.info("Using the following oddity conditions:")
        for oddity in METRIC_ODDITIES_UPDATED:
            logging.info(f" {oddity}")

        return METRIC_ODDITIES_UPDATED

    def process_cmd(self, cmd):
        # capture exit code and log error if non-zero
        logging.info(f"Running command: {cmd}")
        exit_code = os.system(cmd)
        if exit_code != 0:
            logging.error(f"Command failed with exit code {exit_code}")
            sys.exit(exit_code)

    def combine_statistics_outputs(self):
        """
        Combine AGAT and Mikado statistics into single YAML, JSON and TOML files, but under different names AGAT and Mikado
        """
        combined_yaml = os.path.join(
            self.args.output_prefix + ".AnnoOddities.combined_statistics.yaml"
        )
        combined_json = os.path.join(
            self.args.output_prefix + ".AnnoOddities.combined_statistics.json"
        )
        combined_toml = os.path.join(
            self.args.output_prefix + ".AnnoOddities.combined_statistics.toml"
        )

        # load AGAT yaml
        with open(self.agat_statistics_yaml) as agat_yaml_file:
            agat_data = yaml.safe_load(agat_yaml_file)

        # load Mikado yaml
        with open(self.mikado_stats_yaml) as mikado_yaml_file:
            mikado_data = yaml.safe_load(mikado_yaml_file)

        # combine
        combined_data = {"AGAT": agat_data, "Mikado": mikado_data}

        # write combined yaml
        with open(combined_yaml, "w") as combined_yaml_file:
            yaml.dump(combined_data, combined_yaml_file)

        # write combined json
        with open(combined_json, "w") as combined_json_file:
            json.dump(combined_data, combined_json_file, indent=4)

        # write combined toml
        with open(combined_toml, "w") as combined_toml_file:
            toml.dump(combined_data, combined_toml_file)

        logging.info(f"Combined statistics YAML written to '{combined_yaml}'")
        logging.info(f"Combined statistics JSON written to '{combined_json}'")
        logging.info(f"Combined statistics TOML written to '{combined_toml}'")

    def _get_gff_attribute(self, attrib, field):
        """Get Note field from GFF3 attributes"""
        target_prefix = f"{field}="
        extracted_value = [
            field for field in attrib.split(";") if field.startswith(target_prefix)
        ]
        if extracted_value:
            return extracted_value[0].replace(target_prefix, "")
        return ""

    def extract_oddity_transcripts(self):
        """
        Extract transcripts with oddities from GFF3 file.
        Each oddity will have its own GFF3 file, which means there will be duplicate transcripts across multiple oddity GFF3 files if they match multiple oddities.
        Also create a combined GFF3 file with all transcripts and their oddities noted in the Note field.
        """

        # make a oddity_files folder to store all oddity gff files
        oddity_dir = os.path.join(
            os.path.dirname(self.args.output_prefix), "oddity_files"
        )
        os.makedirs(oddity_dir, exist_ok=True)

        # remove existing oddity gff files if any
        for oddity_file in os.listdir(oddity_dir):
            if oddity_file.startswith(f"{self.args.output_prefix}.AnnoOddities."):
                os.remove(os.path.join(oddity_dir, oddity_file))

        output_oddity_gff = f"{self.args.output_prefix}.AnnoOddities.gff"

        with (
            open(self.sanitised_gff) as gff_in,
            open(output_oddity_gff, "w") as out_gff_fh,
        ):
            for line in gff_in:
                if line.startswith("#"):
                    out_gff_fh.write(line)
                x = line.strip().split("\t")
                if len(x) < 9:
                    continue
                attrib = x[8]

                tid = self._get_gff_attribute(attrib, "ID")
                pid = self._get_gff_attribute(attrib, "Parent")
                # not every line will have both ID and Parent, so only proceed if at least one exist when writing oddity specific gff
                if not tid and not pid:
                    out_gff_fh.write(line)
                    continue
                skip_main_id = False
                # here I am for gene level oddities
                if tid and tid in self.oddity_gid_info:
                    # for gene level oddities, write all transcripts of the gene
                    gid = tid  # here tid is actually gene id
                    for oddity in set(self.oddity_gid_info[gid]):
                        self.write_oddity_to_gff(oddity_dir, line, oddity)
                # here I am for transcript level oddities
                elif tid and tid in self.oddity_tid_info:
                    note_value = self._get_gff_attribute(attrib, "Note")
                    for oddity in (
                        self.oddity_tid_info[tid].get("matched_oddities", "").split(";")
                    ):
                        self.write_oddity_to_gff(oddity_dir, line, oddity)
                        # update the Note field
                        if note_value:
                            note_value += f"|{oddity}"
                        else:
                            note_value = f"{oddity}"
                    # now update the line with new Note field
                    if "Note=" in attrib:
                        # replace existing Note value
                        new_attrib = re.sub(
                            r"Note=[^;]+",
                            f"Note={note_value}",
                            attrib,
                        )
                    else:
                        # append Note field
                        new_attrib = attrib + f";Note={note_value}"
                    # reconstruct the line
                    new_line = "\t".join(x[:8] + [new_attrib]) + "\n"
                    out_gff_fh.write(new_line)

                    skip_main_id = True
                # here I am for parent level oddities that are child of the transcript
                elif pid and pid in self.oddity_tid_info:
                    for oddity in (
                        self.oddity_tid_info[pid].get("matched_oddities", "").split(";")
                    ):
                        self.write_oddity_to_gff(oddity_dir, line, oddity)

                # if no oddity matched, just write the original line
                if not skip_main_id:
                    out_gff_fh.write(line)

        logging.info(
            f"Extracted transcripts with oddities written to '{output_oddity_gff}'"
        )
        logging.info(f"Individual oddity GFF3 files written to '{oddity_dir}' folder")

    def write_oddity_to_gff(self, oddity_dir, line, oddity):
        # metric_name = AnnoOddities._clean_metric_name(oddity)
        output_each_oddity_gff = os.path.join(
            oddity_dir,
            f"{self.args.output_prefix}.AnnoOddities.{oddity}.gff",
        )
        with open(output_each_oddity_gff, "a") as output_each_oddity_gff_fh:
            output_each_oddity_gff_fh.write(line)

    def write_oddities_table(self, oddity_counts):
        """Write oddities table to file"""
        with open(self.oddity_summary, "w") as f:
            f.write(f"AnnoOddities\t{self.args.output_prefix}\n")
            for oddity, count in oddity_counts.items():
                metric_name = oddity.replace("{", "").replace("}", "")
                f.write(f"{metric_name}\t{count}\n")

        logging.info(f"Oddities summary table written to '{self.oddity_summary}'")

    def run(self):
        logging.info("Finding annotation oddities")
        logging.info(f"Command: {executed_command}")

        # run agat to standardise gff3
        self.run_agat_standardise_gff()

        # run agat_sp_statistics
        self.run_agat_sp_statistics()

        # run gffread to generate table
        self.run_gffread_table()

        # run mikado tab-stats
        self.run_mikado_tab_stats()

        logging.info("Update oddities as per user thresholds")
        METRIC_ODDITIES_UPDATED = self.update_oddities()

        # parse mikado tab-stats and gffread table
        logging.info("Parsing Mikado tab-stats file")
        mikado_data = self.parse_mikado_tab_stats(self.mikado_tab_stats)

        # load genome fasta
        logging.info("Loading genome FASTA file, this may take a while...")
        self.genome_db = Fasta(self.args.genome_fasta)

        logging.info("Parsing gffread table file")
        gffread_data = self.parse_gffread_table(self.genome_db, self.gffread_table)

        # combine yaml, json and toml outputs from agat and mikado into combined yaml, json and toml
        logging.info("Combining AGAT and Mikado statistics")
        self.combine_statistics_outputs()

        logging.info("Finding annotation oddities")
        # find oddities
        oddity_counts, self.oddity_tid_info, self.oddity_gid_info = (
            self.calculate_oddities(
                mikado_data,
                gffread_data,
                METRIC_ODDITIES_UPDATED,
                self.all_stats_tsv,
            )
        )

        # use oddities info to extract from gff3 if needed
        logging.info("Separating transcripts with oddities, this may take a while...")
        self.extract_oddity_transcripts()

        logging.info(f"All oddities TSV table written to '{self.all_stats_tsv}'")
        # write oddities table
        self.write_oddities_table(oddity_counts)


def main():
    parser = argparse.ArgumentParser(
        description="Find annotation oddities from GFF3 and genome FASTA files",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    # genome fasta
    parser.add_argument(
        "--genome_fasta",
        type=str,
        required=True,
        help="Provide Genome FASTA file",
    )
    # gff3 file
    parser.add_argument(
        "--gff3_file",
        type=str,
        required=True,
        help="Provide GFF3 file with transcript annotations",
    )
    # add indpendent oddities
    # group oddities = parser.add_argument_group("Oddity Thresholds")
    group_oddities = parser.add_argument_group("Oddity Thresholds")
    # five_prime_utr_length
    group_oddities.add_argument(
        "--five_prime_utr_length",
        type=int,
        default=10000,
        help="Threshold for 5' UTR length oddity [default:%(default)s]",
    )
    # three_prime_utr_length
    group_oddities.add_argument(
        "--three_prime_utr_length",
        type=int,
        default=10000,
        help="Threshold for 3' UTR length oddity [default:%(default)s]",
    )
    # five_utr_num
    group_oddities.add_argument(
        "--five_utr_num",
        type=int,
        default=5,
        help="Threshold for 5' UTR exon number oddity [default:%(default)s]",
    )
    # three_utr_num
    group_oddities.add_argument(
        "--three_utr_num",
        type=int,
        default=4,
        help="Threshold for 3' UTR exon number oddity [default:%(default)s]",
    )
    # min_intron_length
    group_oddities.add_argument(
        "--min_intron_length",
        type=int,
        default=5,
        help="Threshold for minimum intron length oddity [default:%(default)s]",
    )
    # max_intron_length
    group_oddities.add_argument(
        "--max_intron_length",
        type=int,
        default=120000,
        help="Threshold for maximum intron length oddity.\nBelow are some guidelines:\n- For fungi species, consider setting this to 1000 bp.\n- For plant species, consider setting this to 10000 bp.\n- For invertebrates species, consider setting this to 60000 bp.\n- For vertebrates species, consider setting this to 120000 bp.\n[default:%(default)s]",
    )
    # min_exon_length
    group_oddities.add_argument(
        "--min_exon_length",
        type=int,
        default=5,
        help="Threshold for minimum exon length oddity [default:%(default)s]",
    )
    # max_exon_length
    group_oddities.add_argument(
        "--max_exon_length",
        type=int,
        default=10000,
        help="Threshold for maximum exon length oddity [default:%(default)s]",
    )
    # selected_cds_fraction
    group_oddities.add_argument(
        "--selected_cds_fraction",
        type=float,
        default=0.3,
        help="Threshold for selected CDS fraction oddity.\nThis is the proportion of coding sequence to that of the transcript.\nValues range from 0.0 to 1.0 [default:%(default)s]",
    )
    # canonical intron motifs
    group_oddities.add_argument(
        "--canonical_intron_motifs",
        type=str,
        # default=",".join(CANONICAL_INTRON_MOTIFS),
        help=f"Comma-separated list of canonical intron motifs. [default:'{','.join(CANONICAL_INTRON_MOTIFS)}']",
    )
    parser.add_argument(
        "--output_prefix",
        type=str,
        default="output",
        help="Provide sample prefix for the output table [default:%(default)s]",
    )
    # force rerun
    parser.add_argument(
        "--force",
        action="store_true",
        help="Force rerun even if output files exist",
    )
    # verbosity level for logging
    parser.add_argument(
        "--verbosity",
        type=str,
        choices=["debug", "info", "warning", "error", "critical"],
        default="info",
        help="Set logging verbosity level [default:%(default)s]",
    )
    args = parser.parse_args()

    # TODO: will use logging.getLevelNamesMapping() later if needed (Added in version 3.11)
    # log_level = logging.getLevelNamesMapping().get(args.verbosity.upper(), logging.INFO)
    log_level = logging.getLevelName(args.verbosity.upper())
    logging.basicConfig(
        format="%(asctime)s | %(levelname)-8s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        level=log_level,
    )

    # check if gffread, gt, mikado are available
    logging.info("Checking for required tools in PATH")
    for tool in [
        "agat",
        "gffread",
        "mikado",
        "agat_convert_sp_gxf2gxf.pl",
        "agat_sp_statistics.pl",
    ]:
        if not shutil.which(tool):
            logging.error(f"Required tool '{tool}' not found in PATH")
            sys.exit(1)
        else:
            logging.info(
                f"Found required tool in PATH: '{tool}' - {shutil.which(tool)}"
            )

    AnnoOddities(args).run()


if __name__ == "__main__":
    main()
