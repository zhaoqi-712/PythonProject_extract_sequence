"""
Extract the allele of each gene in multiple genomes

"""
import csv
import pysam
import pandas as pd
from collections import defaultdict
import random
import os
import glob
from Bio import SeqIO

################################
# 1. process depth result data
def process_data(input_file):
    """
    Process the depth data from a txt file (input_file).
    :param input_file: path to the analysis result (depth of alignment at each position).
    :return: rows: all rows of the txt file. Example:
    [{'seq_ID': 'NC_007194.1',
     'type': 'mRNA',
     'start': 216,
     'end': 836,
     'strand': '+',
     'depth': '50.0000000',
     'Parent': 'gene-AFUA_1G00100',
     'Name': 'XM_001481640.1',
     'locus_tag': 'AFUA_1G00100',
     'product': 'putative MFS monocarboxylate transporter',
     'transcript_id': 'XM_001481640.1'}, ...]
     """

    # Read and process data
    rows = []

    with open(input_file, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")

            # Basic columns (first 8 columns + attributes + depth)
            base_cols = [
                "seq_ID", "source", "type", "start", "end",
                "score", "strand", "phase", "attributes", "depth"
            ]
            row = dict(zip(base_cols, parts))
            row["start"] = int(row["start"])
            row["end"] = int(row["end"])

            # Parse the 'attributes' column
            attributes = row["attributes"].split(";")
            for attr in attributes:
                attr = attr.strip()
                if not attr:
                    continue
                if "=" in attr:
                    key, value = attr.split("=", 1)
                else:
                    continue
                row[key] = value

            # select column
            keys_to_keep = ['seq_ID', 'type', 'start', 'end', 'strand', 'depth', 'Parent', 'Name',
                            'locus_tag', 'product', 'transcript_id']
            row_new = {}
            for key in keys_to_keep:
                row_new[key] = row[key]

            rows.append(row_new)
    return rows


def select_depth(rows, lower_limit, upper_limit):
    """
    Only keeps the rows with depth between lower_limit and upper_limit.
    :param rows:
    :param lower_limit:
    :param upper_limit:
    :return: filtered_rows, same as rows, a list including
    dictionaries for each mRNA position
    """

    filtered_rows = []
    for row in rows:
        depth = row["depth"]
        depth = float(depth)
        if lower_limit < depth < upper_limit:
            # print(row)
            filtered_rows.append(row)

    return filtered_rows


def dict_rows_transfer(rows):
    """
    Transfer a list into dictionary
    :param rows: a list including dictionaries for each mRNA position
    :return: a dictionary including dictionaries for each mRNA position
    """
    rows_dict = {}
    for row in rows:
        name = row["Name"]
        rows_dict[name] = row
    return rows_dict


def extract_candidate_position_list(filtered_rows):
    """
    Transfer the list to a dictionary, only keep the necessary information
    :param filtered_rows: A list including each rows as a dictionary.
    :return: candidate_data_seq, a dictionary, with key:value =
    sequence_id : [candiate_info1, candiate_info2, ...]

    {'NC_007194.1':
    [{'seq_ID': 'NC_007194.1',
    'start': 29598,
    'end': 29981,
    'depth': '12.0000000',
    'id': 'XM_744675.1',
    'locus_tag': 'AFUA_1G00160'}, ...]
    xxxx: ......
    }

    """
    candidate_data_seq = {}
    for row in filtered_rows:
        if not row["seq_ID"] in candidate_data_seq.keys():
            candidate_data_seq[row["seq_ID"]] = []

        row_data = {
            "seq_ID": row["seq_ID"],
            "start": row["start"],
            "end": row["end"],
            "depth": row["depth"],
            "id": row["Name"],
            "locus_tag": row["locus_tag"]
        }
        candidate_data_seq[row["seq_ID"]].append(row_data)

    return candidate_data_seq


def read_gff(gff_path, keep_type=None):
    """
    Read a GFF3 file and organize it by seq_ID.

    :param gff_path: path to the GFF3 file
    :param keep_type: optional, only keep rows with this feature type (e.g., "exon")
    :return: a dictionary, {seq_ID: [row_dict, ...]} sorted by start position
    """
    gff_dict = defaultdict(list)

    with open(gff_path, "r", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f, delimiter="\t", fieldnames=[
            "seq_ID", "source", "type", "start", "end", "score", "strand", "phase", "attributes"
        ])

        for row in reader:
            if row["seq_ID"].startswith("#"):
                continue  # skip header/comment lines
            if keep_type and row["type"] != keep_type:
                continue

            # parse attributes (GFF3 format: key=value;key=value;...)
            attr_dict = {}
            for attr in row["attributes"].split(";"):
                if attr.strip() == "":
                    continue
                if "=" in attr:
                    key, value = attr.strip().split("=", 1)
                    attr_dict[key] = value
            row_data = {
                "seq_ID": row["seq_ID"],
                "type": row["type"],
                "start": int(row["start"]),
                "end": int(row["end"]),
                "id": attr_dict.get("ID"),
                "locus_tag": attr_dict.get("locus_tag"),
                "transcript_id": attr_dict.get("transcript_id")
            }

            gff_dict[row["seq_ID"]].append(row_data)

    # sort each list by start position
    for seq in gff_dict:
        gff_dict[seq].sort(key=lambda x: x["start"])

    annotation_sorted = annotation_rank(dict(gff_dict))

    return annotation_sorted


def annotation_rank(annotation_sorted):
    for chromosome, annotation_info in annotation_sorted.items():
        i = 1
        for annotation in annotation_info:
            annotation["rank"] = i
            i += 1
    return annotation_sorted


def read_gff_dict(annotation_sorted):
    """
    Convert annotation storage format.
    :param annotation_sorted: a dictionary, {seq_ID: [row_dict, ...]} sorted by start position
    :return: dict, {seq_ID: {transcript_id: annotation, ...}}
    """
    annotation_dict = {}
    for seq, annotation_list in annotation_sorted.items():
        annotation_dict[seq] = {}
        for annotation in annotation_list:
            annotation_id = annotation["transcript_id"]
            annotation_dict[seq][annotation_id] = annotation

    return annotation_dict


def merge_candidate_position(candidate_data_seq, annotation_dict):
    """
    Merge adjacent or only one-position-separated candidate regions into a
    larger candidate region.
    :param candidate_data_seq:
    {'NC_007194.1':
    [{'seq_ID': 'NC_007194.1',
    'start': 29598,
    'end': 29981,
    'depth': '12.0000000',
    'id': 'XM_744675.1',
    'locus_tag': 'AFUA_1G00160'}, ...]
    xxxx: ......
    }
    :param annotation_dict: dict, genome annotation
    {seq_ID: {transcript_id: annotation, ...}}
    :return: candidate_merge, dict
    {
    NC_007197.1:
    [{'seq_ID': 'NC_007197.1',
    'region_name': 'XM_741280.1-XM_741278.1',
    'start_gene': 'XM_741280.1',
    'end_gene': 'XM_741278.1',
    'start': 223151,
    'end': 228584,
    'rank': 72,
    'gene_included': ['XM_741280.1', 'XM_741279.1', 'XM_741278.1'],
    'gene_number': 3
    }, ...],
    seq_ID_xxx: xxx,
    ......
    }
    """
    candidate_merge = {}

    for seq, candidates in candidate_data_seq.items():
        annotation_seq = annotation_dict[seq]
        candidate_merge[seq] = []
        mixed_data = {}

        for i in range(0, len(candidates)):
            candidate_i = candidates[i]
            gene_id_i = candidate_i["id"]
            annotation_candidate_i = annotation_seq[gene_id_i]

            if mixed_data == {}:  # if nothing in mixed data, select the current annotation

                # if there is not a mixed data, import this annotation as the mixed data
                mixed_data = {
                    "seq_ID": seq,
                    "region_name": gene_id_i,
                    "start_gene": gene_id_i,
                    "end_gene": gene_id_i,
                    "start": annotation_candidate_i["start"],
                    "end": annotation_candidate_i["end"],
                    "rank": annotation_candidate_i["rank"],
                    "gene_included": [gene_id_i],
                    "gene_number": 1
                }

            # Then there is already a mixed data from i-n to i
            start_i = mixed_data["start"]
            rank_i = mixed_data["rank"]

            # if there is only one candidate in the genomic sequence
            if len(candidates) == 1:
                candidate_merge[seq].append(mixed_data)
                break

            if i == len(candidates) - 1:
                candidate_merge[seq].append(mixed_data)
                break

            elif i <= len(candidates) - 2:
                # candidate i + 1
                candidate_i_1 = candidates[i + 1]
                id_i_1 = candidate_i_1["id"]
                annotation_candidate_i_1 = annotation_seq[id_i_1]
                rank_i_1 = annotation_candidate_i_1["rank"]

                # compare rank of i+1 to i
                if abs(rank_i_1 - rank_i) <= 2:  # mix
                    mixed_data_new = {
                        "seq_ID": mixed_data["seq_ID"],
                        "region_name": f"{mixed_data["start_gene"]}-{id_i_1}",
                        "start_gene": mixed_data["start_gene"],
                        "end_gene": id_i_1,
                        "start": start_i,
                        "end": annotation_candidate_i_1["end"],
                        "rank": rank_i_1,
                        "gene_included": mixed_data["gene_included"]
                    }
                    mixed_data_new["gene_included"].append(id_i_1)
                    mixed_data_new["gene_number"] = len(mixed_data_new["gene_included"])

                    mixed_data = mixed_data_new

                else:
                    candidate_merge[seq].append(mixed_data)
                    # print(mixed_data["region_name"])
                    mixed_data = {}

    return candidate_merge


####################################
# Generate covered and uncovered genomes at one position.
# Check sequences covered or uncovered on specific genomic region of BAM file
def extract_align_seq(bam_path, seq, start, end, full_cover=True):
    """
    Extract the sequence name that are aligned to one genomic position (chr, start, end) of the BAM file.
    The sequences are from different genome assemblies.
    :arg bam_path: Path to the BAM file.
    :arg seq: Sequence name (in reference genome) of the genomic position.
    :arg start: Start position.
    :arg end: End position.

    :return: align_info
    a dictionary that include list of the sequence names, each in a dictionary.
    Information including NCBI accession + assembly name + sequence name + position in that assembly.
    """

    align_info = {}
    # the position in the reference genome
    align_info["pos_info"] = {"seq": seq, "start": start, "end": end}
    # the list to store the assembly names
    align_info["seq_info"] = []

    bam = pysam.AlignmentFile(bam_path, "rb")

    for read in bam.fetch(seq, start, end):
        # Check if the read overlaps the region of interest
        # read.query_alignment_start / read.query_alignment_end are the aligned
        # positions on the read sequence
        # Coordinates are 0-based, left-inclusive and right-exclusive
        # if read.is_unmapped or read.is_secondary or read.is_supplementary:
        # continue

        # Only reads that cover the interval 100% are retained
        if full_cover:
            if read.reference_start > start or read.reference_end < end:
                continue

        # Create a dictionary to store read alignment information
        dic_align = {}
        name_all = read.query_name
        name_parse = name_all.split("_")
        ncbi_accession = name_parse[0] + "_" + name_parse[1]
        seq_id = name_parse[-1]

        if read.is_reverse == False:
            orientation = "forward"
            # the reads corridinates that align to the start-end of reference.
            read_corridinate_start = find_read_corridinate(read, start)
            read_corridinate_end = find_read_corridinate(read, end)

        else:
            orientation = "reverse"
            # the reads corridinates that align to the start-end of reference.
            read_corridinate_start = find_read_corridinate(read, end)
            read_corridinate_end = find_read_corridinate(read, start)

        # dic_align["name"] = name_all      # Complete read name
        dic_align["Genome_accession"] = ncbi_accession  # Read NCBI accession number
        dic_align["seq_ID"] = seq_id  # ID of the sequence, such as chromosome ID
        # if reverse then
        dic_align["start"] = read_corridinate_start  # Alignment start on read
        dic_align["end"] = read_corridinate_end  # Alignment end on read
        dic_align["orientation"] = orientation  # True if read aligns to reverse strand

        align_info["seq_info"].append(dic_align)
    # print(len(align_info["seq_info"]))
    bam.close()
    return align_info


def find_read_corridinate(read, position):
    """
    finding the corridinate of read on the position of the reference genome.

    :param read:
    :param position:
    :return:
    """
    read_self_start = read.query_alignment_start
    read_reference_start = read.reference_start  # aligned reads on reference genome

    read_cigar_list = read.cigartuples

    insertion = 0
    deletion = 0
    length_reference = 0
    hard_mask = 0
    length_ref_0 = position - read_reference_start + 1

    read_length_hard = 0
    for info in read_cigar_list:
        if info[0] == 5:  # hard mask
            read_length_hard += info[1]
    read_length = read.query_length + read_length_hard

    read_length_2 = 0
    for info in read_cigar_list:
        info = list(info)
        type_seq = str(info[0])
        length_seq = info[1]

        if type_seq == "2":  # deletion
            continue
        else:
            read_length_2 += length_seq

    for info in read_cigar_list:
        info = list(info)
        type_seq = str(info[0])
        length_seq = info[1]
        # print(type_seq,length_seq)
        if type_seq == "5":  # hard mask
            hard_mask = length_seq
        elif type_seq == "4":  # soft mask
            continue
        elif type_seq == "1":  # insertion
            insertion += length_seq
        elif type_seq == "2":  # deletion
            deletion += length_seq
            length_reference += length_seq
        elif type_seq in ["0", "6", "7", "8", "9"]:  # other situations, there are also type > 5, not considered
            length_reference += length_seq

        # what if 3, skipped region?

        if length_reference > length_ref_0:
            break

    read_corridinate = hard_mask + read_self_start + length_ref_0 + insertion - deletion

    if read.is_reverse:
        read_corridinate = read_length - read_corridinate + 1

    return (read_corridinate)


def extract_align_seq_from_dict(bam_path, position_info):
    """
    using function extract_align_seq() with a dictionary that extracted from annotation.
    :param bam_path:
    :param position_info:
    :return: align_info_B
    """
    seq = position_info["seq_ID"]
    start = position_info["start"]
    end = position_info["end"]
    align_info_B = extract_align_seq(bam_path, seq, start, end)
    return align_info_B


def find_aligned_assembly(align_info):
    aligned_genome = []
    for seq_info in align_info["seq_info"]:
        ID = seq_info["Genome_accession"]
        aligned_genome.append(ID)

    return aligned_genome


def find_not_aligned_assembly(aligned_genome, align_info, assembly_path):
    """
    Find the genome assemblies that are not aligned to one genomic position.
    :param align_info: The output of extract_align_seq, with the list (align_info["seq_info"])
    consist of dictionaries of the mapped genome assemblies information at certain position.
    :param assembly_path: the path to a txt file, each row is a NCBI accession number of genome assembly.
    A list consisting of all genome assemblies that were used in the alignment.
    :return: not_align_list, A list including all uncovered genome assemblies at the genomic position.
    """

    # load the list including all genome assemblies
    with open(assembly_path, "r", encoding="utf-8") as f:
        all_genome = [line.strip() for line in f]
    # print(all_genome)
    not_align_list = list(set(all_genome) - set(aligned_genome))
    position = align_info["pos_info"]
    # print(f"Position on reference genome: sequence ID: {position['seq']}, "
    # f"start position: {position['start']}, "
    # f"end position: {position['end']}.\n"
    # f"Overall {len(all_genome)} assemblies, aligned {len(aligned_genome)} assemblies, "
    # f"not aligned {len(not_align_list)} assemblies.")
    return not_align_list


#######################################

# filter the input file it between up and down limits


def extract_candidate_position(filtered_rows):
    """
    Extract the necessary position information from the filtered data.
    :param filtered_rows: A list including each rows as a dictionary.
    :return: candidate_data, a dictionary, including each mRNA candidates and its information as key-value pairs.
    """
    candidate_data = {}
    for row in filtered_rows:
        row_data = {
            "seq_ID": row["seq_ID"],
            "start": row["start"],
            "end": row["end"],
            "depth": row["depth"],
            "id": row["Name"],
            "locus_tag": row["locus_tag"]
        }
        candidate_data[row["Name"]] = row_data
    return candidate_data


####################################
# data analysis
# Select up and downstream positions
def read_gtf(gtf_path, keep_type=None):
    """
    Read a GTF/GFF file and organize it by seq_ID.

    :param gtf_path: path to the GTF/GFF file
    :param keep_type: optional, only keep rows with this feature type (e.g., "exon")
    :return: dict, {seq_ID: [row_dict, ...]} sorted by start
    """
    gtf_dict = defaultdict(list)

    with open(gtf_path, "r", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f, delimiter="\t", fieldnames=[
            "seq_ID", "source", "type", "start", "end", "score", "strand", "phase", "attributes"
        ])

        for row in reader:
            if row["seq_ID"].startswith("#"):
                continue  # skip header/comment lines
            if keep_type and row["type"] != keep_type:
                continue

            # parse attributes
            attr_dict = {}
            for attr in row["attributes"].split(";"):
                if attr.strip() == "":
                    continue
                key_value = attr.strip().split(" ", 1)
                if len(key_value) == 2:
                    key, value = key_value
                    attr_dict[key] = value.strip('"')

            row_data = {
                "seq_ID": row["seq_ID"],
                "type": row["type"],
                "start": int(row["start"]),
                "end": int(row["end"]),
                "gene_id": attr_dict.get("gene_id"),
                "transcript_id": attr_dict.get("transcript_id")
            }

            gtf_dict[row["seq_ID"]].append(row_data)

    # sort each row with the order of start position
    for seq in gtf_dict:
        gtf_dict[seq].sort(key=lambda x: x["start"])

    return dict(gtf_dict)


def find_position_seq(sequence_id, start, end, annotation_sorted, up_num, down_num):
    """
    Retrieve the positional information of the m upstream and n downstream mRNAs for a specified mRNA.
    :param transcript_id: NCBI id, such as "XM_742446.1"
    :param sequence_id: sequence of the transcript, NCBI id, such as "NC_007194.1"
    :param annotation_sorted: the sorted annotation dictionary.
    :param n: number of positions to identify at the up and down-stream positions.
    :return: up_down_locations: a dictionary with the location of input mRNA, its up and down-stream mRNA positions.
    {
        'position': {'seq_ID': 'NC_007195.1', 'start': 4650883, 'end': 4653737},
        'type': 'mRNA',
        'upstream_position': [
            {'seq_ID': 'NC_007195.1', 'type': 'mRNA', 'start': 4643751, 'end': 4644257, 'id': 'rna-XM_750980.1',
             'parent': 'gene-AFUA_2G17380', 'gene_id': None, 'transcript_id': 'XM_750980.1', 'rank': 1555}, ...],
        'downstream_position': [
            {'seq_ID': 'NC_007195.1', 'type': 'mRNA', 'start': 4653899, 'end': 4654978, 'id': 'rna-XM_750985.1',
             'parent': 'gene-AFUA_2G17430', 'gene_id': None, 'transcript_id': 'XM_750985.1', 'rank': 1560}, ...]
    }
    """
    annotation_sequences = annotation_sorted[sequence_id]
    annotation_length = len(annotation_sequences)
    up_position = 0
    down_position = 0

    for i in range(1, annotation_length):
        if annotation_sequences[i - 1]["start"] <= start <= annotation_sequences[i]["start"]:
            up_position = i
        if annotation_sequences[i - 1]["end"] <= end <= annotation_sequences[i]["end"]:
            down_position = i - 1
        if i == annotation_length - 1:
            break

    if up_position > 0 and down_position > 0:
        # overall number of positions in up and down stream positions
        number_upstream = up_position
        number_downstream = annotation_length - down_position - 1

        # Extract
        upstream_position = []
        downstream_position = []
        # extract upstream positions
        if number_upstream >= up_num:
            for j in range(up_position - up_num, up_position):  # m upstream positions
                upstream_position.append(annotation_sequences[j])
        elif up_num > number_upstream >= 1:
            for j in range(0, up_position):  # m upstream positions
                upstream_position.append(annotation_sequences[j])

        # extract downstream positions
        if number_downstream >= down_num:
            for j in range(down_position + 1, down_position + down_num + 1):  # n downstream positions
                downstream_position.append(annotation_sequences[j])
        elif down_num > number_downstream >= 1:
            for j in range(down_position + 1, annotation_length):  # n downstream positions
                downstream_position.append(annotation_sequences[j])
        position_info = {
            'seq_ID': sequence_id,
            'start': start,
            'end': end
        }

        up_down_locations_seq = {
            "position": position_info,
            'type': 'mRNA',
            "upstream_position": upstream_position,
            "downstream_position": downstream_position
        }
        return up_down_locations_seq

    else:
        print(
            f"Warning! sequence {sequence_id},start {start}, end {end} is not found in {sequence_id} of reference genome.")


def find_position_depth(transcript_id, rows_dict):
    # print(transcript_id)
    pos_info = rows_dict[transcript_id]
    pos_depth = float(pos_info["depth"])
    return pos_depth


def filter_up_down_depth(up_down_locations, rows, up_num, down_num, cutoff=10):
    """

    :param up_down_locations:
    :param rows_dict:
    :param up_num:
    :param down_num:
    :param cutoff:
    :return:
    """
    rows_dict = dict_rows_transfer(rows)
    upstream_position = up_down_locations["upstream_position"]
    downstream_position = up_down_locations["downstream_position"]
    sum_depth_up = 0
    sum_depth_down = 0
    for pos1 in upstream_position:
        # print(pos1)
        transcript_id1 = pos1["transcript_id"]
        sum_depth_up += find_position_depth(transcript_id1, rows_dict)
    for pos2 in downstream_position:
        transcript_id2 = pos2["transcript_id"]
        sum_depth_down += find_position_depth(transcript_id2, rows_dict)
    ave_depth_up = sum_depth_up / up_num
    ave_depth_down = sum_depth_down / down_num
    if ave_depth_up >= cutoff or ave_depth_down >= cutoff:
        return True
    else:
        return False


def find_candidate_align(up_down_locations, bam_path, assembly_path):
    """

    :param up_down_locations: a dictionary with the location of input mRNA, its up and down-stream mRNA positions.
    :return: up_down_alignment:
    position_info_align, alignment and un-alignment lists at the position.
    up_gene_align, alignment at upstream positions.
    down_gene_align, alignment at downstream positions
    """
    up_down_alignment = {}
    position_info = up_down_locations["position"]  # one dictionary, analysis both align and unaligned genomes
    upstream_position = up_down_locations["upstream_position"]  # list including dictionaries
    downstream_position = up_down_locations["downstream_position"]  # list including dictionaries

    # at the candidate position
    position_align_info = extract_align_seq_from_dict(bam_path, position_info)
    position_involved_assembly = find_aligned_assembly(position_align_info)
    position_uninvolved_assembly = find_not_aligned_assembly(position_involved_assembly, position_align_info,
                                                             assembly_path)

    upstream_position_involved_assembly = []
    for pos_up in upstream_position:
        position_align_info_up = extract_align_seq_from_dict(bam_path, pos_up)
        upstream_position_involved_assembly.append(position_align_info_up)

    downstream_position_involved_assembly = []
    for pos_down in downstream_position:
        position_align_info_down = extract_align_seq_from_dict(bam_path, pos_down)
        downstream_position_involved_assembly.append(position_align_info_down)

    up_down_alignment = {
        "position_assembly":  # at the candidate position
            {"involved_assembly": position_involved_assembly,
             "uninvolved_assembly": position_uninvolved_assembly
             },
        # aligned assemblies at the up and down positions
        "upstream_position_assembly": upstream_position_involved_assembly,
        "downstream_position_assembly": downstream_position_involved_assembly
    }
    return up_down_alignment


# Function: randomly sample n unaligned sequences from the up_down_alignment of a given gene,
# and then analyze the state at each position for each of them
def random_select_assembly(info_list, assembly_num):
    """
    Randomly sample n unaligned sequences at the interested mRNA position.
    This function is currently not used!
    :param up_down_alignment: 改成对list适用的
    :param n:
    :return: selected_assemblies, list of selected assemblies
    """
    assemblies = list(info_list.keys())
    number_assemblies = len(assemblies)
    if number_assemblies == 0:
        info_selected = {}

    elif 0 < number_assemblies <= assembly_num:
        info_selected = info_list

    elif number_assemblies > assembly_num:
        assemblies_selected = random.sample(assemblies, assembly_num)

        info_selected = {}

        for key, value in info_list.items():
            if key in assemblies_selected:
                info_selected[key] = value

    return info_selected


def find_candidate_status(assembly, reads):
    status = []
    for positions in reads:  # analyze each position
        seq_info_list = positions["seq_info"]
        genome_accessions = [d["Genome_accession"] for d in seq_info_list]
        if assembly in genome_accessions:
            status.append(True)
        else:
            status.append(False)
    return status


def find_assembly_status(selected_assemblies, upstream_reads, downstream_reads):
    assemblies_status = {}
    for assembly in selected_assemblies:  # test each selected assembly in up and downstream positions
        assembly_status = {
            "upstream_status": find_candidate_status(assembly, upstream_reads),
            "downstream_status": find_candidate_status(assembly, downstream_reads)
        }
        assemblies_status[assembly] = assembly_status

    return assemblies_status


# Finding whether the selected assemblies are involved in the positions.
def find_candidate_involvement(up_down_alignment):
    """

    :param up_down_alignment:
    :param
    :return:
    """
    # get the random assembly list to test status
    assemblies_ref_allele = up_down_alignment["position_assembly"]["involved_assembly"]
    assemblies_diff_allele = up_down_alignment["position_assembly"]["uninvolved_assembly"]

    # get the position information
    upstream_reads = up_down_alignment["upstream_position_assembly"]
    downstream_reads = up_down_alignment["downstream_position_assembly"]

    all_status = {
        "ref_allele": find_assembly_status(assemblies_ref_allele, upstream_reads, downstream_reads),
        "diff_allele": find_assembly_status(assemblies_diff_allele, upstream_reads, downstream_reads)
    }
    return all_status


def find_up_down_loci_one_status(one_status, up_down_locations, up_down_alignment):  #
    """

    :param all_status:
    :param up_down_locations:
    :param up_down_alignment:
    :return: up_down_loci
    """
    up_down_loci = {}
    # the locus information related to this locus in th reference genome
    upstream_positions = up_down_locations["upstream_position"]
    downstream_positions = up_down_locations["downstream_position"]
    # the aligned genomes onto the locus
    upstream_align = up_down_alignment["upstream_position_assembly"]
    # print(upstream_align)
    downstream_align = up_down_alignment["downstream_position_assembly"]

    for genome_id, status in one_status.items():

        # print(genome_id,status)
        loci_selected = {}
        upstream_status = status["upstream_status"]
        downstream_status = status["downstream_status"]
        up_number = len(upstream_status)
        down_number = len(downstream_status)

        # incase situations without any findings
        loci_selected["upstream_rank"] = False
        loci_selected["upstream_ref_position"] = False
        loci_selected["upstream_read_sequence_ID"] = False
        loci_selected["upstream_read_position"] = False
        loci_selected["downstream_rank"] = False
        loci_selected["downstream_ref_position"] = False
        loci_selected["downstream_read_sequence_ID"] = False
        loci_selected["downstream_read_position"] = False

        for j in range(1, up_number + 1):  # the i th nearest loci in the upstream positions
            # find upstream true locus
            status_loci_up = upstream_status[j - 1]
            if status_loci_up:  # if status_loci == True
                loci_selected["upstream_rank"] = j
                loci_selected["upstream_ref_position"] = upstream_positions[j - 1]

                # extract the information of the position in the selected assembly
                up_align_list = upstream_align[j - 1]["seq_info"]
                aligned_info_seq = []
                aligned_info_pos = {}
                for genome in up_align_list:
                    if genome["Genome_accession"] == genome_id:
                        # print(genome)
                        sequence_ID = genome["seq_ID"]
                        aligned_info_seq.append(sequence_ID)
                        # Assume that each sequence appears only once
                        aligned_info_pos[sequence_ID] = genome

                # print(aligned_info)
                loci_selected["upstream_read_sequence_ID"] = aligned_info_seq
                # print("aligned_info_seq",aligned_info_seq)
                loci_selected["upstream_read_position"] = aligned_info_pos
                break  # exit the loop once condition is met

        for i in range(1, down_number + 1):  # the i th nearest loci in the downstream positions
            # find downstream true locus
            status_loci_down = downstream_status[-i]
            if status_loci_down:  # if status_loci == True
                loci_selected["downstream_rank"] = i
                # The position information of the loci on the reference genome
                loci_selected["downstream_ref_position"] = downstream_positions[-i]
                ##
                # extract the information of the position in the selected assembly
                down_align_list = downstream_align[-i]["seq_info"]
                aligned_info_seq = []
                aligned_info_pos = {}
                for genome in down_align_list:
                    if genome["Genome_accession"] == genome_id:
                        # print(genome)
                        sequence_ID = genome["seq_ID"]
                        aligned_info_seq.append(sequence_ID)
                        # Assume that each sequence appears only once
                        aligned_info_pos[sequence_ID] = genome

                # print(aligned_info)
                loci_selected["downstream_read_sequence_ID"] = aligned_info_seq
                loci_selected["downstream_read_position"] = aligned_info_pos

                break  # exit the loop once condition is met

        up_seq = loci_selected["upstream_read_sequence_ID"]
        down_seq = loci_selected["downstream_read_sequence_ID"]

        loci_selected["seq_chromosome"] = False
        if up_seq and down_seq:
            common_seq = set(up_seq).intersection(down_seq)

            if len(common_seq) > 0:
                loci_selected["seq_chromosome"] = common_seq

                # Assume there is only one shared seq in common_seq
                seq = list(common_seq)[0]

                # calculate interval between two reference genes in read genome
                read_up_position = loci_selected["upstream_read_position"][seq]
                read_down_position = loci_selected["downstream_read_position"][seq]
                orientation = read_up_position["orientation"]

                read_interval = abs(read_down_position["end"] - read_up_position["start"])
                loci_selected["read_interval"] = read_interval
                # print(read_up_position,'\t',read_down_position, interval)

                # calculate interval between two reference genes in reference genome
                ref_up_position = loci_selected["upstream_ref_position"]
                ref_down_position = loci_selected["downstream_ref_position"]
                ref_interval = ref_down_position["end"] - ref_up_position["start"]
                loci_selected["ref_interval"] = ref_interval

                interval_diff_pct = abs(ref_interval - read_interval) * 100 / ref_interval
                loci_selected["interval difference(%)"] = interval_diff_pct

        up_down_loci[genome_id] = loci_selected
        # print(loci_selected)

    return up_down_loci


def filter_up_down_loci(up_down_loci, interval_difference=250):
    filtered_up_down_loci = {}
    for genome_id, loci_selected in up_down_loci.items():
        skip_loop = False
        for key, value in loci_selected.items():
            if value == False:
                skip_loop = True
        if skip_loop:
            continue

        if loci_selected["interval difference(%)"] >= interval_difference:
            continue

        filtered_up_down_loci[genome_id] = loci_selected

    return filtered_up_down_loci


def find_up_down_loci(all_status, up_down_locations, up_down_alignment):
    ref_allele_status = all_status["ref_allele"]
    diff_allele_status = all_status["diff_allele"]

    ref_up_down_loci = find_up_down_loci_one_status(ref_allele_status, up_down_locations, up_down_alignment)
    ref_up_down_loci_filtered = filter_up_down_loci(ref_up_down_loci, 250)
    diff_up_down_loci = find_up_down_loci_one_status(diff_allele_status, up_down_locations, up_down_alignment)
    diff_up_down_loci_filtered = filter_up_down_loci(diff_up_down_loci, 250)
    all_up_down_loci = {
        "ref_up_down_loci": ref_up_down_loci_filtered,
        "diff_up_down_loci": diff_up_down_loci_filtered
    }

    return all_up_down_loci


def count_aligned_reads(all_up_down_loci):
    # count the number of reads completely aligned to the interval.
    ref_up_down_loci = all_up_down_loci["ref_up_down_loci"]
    diff_up_down_loci = all_up_down_loci["diff_up_down_loci"]
    ref_assembly_number = len(ref_up_down_loci.keys())
    diff_assembly_number = len(diff_up_down_loci.keys())
    aligned_reads_number = {
        "ref_assembly_number": ref_assembly_number,
        "diff_assembly_number": diff_assembly_number,
        "all_assembly_number": ref_assembly_number + diff_assembly_number
    }

    return aligned_reads_number


def analyze_all_candidate_position(candidate_data, annotation_sorted, rows, bam_path, assembly_path, up_num, down_num,
                                   lower_limit):
    """
    Analyze the position data for each of the candidate positions.
    :param candidate_data: the merged candidate data.
    :param annotation_sorted: genome annotation.
    :return: candidate_data_positions: a list, including analysis results of the position data
    [{
        "region_name": XM_750984.1,
        "candidate_data": candidate, # information of this candidate, {seq: [merged_info1, ...], ...}
        "position_info": up_down_locations, # find upstream and downstream positions, dict
        "align_info": up_down_alignment, # find reads aligned to these positions,dict
        "status_info": all_status, # select assemblies and find whether they are present at these positions
        "up_down_loci": all_up_down_loci, # the selected locations
        "aligned_reads_number": aligned_reads_number
    }, ......]
    """
    candidate_data_summary = []
    for seq_id, candidates in candidate_data.items():

        for candidate in candidates:
            sequence_id = seq_id
            start = candidate["start"]
            end = candidate["end"]
            # finding the up and downstream positions
            up_down_locations = find_position_seq(sequence_id, start, end, annotation_sorted, up_num, down_num)

            if up_down_locations == None:
                continue

            # check the mean depth of the positions

            depth_status = filter_up_down_depth(up_down_locations, rows, up_num, down_num, lower_limit)
            if depth_status == False:
                continue

            # finding the genes aligned to the positions
            up_down_alignment = find_candidate_align(up_down_locations, bam_path, assembly_path)

            # identify status of random genome assemblies at the positions
            all_status = find_candidate_involvement(up_down_alignment)

            # identify the most distant upstream and downstream loci
            all_up_down_loci = find_up_down_loci(all_status, up_down_locations, up_down_alignment)

            aligned_reads_number = count_aligned_reads(all_up_down_loci)

            # Filter the complicated genomic region where not many reads completely aligned to
            all_aligned_reads_number = aligned_reads_number["all_assembly_number"]
            if all_aligned_reads_number < 15:
                continue

            # Check whether both alleles exist
            ref_allele_info = all_up_down_loci["ref_up_down_loci"]
            diff_allele_info = all_up_down_loci["diff_up_down_loci"]

            if len(ref_allele_info) == 0 or len(diff_allele_info) == 0:
                continue

            summary = {
                "region_name": candidate["region_name"],
                "candidate_data": candidate,
                "position_info": up_down_locations,
                "align_info": up_down_alignment,
                "status_info": all_status,
                "up_down_loci": all_up_down_loci,
                "aligned_reads_number": aligned_reads_number
            }

            candidate_data_summary.append(summary)

    return candidate_data_summary


#############################################################
# extract the sequence of the allele genomic region for each genes.

def extract_allele_sequence(genome_assembly_path, candidate_gene, genome_accession, seq, start, end, orientation,
                            output_path):
    """
    Extract the sequence of specific position of a genome assembly
    :param assembly_dir:
    :param candidate_gene:
    :param genome_accession:
    :param seq:
    :param start:
    :param end:
    :param orientation:
    :param output_path:
    :return:
    """
    # extract sequence
    target_seq = None
    for record in SeqIO.parse(genome_assembly_path, "fasta"):
        # print(record.id)
        if seq in record.id:
            # print("find")
            target_seq = record.seq[start - 1:end]  # Biopython is 0-based
            break

    if target_seq is None:
        warning = f"Warning: Sequence {seq} not found in {genome_accession}"
        print(warning)
        return warning

    if orientation.lower() == "reverse":
        target_seq = target_seq.reverse_complement()

    # output path
    output_dir = os.path.join(output_path, candidate_gene)
    os.makedirs(output_dir, exist_ok=True)

    output_file = os.path.join(output_dir, f"{candidate_gene}_{genome_accession}_{seq}-{start}-{end}.fa")

    # 4. save results
    with open(output_file, "w") as f:
        f.write(f">{genome_accession}_{seq}:{start}-{end}\n")
        f.write(str(target_seq) + "\n")

    print(f"Sequence saved to {output_file}")
    return output_file


def find_genome_assembly_path(assembly_dir, genome):
    assembly_pattern = os.path.join(assembly_dir, f"{genome}*.fna")
    matches = glob.glob(assembly_pattern)
    if not matches:
        warning = f"Warning: No genome assembly found for {genome}"
        print(warning)
        # return warning
    genome_assembly_path = matches[0]
    return genome_assembly_path


def extract_region_seq(allele_info, region_name, assembly_dir, output_path, extend, assembly_num):
    """
    Extract the sequence of specific position of a genome assembly
    :param allele_info:
    :param region_name:
    :param assembly_dir:
    :param output_path:
    :param extend:
    :param assembly_num:
    :return:
    """

    allele_info_selected = random_select_assembly(allele_info, assembly_num)

    for genome, info in allele_info_selected.items():
        seq = info["seq_chromosome"]
        seq_info = list(seq)[0]  # assume there is only one seq

        upstream_read_position = info["upstream_read_position"][seq_info]
        downstream_read_position = info["downstream_read_position"][seq_info]

        start_read = min(upstream_read_position["start"], downstream_read_position["start"])
        end_read = max(upstream_read_position["end"], downstream_read_position["end"])
        orientation = upstream_read_position["orientation"]

        # if end_read - start_read > 100000:
        # continue

        # extend the start-end interval
        start_read = max(1, start_read - extend)
        end_read = end_read + extend  # require modified

        # finding genome assembly path
        genome_assembly_path = find_genome_assembly_path(assembly_dir, genome)

        extract_allele_sequence(
            genome_assembly_path,
            region_name,
            genome,
            seq_info,
            start_read,
            end_read,
            extend,
            orientation,
            output_path
        )

    return True


def find_allele_sequence_inbetween(assembly_dir, candidate_data_summary, output_path, extend, assembly_num):
    """

    :param reference_genome:
    :param assembly_dir:
    :param candidate_data_summary:
    summary = {
                "region_name": candidate["region_name"],
                "position_info": up_down_locations,
                "align_info": up_down_alignment,
                "status_info": all_status,
                "up_down_loci": all_up_down_loci,
                "aligned_reads_number": aligned_reads_number
            }
    :param output_path:
    :return:
    """
    for summary in candidate_data_summary:  # the summary information of each candidate gene
        region_name = summary["region_name"]

        ref_allele_info = summary["up_down_loci"]["ref_up_down_loci"]
        diff_allele_info = summary["up_down_loci"]["diff_up_down_loci"]

        ref_extract = extract_region_seq(ref_allele_info, region_name, "ref_allele", assembly_dir, output_path, extend,
                                         assembly_num)
        diff_extract = extract_region_seq(diff_allele_info, region_name, "diff_allele", assembly_dir, output_path,
                                          extend, assembly_num)

    return True


def extract_reference_allele(candidate_data_summary, reference_genome, annotation_sorted, output_path, extend,
                             ref_assembly):
    """

    :param candidate_data_summary:
    :param reference_genome:
    :param annotation_sorted: {seq_ID: [row_dict, ...]}
    :param output_path:
    :return:
    """
    for summary in candidate_data_summary:
        region_name = summary["region_name"]
        up_down_locations = summary["position_info"]
        print(up_down_locations)

        """"""
        seq_info_ref = up_down_locations["position"]["seq_ID"]
        start_ref = up_down_locations["upstream_position"][0]["start"]
        end_ref = up_down_locations["downstream_position"][-1]["end"]

        # calculate the start and end position of the extraction region
        # start position
        start = max(1, start_ref - extend)

        # end position
        annotation_info = annotation_sorted[seq_info_ref]
        annotation_end = annotation_info[-1]["end"]
        end = min(end_ref + extend, annotation_end)

        # extract the sequence from the reference genome
        extract_allele_sequence(
            ref_assembly,
            f"{region_name}_reference_genome",
            reference_genome,
            seq_info_ref,
            start,
            end,
            "forward",
            output_path
        )

        # --- Read and filter GFF annotations ---
        # annotation_sorted {seq_ID: [row_dict, ...]}
        extract_annotation = []
        for annotation in annotation_info:
            start_anno = annotation["start"]
            end_anno = annotation["end"]

            if start_anno >= start and end_anno <= end:
                extract_annotation.append(annotation)

        # create output path for gff3 file
        output_dir = os.path.join(output_path, region_name)
        output_file = os.path.join(output_dir,
                                   f"{region_name}_reference_genome_{reference_genome}_{seq_info_ref}-{start}-{end}.gff3")

        # write in gff3 formate file
        with open(output_file, "w", encoding="utf-8", newline="") as out:
            out.write("##gff-version 3\n")

            writer = csv.writer(out, delimiter="\t", lineterminator="\n")

            for ann in extract_annotation:
                seq_ID = ann.get("seq_ID", ".")
                source = ann.get("source", "extract")
                type_ = ann.get("type", ".")
                start = str(ann.get("start", "."))
                end = str(ann.get("end", "."))
                score = ann.get("score", ".")
                strand = ann.get("strand", ".")
                phase = ann.get("phase", ".")

                # 拼接 attributes 字段
                attrs = []
                if ann.get("id"):
                    attrs.append(f"ID={ann['id']}")
                if ann.get("parent"):
                    attrs.append(f"Parent={ann['parent']}")
                if ann.get("gene_id"):
                    attrs.append(f"gene_id={ann['gene_id']}")
                if ann.get("transcript_id"):
                    attrs.append(f"transcript_id={ann['transcript_id']}")
                if ann.get("rank") is not None:
                    attrs.append(f"rank={ann['rank']}")

                attributes = ";".join(attrs) if attrs else "."

                writer.writerow([seq_ID, source, type_, start, end, score, strand, phase, attributes])

        print(f"✅ GFF3 file successfully saved: {output_file}")


def find_final_candidates(candidate_data_summary, rows):
    final_candidates = []
    total_number = 0

    for summary in candidate_data_summary:
        region_info = summary["position_info"]["position"]
        region_name = summary["region_name"]
        # print("region_info", region_info)
        region_seq = region_info["seq_ID"]
        region_start = region_info["start"]
        region_end = region_info["end"]

        for row in rows:
            candidate_seq = row["seq_ID"]
            candidate_start = row["start"]
            candidate_end = row["end"]
            if candidate_seq == region_seq:
                if region_start <= candidate_start <= region_end and region_start <= candidate_end <= region_end:
                    candidate_info = {
                        "region_name": region_name,
                        "gene_id": row["locus_tag"],
                        "transcript_id": row["transcript_id"],
                        "gene_info": row
                    }
                    final_candidates.append(candidate_info)

    return final_candidates


def save_final_candidates(final_candidates, output_path):
    """
    Save a list of candidate dictionaries to an Excel file.
    Each dictionary becomes one row, and 'gene_info' (a nested dict)
    is converted into a "key:value" comma-separated string.
    """

    # Ensure the output directory exists
    os.makedirs(output_path, exist_ok=True)

    processed_data = []

    for item in final_candidates:
        # Copy the item to avoid modifying the original list
        row = item.copy()

        # Convert 'gene_info' dict to a readable "key:value" string
        if isinstance(row.get("gene_info"), dict):
            row["gene_info"] = ", ".join([f"{k}: {v}" for k, v in row["gene_info"].items()])

        processed_data.append(row)

    # Create a DataFrame (keys automatically become column headers)
    df = pd.DataFrame(processed_data)

    # Build the output file path
    output_file = os.path.join(output_path, "final_candidates.xlsx")

    # Save to Excel (requires 'openpyxl' installed)
    df.to_excel(output_file, index=False)

    print(f"Final candidate data (genes) saved to: {output_file}")


# file paths
# The modified csv file that contains the mRNA with depth of interests.
depth_path = "/lustre/BIF/nobackup/leng010/test/Asp_fumigatus/check_coverage/test3/all51_to_GCF_000002655.1_meandepth.txt"
# Reference genome annotation of the BAM file
gtf_path = "/lustre/BIF/nobackup/leng010/test/Asp_fumigatus/reference/GCF_000002655.1_ASM265v1_genomic.gtf"
gff_path = "/lustre/BIF/nobackup/leng010/test/Asp_fumigatus/reference/GCF_000002655.1_ASM265v1_genomic.gff"
# Reference genome assembly
ref_assembly = "/lustre/BIF/nobackup/leng010/test/Asp_fumigatus/reference/GCF_000002655.1_ASM265v1_genomic.fna"

# Bam file, used multiple genome assemblies align to the reference genome
bam_path = ("/lustre/BIF/nobackup/leng010/test/Asp_fumigatus/multi_align_2/40_assembly/normal_align/"
            "all51_to_GCF_000002655.1_asm5.sorted.bam")
# The txt file that includes the genome assemblies used in alignment.
assembly_path = "/lustre/BIF/nobackup/leng010/test/Asp_fumigatus/multi_align_2/40_assembly/accessions.txt"
# path to the 53 assemblies
assembly_dir = "/lustre/BIF/nobackup/leng010/test/Asp_fumigatus/multi_align_2/40_assembly/genome_assemblies"
# output base path of the fasta files
output_path = "/lustre/BIF/nobackup/leng010/test/Asp_fumigatus/check_coverage/test4_extract_seq/extract_allele"

# settings
up_num = 5
down_num = 5
assembly_num = 7
lower_limit = 10
upper_limit = 40
extend = 5000
reference_genome = "GCF_000002655.1"

# processes the input candidate mRNAs
# Input and output file paths
rows = process_data(depth_path)
filtered_rows = select_depth(rows, lower_limit, upper_limit)  # the candidate mRNA data
# print("filtered_rows",len(filtered_rows))
candidate_data_seq = extract_candidate_position_list(filtered_rows)
annotation_sorted = read_gff(gff_path, keep_type="mRNA")
annotation_sorted_dict = read_gff_dict(annotation_sorted)

annotation_info = annotation_sorted["NC_007195.1"]
for annotation in annotation_info:
    for key, value in annotation.items():
        print(f"{key}: {value}")
        continue

candidate_merge = merge_candidate_position(candidate_data_seq, annotation_sorted_dict)

# test the main code
candidate_data_test = dict(list(candidate_merge.items())[0:5])
# print(candidate_data_test)
candidate_data_summary = analyze_all_candidate_position(candidate_data_test, annotation_sorted, rows,
                                                        bam_path, assembly_path, up_num, down_num, lower_limit)

# gene_between = find_allele_sequence_inbetween(assembly_dir,candidate_data_summary,output_path, extend, assembly_num)

# find and save final candidate genes and related information
# final_candidates = find_final_candidates(candidate_data_summary, rows)
# save_final_candidates(final_candidates, output_path)

extract_reference_allele(candidate_data_summary, reference_genome, annotation_sorted, output_path, extend, ref_assembly)
# p
