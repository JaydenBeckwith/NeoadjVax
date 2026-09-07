#!/usr/bin/env python3
"""Strict adapters for SpliceAI VCF and strand-aware LeafCutter counts."""
import argparse
import csv
import gzip
import json
import math
import re
from pathlib import Path

FIELDS = ['CHROM', 'POS', 'REF', 'ALT', 'FILTER', 'ALLELE', 'SYMBOL',
          'DS_AG', 'DS_AL', 'DS_DG', 'DS_DL', 'DP_AG', 'DP_AL', 'DP_DG', 'DP_DL']
JUNCTION = re.compile(r'^([^:\s]+):([0-9]+):([0-9]+):(clu_[0-9]+_([+-]))$')


def text_open(path):
    return gzip.open(path, 'rt', encoding='utf-8') if str(path).endswith('.gz') else open(path, encoding='utf-8')


def spliceai_table(source, output, qc_path, threshold=0.5):
    if not math.isfinite(threshold) or not 0 <= threshold <= 1:
        raise ValueError('SpliceAI threshold must be in [0,1]')
    qc = dict(records=0, non_pass=0, missing_annotation=0, annotations=0,
              missing_scores=0, below_threshold=0, retained_annotations=0)
    info_header = False
    columns = False
    with text_open(source) as handle, open(output, 'w', newline='', encoding='utf-8') as out:
        writer = csv.writer(out, delimiter='\t', lineterminator='\n')
        writer.writerow(FIELDS)
        for line_no, line in enumerate(handle, 1):
            if line.startswith('##INFO=<ID=SpliceAI,'):
                info_header = True
            if line.startswith('#CHROM\t'):
                columns = True
            if line.startswith('#') or not line.strip():
                continue
            parts = line.rstrip('\r\n').split('\t')
            if not columns or len(parts) < 8:
                raise ValueError(f'Invalid VCF at line {line_no}')
            chrom, pos, _, ref, alts, _, filt, info = parts[:8]
            if not pos.isdigit() or int(pos) < 1:
                raise ValueError(f'Invalid VCF position at line {line_no}')
            qc['records'] += 1
            if filt != 'PASS':
                qc['non_pass'] += 1
                continue
            entries = [x[9:] for x in info.split(';') if x.startswith('SpliceAI=')]
            if not entries or entries == ['.']:
                qc['missing_annotation'] += 1
                continue
            if len(entries) != 1:
                raise ValueError(f'Repeated SpliceAI INFO at line {line_no}')
            for annotation in entries[0].split(','):
                values = annotation.split('|')
                if len(values) != 10 or values[0] not in alts.split(',') or not values[1]:
                    raise ValueError(f'Malformed SpliceAI annotation at line {line_no}')
                scores = []
                for score, offset in zip(values[2:6], values[6:10]):
                    if score == '.':
                        scores.append(None)
                        continue
                    value = float(score)
                    if not math.isfinite(value) or not 0 <= value <= 1 or not re.fullmatch(r'-?\d+', offset):
                        raise ValueError(f'Invalid SpliceAI score/offset at line {line_no}')
                    scores.append(value)
                qc['annotations'] += 1
                present = [v for v in scores if v is not None]
                if not present:
                    qc['missing_scores'] += 1
                elif max(present) < threshold:
                    qc['below_threshold'] += 1
                else:
                    # Bind each annotation to its ALT, never the entire VCF ALT list.
                    writer.writerow([chrom, pos, ref, values[0], filt] + values)
                    qc['retained_annotations'] += 1
    if not info_header or not columns:
        raise ValueError('Expected a standard VCF with SpliceAI INFO and #CHROM headers')
    Path(qc_path).write_text(json.dumps(qc, indent=2) + '\n', encoding='utf-8')
    return qc


def leafcutter_sample(source, sample, output, evidence, min_reads=3):
    if min_reads < 1:
        raise ValueError('Minimum RNA junction reads must be >= 1')
    with text_open(source) as handle, gzip.open(output, 'wt', encoding='utf-8', newline='') as selected, open(evidence, 'w', newline='', encoding='utf-8') as out:
        header = handle.readline().split()
        if len(header) < 2 or header[0] != 'chrom' or len(set(header)) != len(header):
            raise ValueError('Expected LeafCutter perind.counts header: chrom SAMPLE... (not perind_numers)')
        if not sample:
            if len(header) != 2:
                raise ValueError('Cohort counts require leafcutter_sample matching one exact header column')
            sample = header[1]
        if sample not in header[1:]:
            raise ValueError(f'LeafCutter sample {sample!r} is not in the counts header')
        index = header.index(sample)
        selected.write('chrom ' + sample + '\n')
        writer = csv.writer(out, delimiter='\t', lineterminator='\n')
        writer.writerow(['junc_id', 'rna_junction_reads', 'rna_cluster_reads', 'rna_cluster_ratio', 'leafcutter_cluster'])
        seen = set()
        for line_no, line in enumerate(handle, 2):
            parts = line.split()
            if not parts:
                continue
            if len(parts) != len(header):
                raise ValueError(f'Wrong LeafCutter column count at line {line_no}')
            match = JUNCTION.fullmatch(parts[0])
            if not match:
                raise ValueError(f'Expected strand-aware regtools LeafCutter junction at line {line_no}: {parts[0]}')
            chrom, start, end, cluster, strand = match.groups()
            if int(start) < 1 or int(end) <= int(start):
                raise ValueError(f'Invalid junction coordinates at line {line_no}')
            if not re.fullmatch(r'\d+/\d+', parts[index]):
                raise ValueError(f'Expected numerator/denominator count at line {line_no}')
            reads, total = map(int, parts[index].split('/'))
            if reads > total:
                raise ValueError(f'Junction count exceeds cluster count at line {line_no}')
            # Coordinates from leafcutter_cluster_regtools.py are already 1-based EXON bases.
            junc_id = f'{chrom}:{start}-{end}:{strand}'
            if junc_id in seen:
                raise ValueError(f'Duplicate junction in selected LeafCutter sample: {junc_id}')
            seen.add(junc_id)
            if reads >= min_reads:
                selected.write(parts[0] + ' ' + parts[index] + '\n')
                writer.writerow([junc_id, reads, total, reads / total, cluster])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='mode', required=True)
    vcf = sub.add_parser('spliceai')
    vcf.add_argument('--input', required=True)
    vcf.add_argument('--output', required=True)
    vcf.add_argument('--qc', required=True)
    vcf.add_argument('--threshold', type=float, default=0.5)
    leaf = sub.add_parser('leafcutter')
    leaf.add_argument('--input', required=True)
    leaf.add_argument('--sample', default='')
    leaf.add_argument('--output', required=True)
    leaf.add_argument('--evidence', required=True)
    leaf.add_argument('--min-reads', type=int, default=3)
    args = parser.parse_args()
    try:
        if args.mode == 'spliceai':
            spliceai_table(args.input, args.output, args.qc, args.threshold)
        else:
            leafcutter_sample(args.input, args.sample, args.output, args.evidence, args.min_reads)
    except (ValueError, OSError) as exc:
        parser.exit(1, f'Splicing input error: {exc}\n')


if __name__ == '__main__':
    main()
