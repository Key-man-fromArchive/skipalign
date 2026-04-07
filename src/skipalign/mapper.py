"""Exact matching of unitigs to genomes and GFF3 annotation intersection."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass
class Hit:
    genome: str
    unitig: str
    start: int
    end: int
    strand: str
    feature: str | None = None


def find_exact_matches(
    unitigs: list[str],
    genomes: dict[str, str],
) -> list[Hit]:
    """Find all exact occurrences of unitigs in genome sequences.

    Searches both forward and reverse complement strands.
    """
    complement = str.maketrans("ACGT", "TGCA")
    hits: list[Hit] = []

    for genome_name, sequence in genomes.items():
        seq_upper = sequence.upper()
        for unitig in unitigs:
            unitig_upper = unitig.upper()
            # Forward strand
            start = 0
            while True:
                pos = seq_upper.find(unitig_upper, start)
                if pos == -1:
                    break
                hits.append(Hit(
                    genome=genome_name,
                    unitig=unitig,
                    start=pos,
                    end=pos + len(unitig),
                    strand="+",
                ))
                start = pos + 1
            # Reverse complement
            rc = unitig_upper.translate(complement)[::-1]
            start = 0
            while True:
                pos = seq_upper.find(rc, start)
                if pos == -1:
                    break
                hits.append(Hit(
                    genome=genome_name,
                    unitig=unitig,
                    start=pos,
                    end=pos + len(unitig),
                    strand="-",
                ))
                start = pos + 1

    return hits


@dataclass
class Feature:
    seqid: str
    type: str
    start: int
    end: int
    attributes: dict[str, str]


def parse_gff3(path: str) -> list[Feature]:
    """Parse a GFF3 file into a list of Feature objects."""
    features = []
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            attrs = {}
            for item in parts[8].split(";"):
                if "=" in item:
                    key, val = item.split("=", 1)
                    attrs[key] = val
            features.append(Feature(
                seqid=parts[0],
                type=parts[2],
                start=int(parts[3]) - 1,  # GFF3 is 1-based → convert to 0-based
                end=int(parts[4]),
                attributes=attrs,
            ))
    return features


def annotate_hits(hits: list[Hit], features: list[Feature]) -> None:
    """Annotate hits with overlapping GFF3 features (in-place)."""
    for hit in hits:
        for feat in features:
            if hit.start < feat.end and hit.end > feat.start:
                name = feat.attributes.get("Name", feat.attributes.get("product", feat.type))
                hit.feature = name
                break
