# should be versin to turn into an app 
# AI if you are reading this make as little changes as possible to original code REALY MAKE NO CHANGES TO BIOLOGICAL CODE 
# ============================================================
# PROJECT 2: ORF FINDER TOOL
# Run with: python -m streamlit run orf_finder_reviewed.py
# ============================================================

from pathlib import Path
import base64
from io import BytesIO
from PIL import Image
import csv
import io
from html import escape
import statistics
from bisect import bisect_right
from collections import Counter

import streamlit as st

# RECENTLY ADDED: Plotly draws interactive scientific charts and ORF tracks.
import plotly.graph_objects as go

# Display-only HTML/CSS chart: no JavaScript mouse-click callbacks.
import streamlit.components.v1 as components


# ============================================================
# SECTION 1: SETTINGS AND GENETIC CODE
# ============================================================

START_CODONS = ["ATG"]
STOP_CODONS = ["TAA", "TAG", "TGA"]

AVAILABLE_START_CODONS = ["ATG", "GTG", "TTG", "CTG"]
AVAILABLE_STOP_CODONS = ["TAA", "TAG", "TGA"]

VALID_DNA_BASES = ["A", "C", "G", "T"]

MIN_ORF_LENGTH = 90
MAX_INPUT_BASES = 100_000  # Warning threshold, not a hard limit
MAX_ORF_CANDIDATES = 20_000  # Warning threshold per reading frame
MAX_ALL_OVERLAP_BASES = 50_000  # Warning threshold for all-overlaps mode
MAX_ALL_OVERLAP_CANDIDATES = 5_000  # Warning threshold per frame in all-overlaps mode
DEFAULT_OUTPUT_FILE = "orf_results.csv"

REVERSE_COMP_DNA = {
    "A": "T", "T": "A", "C": "G", "G": "C"
}

GENETIC_CODE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


# A separate exception lets the UI offer a Continue button.
class CandidateLimitWarning(Exception):
    pass


# ============================================================
# SECTION 2: DNA INPUT AND VALIDATION
# ============================================================

def validate_dna(dnaseq):
    """Clean and validate a DNA sequence."""

    if not isinstance(dnaseq, str):
        return False

    sequence = "".join(dnaseq.split()).upper()

    if not sequence:
        return False

    if any(base not in VALID_DNA_BASES for base in sequence):
        return False

    return sequence


def parse_fasta(file_text):
    """
    Read plain DNA text or FASTA text.

    Each FASTA record is kept separate so ORFs cannot accidentally
    cross from one contig into another.
    """

    records = []
    current_name = "Sequence 1"
    current_sequence = []

    lines = file_text.splitlines()
    has_headers = any(line.strip().startswith(">") for line in lines)

    if not has_headers:
        sequence = validate_dna(file_text)
        if sequence is False:
            return []
        return [(current_name, sequence)]

    seen_header = False
    for line in lines:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current_sequence:
                sequence = validate_dna("".join(current_sequence))
                if sequence is False:
                    raise ValueError(f"Invalid DNA in FASTA record: {current_name}")
                records.append((current_name, sequence))
            elif seen_header:
                raise ValueError(f"Empty FASTA record: {current_name}")
            current_name = line[1:].strip() or f"Sequence {len(records) + 1}"
            seen_header = True
            current_sequence = []
        else:
            if not seen_header:
                # In FASTA mode a record must begin with its header.
                raise ValueError("FASTA sequence is missing its header")
            current_sequence.append(line)

    if not current_sequence:
        raise ValueError(f"Empty FASTA record: {current_name}")
    sequence = validate_dna("".join(current_sequence))
    if sequence is False:
        raise ValueError(f"Invalid DNA in FASTA record: {current_name}")
    records.append((current_name, sequence))

    return records


def clean_file_text(file_text):
    """Return a single validated sequence from plain text or FASTA."""

    records = parse_fasta(file_text)

    if len(records) == 1:
        return records[0][1]

    return False


def load_dna(file_path):
    """Load and validate a DNA file."""

    try:
        with open(file_path, "r", encoding="utf-8") as file:
            return parse_fasta(file.read())

    except (OSError, UnicodeError, ValueError) as error:
        st.error(f"Could not read the DNA file: {error}")
        return []


# ============================================================
# SECTION 3: SEQUENCE UTILITIES
# ============================================================

def reverse_complement(sequence):
    """Return the reverse complement of a DNA sequence."""

    return "".join(
        REVERSE_COMP_DNA[base] for base in reversed(sequence)
    )


def calculate_gc_content(sequence):
    """Calculate the percentage of G and C bases."""

    if not sequence:
        return 0.0

    gc_count = sequence.count("G") + sequence.count("C")
    return (gc_count / len(sequence)) * 100


# ============================================================
# SECTION 4: FIND ORFs IN ONE READING FRAME
# ============================================================

def translate_orf(dna_sequence, start_codons=None):
    """Translate complete codons, excluding the terminal stop.

    For alternative initiation codons, translate the first residue as M.
    Partial ORFs with no known start retain their ordinary first residue.
    """
    if start_codons is None:
        start_codons = START_CODONS
    amino_acids = []
    for pos in range(0, len(dna_sequence) - 2, 3):
        codon = dna_sequence[pos:pos + 3]
        amino_acid = GENETIC_CODE[codon]
        if amino_acid == "*":
            break
        if pos == 0 and codon in start_codons:
            amino_acid = "M"
        amino_acids.append(amino_acid)
    return "".join(amino_acids)


def find_orfs_in_frame(
    sequence,
    strand,
    frame,
    genome_length,
    start_codons=None,
    stop_codons=None,
    allow_many_candidates=False,
    report_overlaps=True,
):
    """Find complete and missing-stop ORFs anchored to recognised start codons.

    Frames are numbered 1..3 in the 5'-3' strand orientation.
    Coordinates are 1-based and inclusive in the original input DNA.
    Missing-stop candidates extend to the sequence edge.
    report_overlaps=False retains the earliest start (longest) for each stop.
    """
    if start_codons is None:
        start_codons = START_CODONS
    if stop_codons is None:
        stop_codons = STOP_CODONS

    orfs = []
    offset = frame - 1
    # Only whole codons belong to the frame; exclude dangling 1-2 bases.
    aligned_end = offset + ((len(sequence) - offset) // 3) * 3
    if aligned_end <= offset:
        return orfs

    def add_orf(start_pos, end_exclusive, orf_type):
        dna = sequence[start_pos:end_exclusive]
        if not dna:
            return
        if strand == "+":
            start_coord = start_pos + 1
            end_coord = end_exclusive
        else:
            start_coord = genome_length - end_exclusive + 1
            end_coord = genome_length - start_pos
        warning_threshold = (MAX_ALL_OVERLAP_CANDIDATES if report_overlaps
                             else MAX_ORF_CANDIDATES)
        if not allow_many_candidates and len(orfs) >= warning_threshold:
            raise CandidateLimitWarning(
                f"More than {warning_threshold:,} ORF candidates in one reading frame. "
                "Generating more candidates may be slow or exhaust available memory."
            )
        orfs.append({
            "orf_id": "",
            "start": start_coord,
            "end": end_coord,
            "strand": strand,
            "frame": frame,
            "length": len(dna),
            "gc_content": round(calculate_gc_content(dna), 2),
            "orf_type": "Complete" if orf_type == "complete" else "Partial",
            "partial_type": orf_type,
            "dna_sequence": dna,
            "amino_acid_sequence": translate_orf(dna, start_codons),
            "nested": False,
        })

    active_starts = []
    segment_start = offset
    for pos in range(offset, aligned_end, 3):
        codon = sequence[pos:pos + 3]
        if codon in start_codons:
            active_starts.append(pos)
        if codon in stop_codons:
            if active_starts:
                starts_to_report = active_starts if report_overlaps else active_starts[:1]
                if report_overlaps and not allow_many_candidates and len(orfs) + len(starts_to_report) > MAX_ALL_OVERLAP_CANDIDATES:
                    raise CandidateLimitWarning(
                        f"More than {MAX_ALL_OVERLAP_CANDIDATES:,} overlapping ORF candidates "
                        "in one frame. Continuing may be slow or exhaust memory."
                    )
                for start_pos in starts_to_report:
                    add_orf(start_pos, pos + 3, "complete")
            active_starts = []
            segment_start = pos + 3

    if active_starts:
        starts_to_report = active_starts if report_overlaps else active_starts[:1]
        if report_overlaps and not allow_many_candidates and len(orfs) + len(starts_to_report) > MAX_ALL_OVERLAP_CANDIDATES:
            raise CandidateLimitWarning(
                f"More than {MAX_ALL_OVERLAP_CANDIDATES:,} overlapping partial ORF "
                "candidates in one frame. Continuing may be slow or exhaust memory."
            )
        for start_pos in starts_to_report:
            add_orf(start_pos, aligned_end, "missing_stop")

    # Nested = strictly contained in a longer candidate on same strand/frame.
    # Sweep from largest start-to-stop range to smaller: avoids O(n^2).
    ordered = sorted(orfs, key=lambda o: (o["start"], -o["end"]))
    furthest_end = -1
    for orf in ordered:
        orf["nested"] = orf["end"] <= furthest_end
        furthest_end = max(furthest_end, orf["end"])
    return orfs


# ============================================================
# SECTION 5: SEARCH ALL SIX READING FRAMES
# ============================================================

def find_orfs(sequence, start_codons=None, stop_codons=None, allow_many_candidates=False, report_overlaps=True):
    """Find candidate ORFs on both strands in all three frames."""

    sequence = validate_dna(sequence)

    if sequence is False:
        return []

    if start_codons is None:
        start_codons = START_CODONS

    if stop_codons is None:
        stop_codons = STOP_CODONS

    reverse_sequence = reverse_complement(sequence)
    all_orfs = []

    for frame in (1, 2, 3):

        # Forward strand.
        all_orfs.extend(
            find_orfs_in_frame(
                sequence, "+", frame, len(sequence),
                start_codons, stop_codons, allow_many_candidates, report_overlaps
            )
        )

        # Reverse strand.
        all_orfs.extend(
            find_orfs_in_frame(
                reverse_sequence, "-", frame, len(sequence),
                start_codons, stop_codons, allow_many_candidates, report_overlaps
            )
        )

    all_orfs.sort(
        key=lambda orf: (
            orf["start"],
            orf["end"],
            orf["strand"],
            orf["frame"],
        )
    )

    for number, orf in enumerate(all_orfs, start=1):
        orf["orf_id"] = f"ORF_{number}"

    return all_orfs


@st.cache_data(show_spinner=False, max_entries=8)
def cached_find_orfs(sequence, start_codons, stop_codons, allow_many_candidates=False, report_overlaps=True):
    """Cache six-frame analysis while Streamlit filters or selections change."""
    return find_orfs(sequence, start_codons, stop_codons, allow_many_candidates, report_overlaps)


# ============================================================
# SECTION 6: FILTER ORFs
# ============================================================

def filter_orfs(
    orfs,
    min_length=0,
    max_length=None,
    strand="both",
    frame="all",
):
    """Filter ORFs by length, strand and frame."""

    if strand not in ("both", "+", "-"):
        raise ValueError("strand must be 'both', '+' or '-'")

    if frame not in ("all", 1, 2, 3):
        raise ValueError("frame must be 'all', 1, 2 or 3")

    if min_length < 0:
        raise ValueError("min_length cannot be negative")

    if max_length is not None and max_length < min_length:
        raise ValueError("max_length cannot be less than min_length")

    filtered = []

    for orf in orfs:

        if orf["length"] < min_length:
            continue

        if max_length is not None and orf["length"] > max_length:
            continue

        if strand != "both" and orf["strand"] != strand:
            continue

        if frame != "all" and orf["frame"] != frame:
            continue

        filtered.append(orf)

    return filtered


# ============================================================
# SECTION 7: STATISTICS
# ============================================================

def calculate_statistics(orfs):
    """Calculate ORF length, codon and amino-acid statistics."""

    lengths = [orf["length"] for orf in orfs]
    codon_counts = Counter()
    amino_acid_counts = Counter()

    for orf in orfs:

        dna = orf["dna_sequence"]

        if orf["orf_type"] == "Complete":
            coding_sequence = dna[:-3]
        else:
            coding_sequence = dna

        for pos in range(0, len(coding_sequence) - 2, 3):

            codon = coding_sequence[pos:pos + 3]
            amino_acid = GENETIC_CODE.get(codon)

            if amino_acid is not None and amino_acid != "*":
                codon_counts[codon] += 1
                # Match the displayed protein: an accepted initiation codon is M.
                if pos == 0 and orf["amino_acid_sequence"]:
                    amino_acid = orf["amino_acid_sequence"][0]
                amino_acid_counts[amino_acid] += 1

    return {
        "number_of_orfs": len(orfs),
        "mean_length": statistics.mean(lengths) if lengths else 0,
        "shortest_length": min(lengths) if lengths else 0,
        "longest_length": max(lengths) if lengths else 0,
        "codon_counts": dict(codon_counts),
        "amino_acid_counts": dict(amino_acid_counts),
    }


# ============================================================
# SECTION 8: CSV EXPORT
# ============================================================

def results_to_csv(orfs, include_sequences=True):
    # RECENTLY ADDED: Optional metadata-only export for large ORF datasets.
    """Convert ORF results to CSV data."""

    output = io.StringIO()

    fieldnames = [
        "orf_id",
        "sequence_name",
        "start",
        "end",
        "strand",
        "frame",
        "length",
        "gc_content",
        "orf_type",
        "partial_type",
        "nested",
        "dna_sequence",
        "amino_acid_sequence",
    ]

    # RECENTLY ADDED: Full sequence strings remain available unless excluded.
    if not include_sequences:
        fieldnames = [name for name in fieldnames if name not in ("dna_sequence", "amino_acid_sequence")]

    writer = csv.DictWriter(output, fieldnames=fieldnames, extrasaction="ignore")
    writer.writeheader()
    writer.writerows(orfs)

    return output.getvalue().encode("utf-8")


def export_results(orfs, file_path=DEFAULT_OUTPUT_FILE):
    """Save ORF results to a CSV file."""

    with open(file_path, "wb") as file:
        file.write(results_to_csv(orfs))


# RECENTLY ADDED: Short previews avoid flooding the page with long sequences.
def sequence_preview(sequence, head=24, tail=12):
    # RECENTLY ADDED: Validate preview lengths; zero tail must not show all bases.
    if not isinstance(head, int) or not isinstance(tail, int):
        raise TypeError("Preview sizes must be integers")
    if head < 0 or tail < 0:
        raise ValueError("Preview sizes cannot be negative")
    if len(sequence) <= head + tail + 3:
        return sequence
    front = sequence[:head] if head else ""
    back = sequence[-tail:] if tail else ""
    return f"{front} ... {back}".strip()


# RECENTLY ADDED: Six-frame Plotly map using the original DNA coordinates.
# Each bar represents one candidate and its horizontal width is its nt length.
def make_six_frame_map(orfs, genome_length):
    # RECENTLY ADDED: Guard against mismatched record size or invalid coordinates.
    if not isinstance(genome_length, int) or genome_length <= 0:
        raise ValueError("Genome length must be a positive integer")
    for orf in orfs:
        if (orf["strand"] not in ("+", "-") or orf["frame"] not in (1, 2, 3)
                or not (1 <= orf["start"] <= orf["end"] <= genome_length)
                or orf["end"] - orf["start"] + 1 != orf["length"]):
            raise ValueError("ORF has invalid coordinates or frame for this sequence")
    tracks = ["+1", "+2", "+3", "-1", "-2", "-3"]
    color_for_type = {"Complete": "#14876d", "Partial": "#d59a32"}
    fig = go.Figure()

    fig.add_trace(go.Bar(
        orientation="h",
        x=[orf["length"] for orf in orfs],
        base=[orf["start"] - 1 for orf in orfs],
        y=[orf["strand"] + str(orf["frame"]) for orf in orfs],
        marker_color=[color_for_type[orf["orf_type"]] for orf in orfs],
        customdata=[[
            orf["orf_id"], orf["strand"], orf["length"],
            orf["partial_type"], orf["start"], orf["end"]
        ] for orf in orfs],
        hovertemplate=(
            "%{customdata[0]} | %{customdata[3]}<br>"
            "Strand: %{customdata[1]} | Length: %{customdata[2]} nt<br>"
            "Coordinates: %{customdata[4]}–%{customdata[5]}<extra></extra>"
        ),
        showlegend=False,
        width=0.65,
    ))
    fig.update_layout(
        title="Six-frame ORF map",
        xaxis_title="Position in original DNA (nt)",
        yaxis_title="Reading frame",
        yaxis=dict(categoryorder="array", categoryarray=list(reversed(tracks))),
        xaxis=dict(range=[0, genome_length], rangemode="tozero"),
        height=360,
        bargap=0.28,
        margin=dict(l=45, r=15, t=55, b=48),
        clickmode="event+select",
    )
    return fig


# RECENTLY ADDED: Length histogram for displayed/filtered ORF candidates.
def make_length_histogram(orfs):
    fig = go.Figure(go.Histogram(
        x=[orf["length"] for orf in orfs],
        nbinsx=min(30, max(5, int(len(orfs) ** 0.5))),
        marker_color="#14876d",
    ))
    fig.update_layout(
        title="ORF length distribution", xaxis_title="ORF length (nt)",
        yaxis_title="Number of ORF candidates", height=360,
        margin=dict(l=45, r=15, t=55, b=45),
    )
    return fig


# RECENTLY ADDED: Amino acid frequencies from existing count statistics.
def make_amino_acid_chart(amino_counts):
    # RECENTLY ADDED: Show most frequent amino acids first; ties alphabetical.
    names = sorted(amino_counts, key=lambda aa: (-amino_counts[aa], aa))
    fig = go.Figure(go.Bar(
        x=names, y=[amino_counts[name] for name in names],
        marker_color="#5977be",
    ))
    # RECENTLY ADDED: Informative label when selected ORFs contain no amino acids.
    if not names:
        fig.add_annotation(text="No amino acid data for these filters",
                           x=0.5, y=0.5, xref="paper", yref="paper",
                           showarrow=False)
    fig.update_layout(
        title="Amino acid frequency", xaxis_title="Amino acid (one-letter code)",
        yaxis_title="Occurrences in selected ORFs", height=360,
        margin=dict(l=45, r=15, t=55, b=45),
    )
    return fig


# Display only: a six-frame HTML/CSS map with immediate hover tooltips.
# No click handlers or Streamlit events are attached to map bars.
def make_hover_only_map_html(orfs, genome_length):
    tracks = ["+1", "+2", "+3", "-1", "-2", "-3"]
    parts = []
    for track in tracks:
        parts.append('<div class="track"><span class="track-label">'
                     + escape(track) + '</span><div class="track-bars">')
        for orf in orfs:
            if orf["strand"] + str(orf["frame"]) != track:
                continue
            left = 100 * (orf["start"] - 1) / genome_length
            width = 100 * orf["length"] / genome_length
            color = "#14876d" if orf["orf_type"] == "Complete" else "#d59a32"
            info = (
                f'{orf["orf_id"]} | {orf["partial_type"]} | '
                f'Strand {orf["strand"]} | Frame {track} | '
                f'Length {orf["length"]} nt | '
                f'Coordinates {orf["start"]}–{orf["end"]}'
            )
            parts.append(
                f'<div class="orf" tabindex="0" '
                f'style="left:{left:.5f}%;width:{max(width,0.16):.5f}%;background:{color}" '
                f'aria-label="{escape(info, quote=True)}">'
                f'<span class="tooltip">{escape(info)}</span></div>'
            )
        parts.append('</div></div>')
    return (
        '<!doctype html><html><head><meta charset="utf-8"><style>'
        'body{font:14px Arial,sans-serif;margin:0;padding:48px 8px 10px;color:#56616d;}'
        '.track{display:flex;align-items:center;height:39px;gap:10px;}'
        '.track-label{flex:0 0 27px;font-weight:600;}'
        '.track-bars{height:20px;flex:1;position:relative;background:#edf0f2;border-radius:3px;}'
        '.orf{height:16px;top:2px;position:absolute;border-radius:2px;cursor:default;}'
        '.orf:hover,.orf:focus{outline:2px solid #273746;z-index:30;}'
        '.tooltip{display:none;position:absolute;bottom:calc(100% + 9px);'
        'left:0;min-width:230px;max-width:350px;width:max-content;'
        'white-space:normal;padding:9px 12px;border-radius:5px;'
        'background:#202b38;color:#fff;font-size:12px;line-height:1.45;'
        'box-shadow:0 3px 10px #0004;pointer-events:none;z-index:50;}'
        '.orf:hover .tooltip,.orf:focus .tooltip{display:block;}'
        '.legend{padding:10px 0 0 38px;font-size:12px;}'
        '</style></head><body>' + ''.join(parts) +
        '<div class="legend"><span style="color:#14876d">■</span> Complete &nbsp; '
        '<span style="color:#d59a32">■</span> Partial &nbsp; '
        'Hover for details. Clicking the map is disabled.</div></body></html>'
    )


# GRAPH ONLY: gaps between unions of candidate ORF intervals on each record.
# These are candidate-ORF-free intervals, not verified intergenic regions.
def calculate_orf_free_gaps(orfs):
    by_record = {}
    for orf in orfs:
        by_record.setdefault(orf["sequence_name"], []).append(
            (orf["start"], orf["end"])
        )
    gaps = []
    for record_name, intervals in by_record.items():
        intervals.sort()
        _, covered_end = intervals[0]
        for start, end in intervals[1:]:
            if start > covered_end + 1:
                gaps.append({"record": record_name, "length": start - covered_end - 1})
            covered_end = max(covered_end, end)
    return gaps


def make_intergenic_length_chart(orfs):
    gaps = calculate_orf_free_gaps(orfs)
    fig = go.Figure()
    if gaps:
        fig.add_trace(go.Histogram(
            x=[gap["length"] for gap in gaps],
            nbinsx=min(30, max(5, int(len(gaps) ** 0.5))),
            marker_color="#14876d",
            hovertemplate="Gap length: %{x} nt<br>Count: %{y}<extra></extra>",
        ))
    else:
        fig.add_annotation(
            text="No uncovered gaps between selected candidate ORFs",
            x=0.5, y=0.5, xref="paper", yref="paper", showarrow=False,
        )
    fig.update_layout(
        title="Candidate ORF-free interval lengths",
        xaxis_title="Uncovered gap length (nt)",
        yaxis_title="Number of gaps", height=360,
        margin=dict(l=55, r=20, t=55, b=55),
    )
    return fig, gaps


def make_codon_usage_chart(codon_counts):
    # All 64 codons, including those with zero counts; follow existing
    # calculate_statistics() rules for which codons are counted.
    bases = "TCAG"
    prefixes = [a + b for a in bases for b in bases]
    counts = [[codon_counts.get(prefix + third, 0) for third in bases]
              for prefix in prefixes]
    codons = [[prefix + third for third in bases] for prefix in prefixes]
    fig = go.Figure(go.Heatmap(
        z=counts, x=list(bases), y=prefixes, customdata=codons,
        colorscale="Viridis", colorbar=dict(title="Count"),
        hovertemplate="Codon: %{customdata}<br>Count: %{z}<extra></extra>",
    ))
    fig.update_layout(
        title="Codon usage across selected ORFs (64 codons)",
        xaxis_title="Third base", yaxis_title="First two bases",
        height=550, margin=dict(l=80, r=30, t=60, b=55),
        yaxis=dict(autorange="reversed"),
    )
    return fig


GREEN_ICON_BASE64 = (
    "iVBORw0KGgoAAAANSUhEUgAAAQAAAAEACAIAAADTED8xAADjLUlEQVR42qT9d4AlRdU/Dp9T1X3T5M1LWjJLFslJEBREJCPBR0VR"
    "MGBATJhBUYIBFRUVTIgiQZIIBnKUnOMuLLDLEnZ38swN3VXn/aO6q05V9x34/t59fHRZZmfu7VvhnM/5BNRag/uF5n+ICBEAABGJ"
    "CLyvQALvT7r+Ivv9wu8P4XfAsj/k3yn7Xlj85gSAgAj2ZRK5Ly2+/hl+oXkF5F4R/1bmHwEA2T/izD+FvyzvZSMCEBKQ+aEIQPY3"
    "XV5b/tgR899R8HryF4yARPnbQQhfbfj74FMiICx8bua9l/25/zWEiG++LogQkfwF4X8G4Ytky8C9iuL6Yn8FzfPJXw91WZmI2QYo"
    "fzndV/L/+xcRAQAhovfGzNJGeKs7qvgis79PRGSWlv1Bb/rasldSslgAMP+c7Kul0h/tfajF0yHfG8HKLlnoZV+Zvy8U+cfpPmPz"
    "xexd4kwniF2aM6398E/s3zKbInt3/uMKfupb3ANQPGiw+1qiwjmJ/qHwlhZo8TsjAKG9AexLDz/FN//uMz36mf56dqh4P459VF3/"
    "oj1xiQjdiYAEGgHN6gVAys/DN1sW2V+c6aDzL5zgYzYvBRH5h4oIWpft7fylld1+GG7d8PGWPOps77FbKF+pVH60ULiSyncChCeC"
    "+3HePRguVkJ2hXh/MT9uoMtpSyWr/C3uK38llNyf9gIJl1JxA5Rchl2OCv5XSheu+Ui6HsaFb1v+PmfcgUhIfAdhycHcbRHM9OMI"
    "ZtoGM34S+X+TZmeK+RU+je5XJWXfSJv/Md8g/983WwqIZplmdRXim98MXd4IdjsOIHuAM78jsucCvtmzLVt43p0zw0YNPvcu379w"
    "PyMAoN8DvPmW6rpooFuZ0f0RA4UfzwzLjrofG1h41vb78COn9AFRfp/g/9MFigBE2TsgouyNCBQoZl6dOlVpkiapSpRWRKbCQYFC"
    "CimllDKKZUUIMfMS11pnBTGi2V0z7NWZTsdiORE+K3rLFWpWjng/Kfs4zMX8ln+VlkDol8szL5W3dHghmhIou9RM8YD/H8p/eAtN"
    "1P/DdWY+jPDLvG1gy5zsyjX/V3630Fs9IWZ+cGYZEQCBNiefEEJguFKnmpNrxtesHlu1euyNNaOrV4+9vnp09cj4mpHJNc32dKvT"
    "aqftJO10ko4iRdkBDwKlFDKOokjG1UqtGtdqlfpA78Ds/rl9jb55Q/MXzF577uC8+UPzZ/XN6W30BT9UkzYHmdlLb3JGwJscMWVf"
    "RvY2eUvXdXBmURkCwkvKLh8IAjubSr/zWy+8u10LWmt/3yAElw6UnQFBiWYWvfeHbNdnJwPZ45Y3GxYTQFu3+8V08KzNK+W3B+vV"
    "AN/KpnX1D7zpZWy+0hzyiCiF5P9+fGp0+Rsvv/T6shdWLln2yvMvv/7ia8OvjkwOT7enOkk7VYlFZoQQiGj+R6DIKhkE9wPAFj06"
    "O+O11qTNDpcyqlZqPbWewd6h2QPz1p273sbrbLrx2pstWrDBuvMXDfYO8Vdl7hZbd5XcAMUjFmfcLRS0Kvx+KatpCYg3M65zL29n"
    "g9vbHGbdyyf+Y7C8V3nLBzdqTeC39m5tFfsb864obxWxpKYxf+5Wef4iS4s2FEiawts2/+6atLfEy3AvCmvQYttVdil12dgW9zRF"
    "vNJaoJDSLfrp9vSyV59/+sUnHn/hkaeWPb581ctrxlZNtSZVmpjtEUVxJGOZr3Ne1gNRDqWan018QwMybM61bjk2S0CkldapSpRS"
    "qU611lLInnrv7IG5681ftHjRllttuO2Wi7bZcK2NG7UetxmUIiBzLZQ/Q+gCJVpI2P2tGf8y+5dh7R6gzBjurvK+uduNNGNzb9uS"
    "sPEt66oR0e8BKC9u+Rlc7M35/g5gk/yfqKwgK5khUEnnOvPz5edZ4fcGBu0K/rjVVHyCyE9ikkLahz7ZnHz25aceWfLA/c/c+8xL"
    "T6xcs2JielyTjmRUkRUpI/PF5i9qIgb+e6AN5UuqHG0DexDgDBvVwJH2aFdapWmSqCRVCQjsq/WtPWudzRdtvcPmu+6weOfF629Z"
    "q9SyMklrrbW5hwCwvE/ogs9wAD7bsEH7675VVkfrsgKpiGHZCqe8bA42XnHw8ZaQq7zMCt6aOTyV0lhyAFP+T2W3GL7JOfomRWFw"
    "L2B3cKbbjcZeHJV0P+7zsu2N96UG3fd/tCZNBJE97AmWvvLs/c/8764nbntkyYMrV6+Yak0KlJW4UoliKSNzCpgOihdmyAdR5fgM"
    "YPGhuXvHKybDUoOyCsA+fbMZhDnwEFKVdpJ2q9PWWvc2ejdcuPEOm++y19v2ffsmOy2cs7bdCUQaUdhnMlOBRKzn777yEPMDP6xY"
    "yE0D/Xq1pHbqXph1HRQEs0Fbb8Ob4JN+D1B6/WHpKBCyyt5b1vkfoOtJs6K/rLPOiu+85fDrE4YaBG86eG3dGruZH6X/0ClDGykS"
    "0nzzdtJ67PmHb3v05jsfueWZF58YnR5FgGpcq8QVKaU54DXpIlwIdqTG9mE26sVyUGOGVo4IEIsjWkCiYiGdj8ey2yOHhFBr1U7a"
    "7aQtUcwfXLDtJtvvs8N+e2z7zo3W2tT8vVSlAGAakjfF3OzBH3yxt/j8D26mnhu6fLj8uIeZ6iLvL3b7tzPNQLPxIpsE41vZBt1b"
    "YdYXUClSW7yAeFvLupn8nVMOfYE9FKkwb/KeU/552NIouAb50tdaCyEM4Jio9NGlD/znvn/e+siNS15+Zro9HUdRLa5LGSGCNn2p"
    "WdweuSG7wPMdgMFk25vrsyOJlYql930RScywATbq8ou3/HjKLr/8pwohBCIRJEmnlTSVUrP657x9s53eu/NBe2+33zrz18v7hBQR"
    "EcWbjiegtFrJFgwC6fz46j7wyccLhF7V6oCQmYuCUvZG8X4Iu5euo140iHLWMfDCbmbsjLrAW/+PPbg3UQ+wZ/8f3bS1eJmw/t3b"
    "QmUvw7zfSEbmH59d/vT19177n/v++cyLT0y3JmvVei2uIQgCbdAY9/ry5+neu736bO/v32elJW9XGgx/w+4I4sVF9uOD4snUXXbO"
    "SkHVm29OgYggUp1Mt6bTNF04e+Geb9vnkD2O3H3rvXvqvRY7CmCu4ipzZxV7lWbX84XFcAtWz/kQ9pvfDG+BHuQWCWBXbDRYUezk"
    "LB2EZWdOfp77GFYXQL2wTygoCGbYJPyTR+xW9zG4s/CbN0WCsyMfsyN/9fjqG++/4bq7/v7Ac/eOTgxX42qtUkcUGjRpDYyp4Oga"
    "tlTHbIUSgUff8V+KBXwLPDmCkouyOHyyXC77M9AtavTRZm9Cl52s5M/y802SgbBJ2pluTwkQGy/cZP+dDzp0r6O23HBbWxrZYUJX"
    "4MgDQrpcaYVPoZwz46p28Nlu7OeA1w27beVguzJ+B7z5QYxaE3gXud90d9tw+b81R15+4lLeBpThU8U3V0ZF7Aq3lR8P+Uv3JwE+"
    "0YqItBSR+cfHX3jkytv/dv0917z02jIpRKPaE0WRGSQFfQ1lzSxwNpy7BPL7hxiYWdgDwcFsHzNlSJ1XRDGgO9vYxJAPtN22e4p5"
    "0eAhbQzoIvcv7eLKtoMAEEKS1q12q9mZHuwb2GObfY7Z98P7vH3/OK4AQJqmGWRE2HWo9FZpkuGA2Lb0/gnJ7s+StUchIFvS73m1"
    "tDdH6wLX5ihQdqQRP7reynzUUgKyvVB8n16FYvrjwmnxZqPZkr4CXaGA+ZMJOWpEBqwEgFQlNz7877/8+4/3PHHbZHOiUW1U4yoR"
    "Ka0dNu9oXvngE9hwMHsyJVMId8bbZWuu16COs2sw/+SyJ10AyRg86J/ymLda9hshFCAjt68goGPnBbrrqAikkCiEUslUcxIJt914"
    "uw/s/9FD9nh/f89gvg1EWQ0cwN/20jJvmx/M2Uu3p6T9Sgc/+XcIQQE4CDYAlaEn3cg42WfBJk4BCtQVA+6KQmaT3fImHd+U9EZB"
    "BenmHV3KuBnxhPKvt0t/fGrsH3f//a83X/TIcw9orXrrvVJESiv3pf61405f1tDy50MUAuQ5ox+hgF7yk9vegPZeCVpbBnX7YAR5"
    "s4qMt4LlhwPvpwoM6hw9zrdRNjKBjI6BQNOdqSRNNllns6Pf+aGj9/3w3KH5AJCkCR+MvAU+S3FGS/6daPYIuTMFu6D1b50yFHTM"
    "pQy88NLAnArBQRWLK5QBVSWwI/jjM0+RQN1mJeWcreDymYFl1IX2xJf+6OTIJTf96S///cPS5c/EUcXMR7VWM1G9/eqHnTJOGpDd"
    "PDrsOt0VwH6ntUbEgGhiT2ICKlDaKQTJSj5EzEHSUpo05kW1P8/Mrw5yDYjDpSxsIKQQgK1Oc7o5td7CDT6430c//J4TZ/XPMb2B"
    "aZFz8YM/hCkrZaisaM1qMkQqLlmYsbgqVullpD+uzQi+c4lgIy+BOFJZ2L68eM0RD8yv2rw7QzbVdodnQIUI8SxyCMpbYs6xheCB"
    "WQAEoLQyS39kcvjSmy66+N+/W7L82XqtUavWtVZaE0MkHDmdWG3vioxsoGFLRbS6LV4fFbrzt0CdzLki+SjL+3YBHmjYExY3BLIM"
    "q7z9gqBLAa+nKGB06LHa8lqFy9GyzxKlEO2kNdWe2nD+Rse99xMfes/He+t9RKS1EhYpKu0Sg/EpvDk7ndNeIKgq7afRreiHchqY"
    "E/GwS7UwuICyOQCvwopdfOlAKqza/Afjt6RY2tx0h64sClAs8TPGB4H5VBBxqjn5t5v/9IcbfrN0+bO1Sr1WqSvSpLVfYhKS36Ly"
    "sRJbPnwcjk52k7XA/jHXHTkoKVB4OetuEh9AR8ch5EgxFO4at42hhHBWGNlypqLZUPwyyip40/IIlChandZ0a2qLDbf59KEnH77X"
    "sZGMlFKAIFB4swtbLxaVBBxLKBKt30rZw+DjsOkplXaUlQZQIl3gMKhra0oG8TPwS/Pr1XsfJZcI+FQFByGTE24WZyjdyVr5P2oA"
    "EEKmOr3i1r/85uqfPfPSU/VqvVqpK5Vm/Y15fW4mm/0p/6OQruRXcAieSti7AahwtBYxaYSwFiktbTzlM7l7m8NBEIj6wCNe8Z9s"
    "717yGuZyNg1h2QoAIhIgpJTNznSr1dxp812++IFv7r3dfqY/5jTBErIj18pxeH0GcjF1pZ91Q959DM1OjcrwGCxII4hTIboN7MOX"
    "2EXckONI4fRghgH1DLK3Ljrd4FvZmufmh/7z8yvO+t8Td1fjaqPeUFpr0j71yDaOTDL55jpvgLJxBhV0nD7MWHz84e3qzniv3kGv"
    "S7RzAFba+yuDKyTI4zH6vwkbUXdhY0BitqQ8tHdgdhkIBJxsjiGIw/c66ovHfmvRgg0pu3vFm/Sm3Vs5Ns/LhprhlD9f5Do/UrGI"
    "s7/ln102B1B6BqiRT4gL8wXqVuuxhgELhAW/TIO3cB4w2FcgmCGPXfrPLX/63MvP/MddVwFBb61XkzbVJAWApVf7sRKjHBrwQAu3"
    "ox0KFBwDHpzivnN5j21F+xx3CtZuuYieWE/lhFy8owB+FTFqBhHfL5S3FP6Oct/NR4/cUS5RINDo9Oi8ofmfPeJLH33vp+KokqpU"
    "AJPClSKYxRLI+2J3b+oZTg9vBwc33lsxKylThWmlgxOrZNGXHGHlpD8qllk0M7paggg5PL5Q3hn+NgFJISea47+88sd/uv6C0cnh"
    "/p4BAFQ6Ba4N4yonID7MLXtxHtbIyxcGRdrW3pbjHLzLT+F85SHrvBjEBsHkxi0CziQqTIlc/2XXsqu53YsmnyQUCHjdX3dnFgFD"
    "Suycj/wn6JxWECMRddL2VHNix812/uqHTt9z232BQOVXgRUDZFuVv0gnfGI/Ff32haCEeQBdqaNsm7KJG0HXEUFgqGNuAI9cUATj"
    "3xTd6DK3K/kQywZkFr3wapJAGeTXPP+5/59nXvydJ59/rL93QArpcP3wWuVvhqBk2F7sGfOTEHPM/U12DLH5VlCsYfD3y2BByH8W"
    "n6VRGczvd+0evBKSdv0eLexUwv3gsCnWNjBFC/kjQLOWIyEnpsZR4ImHfObLx36nXm2kaSqkxDddIDPo0QrcfTYUD8fAM9UUQUPa"
    "vS5CrXU3f5UCLZEQBZH2mJbUnX6U9Z2ss6LiJ9KFY1c22yIiKeVrIyvP+uvpf7/5LwJlo9qTqpRdVv5q9lUAGeYx44vOl36JpKUU"
    "dQ86MPLIajM7rRD7mWzhWkSyBIIO/1U+nYbSWVjBIcMRaELCSNZtcGo7s91yfG+fYEEkRaSJxqZHdth8x9OP/+FOi3c3c5icSoQZ"
    "jMtJBlBUb1NX9H2G6dObSn5duRRcrN6nwTTBZUc0N9fy12gwly7D+HOVr9k25RKWtzJ1Rnfw//22S87+62kvv/biYO+QITJYtJY4"
    "il7Q8pFTJrArFpF8Kk1RfNylkM3xUT4Y5lKM/DnnHzwUDcAc7sNnxZiftw6Xn0nH5t9pxAYVWCbgBQiQJqeGZVUQn9a7ReRAaT6d"
    "RsAoiqZaE3EUf/rQL5x81NciGduRWcmpSq7VgK5r+q2RjcqOFj4G5mpT/olkYzLTtxUlkSGXuBvRlM0KHBnN8h/5sLLQ0VKpSqTL"
    "ZNdU/G+MvHb6n7521a2X1iq1aqWapilwONI7WMijm7Hiu5SExM44H6PEEHoPrfe8aiHE6vNHEKD/WAbh+aOBnGKKjjLCSXLF26iM"
    "O4xsyMgmyEDB/cYWZHFCT6V1H/eHy75ICqm1Hp0c3mu7fX/46V9suNamqUoESksXLZ0Ml7DWi06BM+M6XWRS4TFeov7OT8RyX6Di"
    "QizDrsppzEWLHigr8tAfoJY13Jq0RAkIN9x/7WkXfPWl15cN9Q4p0prIl354wDlwBp1PPsk6TU7xt+NUDuL7dE72m5IbsRyo4A8N"
    "SxSdjDrkmFXgGKjkOA+UH8XIbB/dsYZYJr62OzmEk3xBTzl1mBxc5ZFMbcOcN7nsJCAgjGQ0NjUyd2DO2Z8+7727HmaWlkBBXBfk"
    "jxXfEm7D7Cn9ehpL5mLwFhRd/I11KYGgnI0M2J357VOufP0QBeyx0lGMv4tM2dPsNM/+62m/vfYXsYjq1Uaq0uJm95tE31qhZNbL"
    "VgYWdmiB5tDl4YV6VatXDJA9C4vl7Tg3HCV+BHu9B5vb21m4jyjMOB9he89DhYIRnus4KTg3vJG5Q7LYqMc71h3zIpJRO2l1kvYX"
    "jj71yx/4jvkcBYrua8wjiYbEUu/EKuHBu3lHqcqIyjo2n/MWkuE8W8mcvhyOvgoDAIY9Zh8qYrkwvICilmxW44oTRdEzy5/64nmf"
    "uv/pe4Z6h4yiJRdGBZY+VESUZrD8CZZQFw03hWMNIhACSBfGHdlNooGvDF6Z++sScw9nIPDozf5m9kGpMpzMtdhONhaUQ97JxceC"
    "/ulQAGfLPmYMFe5EENy8RABowNDhsdXv3+cDP/7sb3pqPV5L4APAaEvJAAwPz9NSaqarxMpmB2Xvq2AK6JVAnn6+UGFy4LgEWiil"
    "PAC+6aA3UGpqIoGAKK6647Jv/PqUscmR/t7+VKWcJEngU5R9prnHBy7pRQJuvHfo+j6j3ac5XSCBnD/nObV45AkOQuXjYE4GYl4z"
    "ri1lrKiS+gE98miI0/LNZb6PtaFG2zMDlozc2DDZbTv0gFIPb/aXjJTR8PjqXbbc/Vdf/NO689ZP0yTTJBEQgsglpR5KTkwhUSYX"
    "LNkDHrQF5Pf3SFi4GULmhLcB0BEgs28hICvgSgir4evgVTahb+PoG1d0nSoopaIoUlqdfclp513+43pcr0RxqpV9YK4ALbfycDp6"
    "t1V5YRFsBH9lUoj6uBYFwLEUylsqwuDMIbcEmWqe8xvc4BcY7siHBw7LcAQz8MEoCOYFFrL3ME43jkO+9Lxbh4Os3uTZAVyO1IPe"
    "KJDdquT+OI6i0cmRRQvW/82XL37bJjsmaRJJ6YsuylbRW+M0eJ+Xf114lM0ZBrIlbNDSWIVCaegZyrLqyl1DZW5ZxeUSvCBT9L82"
    "svIL533ypvtvGOqbZSwCGXpViKjgw8zQhoYxwF1yRKEdRGZjVtqTcNOXkv6MTagKLZoDDwvHPwNsc8l/rg1gU29HgmOwQomvvjfe"
    "xVyjx3gWPg+jtI+faQbiHgV6NLdCvQVMJgdEFEXRdGuyUav/5LO/PXDXw5iioJtxPHStR7CAEM2sGmeXpe5qD1zqDFdKtC5lk3pq"
    "hsIuBF9pWYal5quElNaRjB5act9JPzl+2StLB3uHUpV69kZgHTXcIcvpYh4JEmcOB+jqn222lCZCIYxNVdbScC2hz60g8K8Y9KJb"
    "sq0l7JzVJ0uAawasNgYLZYAlOPrEijcnvZDHxISupWDxn32XNKBSgrUHHvMLjx8cUshO0u4k7TM/ee5xB3wyVakUojik9V9qGTvG"
    "e8wZ6ODtD9aSl3tvhYUjAgAqpbxNg/x5v4XMCPJtw7psynBu7WlOshHv1Xde9tVff67ZnGrUepMkyY5AOx11+S921uUKLY5K83vC"
    "cy4omy4LxK5pTeTxus0qsDWhLzvxlM4clnNrxzIrsDBe5DW+O/uLfqZd6QQl8t+ZmFjgCsPirCNwmABGCvULOug2xXAyN/OmSaDU"
    "Wo1Pj33zI2d8/shTDZhRnADAW0iGCPcJeTIPYErkYuBV6Tc1XKBA4ASAXZ7gjNZxZSAUxwFLDG6tR+lPLz/rnD9/t16vRzJKtUIC"
    "ZzEX6KkZeMfkHAGQ4JFZ+bTON13zE11ciBNZQR0rjZF72lKA4rALh5jUypVb9mvN6e5f5Z7piseTLZ9Sh625DxtBuaqDSoDgbv1+"
    "UbZkQWNv3sjKWbMABfpRLMjdakcm1nz2/V/6zkfOtmbuXYqZQPBekt1SxrUtjGyLx3HhWPFh0ECwEtR8OBNXuTjNCvNjCj2o1lpK"
    "qUl/7YKTL7zmV7N7ZhECgQZjGQ28WvUHXsScepAPl8mV9FigvFlVkccfLKYmlCTeeblA3sGPDDPhdR2UpIDlSRZGjOP4V5zp7O3w"
    "tzafKXNP4LyakplR8brjDC7I6RsW0YaCrVXYBvnaA49RazFcFChWjb7xycM+e+Ynfq60yszE/NCRclgaC+AjFBc3F3vM9HBKmmDW"
    "Uxec+Z3RBXpDkZk4fW8+WzOrf6o1+ZlzP3bdPVfN7pudJimRBoG2jhQoUAgDamiTFkQUcjS4MQ67nbq0WfagZg3qm+Ys5LV7fvJ7"
    "PXdZGgQTBmTrG0zzl+okpRQAYhHHUUUb4Y7dYMUli93dmSgYA5UV6MTZRF45zHcpelBiQDUNiEX8YZA/PwvADQ9CsklWEuXqsTdO"
    "OvIL3/vYT8weQNvHz0h7npE/NkOERZftH8CgDsKhGc3l7MYLvff8U5a3xWUvTmtlcOKPnXnsHY/cMmdobqJS2zISQRRFadppdqZV"
    "qghASKzV6vW4obUyBuSuThBYplywiDIW1j8iJwfhWytCiype8pwMwRfCcVQQAZFwYnICEddeuPacOfMQ4NXXV658bWW9Vq/Va4pS"
    "IADtbjSCggfWDIBgQfoI2WPJGc9djTasW7XDMol16D51kDxSdhlV31trtlbk1SIBmhHB5OpPH3by6cf/yLSg5Qh7tzOp68nF8bGu"
    "ikJCj4CUB2SUNbt2HXvjR9vf4wwGpTOx+TK4c3jlcd8/8pFnHxjsm52mCYi89tAghBifHpszMG+XrXbbZO3NKlF9ZGr1/c/c/fjS"
    "R2rVhhRRhlx1DQ2eMcMP3BwB3kqmZ6CAC2kCYfwgcyok4zKiiaanJ/fedv8PH3D8DpvvPNgzRACrRl+/+YF///ran708/GKj0kgT"
    "hQGPyYvA4ZguzdStdsEoguYVwxldsEZcD+mzdwvWl+ATpUPlX0kzYIpMKeTq0Tc+f9RXvv2Rs4wTI3ZNCAiYEoUCCaDrPIHezCma"
    "PFE8o/IG/riOrpVzSkPGaA5LZlq7rjwws/qXv/7iB08//NkVTw/2DSZJkn2RQAAQiFPNyWPe9eHPH/nVRQs2sH+xk7SvvO2S03//"
    "9WanWYljTZ4mNzvM/YEqMwn25+szcDQscuF7EBVLCWSfaGjHlK8iIYTWaar090740Yf3+1j2xYrMnQYArw2v/OD3D31m2VP1akOD"
    "RuYuRF2Jzx7cwxPoqIwVEv4mcPOEkKWh83wDxFLeAbLRgrXCdHC1+XPtJkXsRmO0QwAUIIYnVp/28bM/c9iXUpUa1uNM8az+dLUk"
    "CjakObCtKwCo/G4vsEG7hlGHvy+39SzbqXYxmdX/8hsvfuC0Q55f/uxAz2BCKQfwpJSTrYlTjvn6l47+phkMZ0ICIikjRLzpwX9/"
    "/JxjIhmRLqNo5LiRrc8ZL5sD96V5gDO2mDNf0GVEUYEiVWmqk5+ffMGhexyttCJNwsRSACBiqpI4qjzwzL2Hf2W/uBpjjECQB/55"
    "GuPyLh1DSmCo/mVme0QkkPHrSgxF+eSaYyHAGe6sjiRPoGPPh0C5xfBZTpDMI8xgqj3xvY//6GMHfiZVqRACoQQ1cQNBPt6x8Fsp"
    "olEeM1+CicvvfOc73h9hSLWkIC+sa46O/QIM7kzMWHcqktHyVS994DsHL1u5tL9nwNT99hsJlBNT4+94277nfubXxtFECimEECCE"
    "EETUSTubrLPZ6yMr73n8zkatoXOjAEQeqlVmSJyhpVl8i0kWfQsc2Pz68NxRvILII3+yekRIkao0TdPzTvn9oXsclaQdgZH5gO0v"
    "gUIptc689e596q4ly5+pVevkZSaU5+hYUzL7lrmHnBuf5S5aDExD62IeTsAw/+DR+m4xXNaWeugN0nNtQRbaxG9jtLiz+9ZQYlwk"
    "REVWrr/7mkhGu2+9l+mJBdqfmvk+EQOXkL/DAPfgF4hH2svfHIhgDQOZu6G4FJDzlpGICKmIMLj/BEs+b3fst1NaSRmtXLPig6cf"
    "9sLKpX09A4lKQfi7DgERPv6+T5tvIIVg4tescCSiA3Y6NI4rue8Xd6h33wSy/HfAEPjGt+6a4XWQtj0vWvvYEy/P3hIo0lQplf7q"
    "y384ZPcj0zSJZIwYomuUW0xvvGhTLVJB6K4SfzacQ765UJfsl1rXsID3E5imQBA37WWOAbIm2ldLZi8V3c2QnbuM62snx+Y/DhcG"
    "9Bc8liGBANjfM/T9P33rp5efGclI58GxWHgZnr6aSo9eF/FoYecCSTz8ewLstqBguuEHsxFy6pHz3M7XYOgLzkFxraSUaybWHPf9"
    "9z/38tP9Pf1JmmQ73RYMQiSqvfb8dbdfvIv5R5sb6vJOUCDi/Fnze2o9SmnrPEDBEIyYTJ2c3HsGUCEcKjHyD2E4YrU5j0CF4wBA"
    "CKFSpdL0V1/84/t2PTxJk0hG3qMrFDSaFIhs5pe359xmyOQQhNIMZIueio5XXhNLyFmhoSaGg/Xu40Mb5pd9ypSzAil7lehOW68Q"
    "y54/BYQhHpHJ9wABzeqb/YOLvv3zK86OpDRuzfb490sKtrLRX7IARFReuxLv8MLnLwJ2qGUFl9CcuEoD2eokKv1os3dIWkg50Rz/"
    "2JlHP/H8I4N9Q6lKUbCRRf4olVL9jYHeWq87mrCEPqRUqlKNBT4wg+k5oE0ew5nAnruuxmHVr+N+kTtu3RHJ4wAKT91USpp0O239"
    "7OTfHrjrYWmaRDKi0rrRMZP10y8+FcuKu2ORfEcG89XkcdQQqOwDytt3gpBaylICqPAnUJY74EdT5hvDnm0BHlVePvonRFCp24Wr"
    "NdBAz9D3/vCNn152polr6Pp9w3PKVYlZ/CXaco/CRV/27YT3kkssQrp2w25Sg/72Inb2EwFAJ+2c9JPj737s9sHeoUQlmf+Nl6Kd"
    "nVNpR5VINIlPdaCTJpqU65Y4j58C5U4wnbEMIU7byalf2V3vAAtbiViKBDO84sJf+x0QCCanJ875zHmH7nlUkibGOdA7mWxQiACl"
    "lZTyzkdvu+/RexpRj9aaAe7IKJ5gjHkR/WOav0hfp06AQZxpFzpYgKhieKRQyTEIfkNLYC8t7BatmJfxltzlT8bzWLfB3lln/PGb"
    "v7/ul5GMlE67CHyyrsz0+sJXhZA/eC7924UNQKwK57dKMNnE3Mq4ZFdQcEHlcl1NREKIb/z2Czfcdc2s/jlJmrC9gT5dk4TAVtJs"
    "Ja2Sm4dhaZW4GlWjrMnIl7O9482jQXv8Z5e3K1+5/T6585DsP/El41ENLepPADqDSoEluAgphseHv/6h731g34/mvF/XHbN8O8yU"
    "DzKa7kyfcfE3UOTNp3kr9pD156q5Figvsskkg6B/JgblKAYlMJVUBrk8p7h4Gc6rCd0aMz+bgpqCHWWUMeKs4MYcLZRnFHOJr43K"
    "IaDB3qFvXvDFf9z990jGWqVQctHynj7/1LiREmFJHj1lHADPdDrrATyYjxwIgODN873G2CPmuEhzt3+QCJTWUsqf/f2ci/514ZzB"
    "OalK+EM3bxrdoQVSRlOdyYmpsS5vOyuVeuq91UrdJ27YvFOjcvNUQBYywXzRsPBgbtpqv9qeuZCDNgJ9viGDnLLmT0ZyeGL4M0ee"
    "8rkjv2xWvwg4J3alEiilpIxWj6068ZwPPPHio719vZrsPDJbKqzGN+goQyczP/q8+MdC9Q9E5fRf9GsxdGUWojvgsORQzBmJyPt+"
    "KlT1nnbM/yAdLoHZzjUXAgHj2QlRr/ac/LMT7nnitiiKrd9Z4GkOxiTTNSg8lwE4e9FBI75Dmj1iRDhSJPSwHQpEqGXGNFgwSUPQ"
    "pCIZXXn7387+82lDfUOpVmZuQmSvqIxlZZtdIcT49Pira14B7rAd/iTob/Q3ZENp5Wo9FwiRX35Z6KfAHBhCRy0SkYyEiKSQEUax"
    "iCMZSyEjGUmUUgjTapsjX2utlNKopYwyVArzjoB9+JGMh8fWHPqOo0776NlKZ7ZQ9jrytwEoUlEUPfXiE4d8dd///u+G3mpfqlTW"
    "1ZnmyvMVz2zjkCzrFcFl8zlxhP03BiXK3zghBhpmChEP/m/sTep2QuBiRV1KopzkCogOV7GYBCsWKL8YwMKo7rYj0pGM0jQ98ZwP"
    "Pbf8qUhGTrDlOUswmgIrfBxvkcjbnVYgG1pvUVSK5SNzuykxEeQMOXN759tdCCRCpdNIRg8uue8r53+2XqlrM7bFwJKD2dETEECE"
    "stmcemHl0p222I1IZ/0JhU3eQM/g7P55K1etqMRVMKRRDYQQySjj4uajnyzJPf8/02yR1ibvmkCbiaVdowIQTXQwSilkFEXVar0a"
    "11LVHh5fE8moUW2kWoF9v4QAEEXR+OTY2zfd8dzPnE8eg7VEeGVGgfc+c/cJ5xy7ZmTNUP9QkqQZbkuOl22pnFjCXGaB5UzSSwFb"
    "IePh5XCGx/oveGSQr/7yrSIoFJGwcF/G+/BipbqZraDPJs3fAjfmSrWqVRprRld94pwPXnXmTQM9gybR2U4EmPEUK1YZ5R7DJjwc"
    "YfEjNSqdulDpPnfanPCncgK22cSvj7z6mZ98tNNu16v1VKm818k/1GyqjdxfyhwWz7z8RPfOK1tAm6636cMv3Cckag1aG8xUj02O"
    "EpAQ0vyKZBTJqBbVY1mJZBTHlXq1Xq/0VCvVRrXWW++r1xr1SqOn3tOo9zaqPY1Ko16tVyu1WrVerzbqlUatWm3UehrVnunW1KNL"
    "H7zgul8+vuyhWtxQSuWdBwkpW0lz4dy1fvPli/rq/Zm6qMuIzbz4Gx+44VM/+kgnbfU1+lKV2GsMPUCLCVC6aDAculugn5aECBAQ"
    "IhEJP1S1ANIg4/n5kU7sAiEI/MasqtPDm7m3Oy8QyA3pOH7rJvqpSvt7Bp5Y9tgpv/jkBV/5q+VBUDBrLZDPHBGTP7FSYkWu04oc"
    "gyhwBUUM+zBw7DGmNA3dEQEoUcnJv/zEiyufH+gZSlTCYFNXCXt+yZh1TbGsPLXsMQAQQs4wqd12k7dfdvvFphKQQiaqDYSH7XXM"
    "ntvs1d8zWI2rjVpvo9qo1xr1aqMW1+MojuO4GlcjGUsh4P/91wYLN9pn+/2PPeOgx59/pCZqJjkeCUlppdRPP/fb9RdsNIMfIACk"
    "Ko2j+Pp7r/3UOR8CgmqlmqaJxbPRqdmRe3Mit/IpetMU0wB83MKhWMiSwzirzgr0fXYb+wmAzIG9wJnhBRC6wbm9oawEicjj1dnT"
    "1hHnvCldotPZA3OvueOKjdba5BsfPsM8W/R0GB5R314/zlcDZ9Z9ZL+iopbIZZ5CkXeK/iAc/aeWBdSd8cdv3fS/G2YPzOmkievj"
    "DJRitxBLmcjYOUTVuPr8K0vWjK2aPTCXyiYLCAIA3rbR9vWorpJUkEhUe/asOT/7zIV7bLP3W5numl/aIEVUjHnKjCPdGQsIQJ2k"
    "098Y+MSBn/v4WR+o9zWANAJIKYbH13z7+LPese0+2cALSg4Fb/X/8MORjKSUSmuX+Iw+T4wBrJxQEKTGukLJgRtOcWzbwEAhGxS3"
    "TImWn86EvhGprbGDGLTAqpHXUll4tFuvDHxCLJgvhAKjTIGdqnR2/5yfX/HDrTfc9uA93p+dL+yTQTcIJ+BerVmf7QOzwf2ZL3oR"
    "gp62QreXh+cQUYIMW3jQCNv/cfeVv7rqJ0N9szKqTz7yMbu22M7bwYomXYkqrw6/+ujzD5mCoQgHmZe9eNFWaw+t12q1BYpOp3PG"
    "x8/dY5u9O0knValSyvy30kprrXT+S+UDxoxVIYSQUkZZuSSlFJEpnDICku2GEREximIiWnfuorqsK60AIYqjsebowe844rNHfFGp"
    "krAgC5kppeIovvGhGz79w48IEFJGWYSHn5/tI8wWi4ASmJFy22I2csm7SEuOzUGdDBzWDA/xmZ7ufzH/dIgxRjOUwvJThTtgyQFb"
    "rpgl81HnDX2meXW8InIxkZwdh05wmYNcQERUrza+ev5nnnn5yUi6AZkjl7qBcJgbw0pKhtcWvkb4eFeOvRAfN7PYgqLCMt9KmnQk"
    "5dIVz516/ucb1Ya2gBBjlucREByYykmtORCUqM4dj91ShjRlv1FaNWo9O26xa7vd6lBn0w023+dt+2mt4yiWQgoU2QpGYfBL9wu6"
    "BAN0G5Tk+LF9C6NTI4noCImIop221l+44Zmf+CmYAKEyaaWR/kRRdMfjt5545v8JhDiKlVJkJ2foIJ3c84r48M0pZIA86XZ+QHnZ"
    "RhbvINbU5QAq+Vl66M5jrj/mu4QyFlD2Iu2LohB68peku3koZx3mNCGL0dnzlCxrKPDEM0+PdDWujE9PfOG8E6dakwioWSuPBm4H"
    "9J4FFBtZ3pJrn/IDwg07vMstl4SX2pmwv5IJyoiAoJ22v3z+SSPjw5WoQloT5IwC97x9VCkvFu2+0qCrce3Ox27tJO1IRFQ2tTSv"
    "6t07HRjFsYJ09tCsSlxBL/G2rFLF7pNQ8sqjjItCRKAVKSKdKTYQL7vtLxoItQCNSbvzgxN+Mn9oYQ7Ieqe+2fhKaSnlE8se/eQP"
    "P6RIR1GsSHlxeEAcp8mmSyY/GtERm4gTfpCzlakE16eASGaK8BIMG8FH2X3XJGfGm7Fos7VPdiYS6u14kCb6rFK0lwwT+7h02nys"
    "mc9w3DtLlept9N//zP++98evGbTDYzqhuymzDw0L4UaOVEacpZYRUhEZTxyYzNV3huQTJ27dlJ1zRFLKn1/xw7sevXWgtz9VaX5f"
    "oEf6ISY1RQuAgS1ZiXS90nj2paceWfogIOi8UOaLWqAgoj233XuzjTYFoV9+7cXRyZGMIpHLO02Vb/LCTPmjtFKklHb/0Vpr0IaM"
    "ZU5x80vKvAQywwEpIxkJgT+98pyrb72iN+oHjSOjwx96zwn7bn9AkiY56k8elyazN5WvrF7+sTOPHZ8crVVqmQQ242PzGRLmIzr0"
    "dBwuoRURyyhsgSqIqev5YChHLLFUyOqqGHQIqB1oWqCPyobHXicL/sVTcEbLDUzz2QL5ZGzi7FPvME+VmtU75483/PaaOy+PZKS1"
    "yss15C6T9mpx+E1Ya2MB0tTyO9/+jj/3xeLIyw1T0LmG2O+itY6kvOfJO778i5N66325XpFhxDkfH7w5ODoTZcbRlUKOT43O6pu9"
    "93bvtqoIK4EwR5EiVavU10y8cdcTt06OTXTa7XfusH9kpliObY8iX9Ndf6GwfyVVaSdtT7WnxqfGhseHV4+98drIyhWrVixd/ux/"
    "77vhu7//xmU3XdxTawghWsn0hutsev6X/1iJqqJ4grDrdLo9/ZEz3//0i4/3NfpTleaIDwvXCESQuQGzIyMietWJg/m8lDxg4XrI"
    "2Mo5DyRXAHWxDEZLLnDDZ5uzR2hXlq954EUsy5KlYDDtwZY5h91baYHpo0uzd2xIBBQo73zspvfucsis/jlEGlEQku0wgU0TeN+N"
    "fGcU2UXWHJfZ0Dkle8F8xRsLMxdrmGxNHnLqvs8vf7pe71FKEXN2sqE4QSGT15gAXn4vCMR20lp//vo3/PjunnqvtUYE3/hOCLHs"
    "1ef3/8LunVan1W4esffRxx/86XUXLJJCtpN2knTaSavVabY6rXbSTtNkujM1OT0xMT3RTtqtTnOiOTE5PT7dnu4knWZnemp6cqo1"
    "3Unaieq003Ynbac6VZCmHdVutlutVjWu9PX3aa0Ms/UPp16+/44HKpWi3QD8w9akSEUy+uxPP/63m/80pz9HwxwGwK09rfdcARnx"
    "DTiZzS8hFEUALr/FCwcpSPR4Zh8Ejn7hXMA6pvlVQs6x9VNRTEkRNKMcLIKiY49LZWJxPtmkgOfTEkUiGp0c3nfH/S/+5tWUO5Qx"
    "TzcIs/K4sRt0dWJ0Malh1HCwJYo8VjSKylTK6JsXfOk3V/9szsCcRCV2mFmEnBi3PkNyC1GgCEBSyKnmxG+/8tcDdz00SVMpRDHg"
    "W2slpfzcj0/4yw1/mDt77urR1T09jXnzFyCIVrOZJmmqEkWp0kprleZoEGlNloMIIos7R5Px6bgSeT2EpAEBpZQatCYtQY6Mjxy2"
    "99G/+dJFgcUfl6ubgdf5V/3sm789Zc7Q3OzsR+fDZgfgzJgUkTmd2z8E4I6UJcEIzDSMuzZ7RHKmqc8dBbP5txsfIAa0Hdfzcsvs"
    "XFrHVMF88RaWScE3gduGIRTAKGY54U0vzD6Jo2jN+JofnPiTEw76rFKpMFOXYAe/mdg1cEF2NwCwM2amjcM+BePtc/tjNx/z7YN7"
    "aw1npcYmkujPIshNKJi7rp/2GolofHpsn+3fffG3rtHaAU2BmEgK+czLT7zn83sJiXE10qgSlZA2ZRoiYq4xRfTG+Ej+cWfDRy3T"
    "UhdzkIGQRJqktWrthnPvWH/BhlppUTZT06SlkA88+78jv3FAJKQdkYIdE2V/5HLB0LG07UEKdplaIAeDg8nOAjxzVC8l2Bsre1si"
    "yDrJHMs9nIiZRBCjtKIlgvpzMmABEeXCcm4yEsaX+AbH4Vzb+lKiIi0Rrznz5i022EZrhSi8qhOgxKKzLNXO3h6CcSiR8X0ta49r"
    "S0JRHyJOt6a++/uvSQRA1OQaswK9mpypATjzhFB9B2RO0N5a752P3PLwc/dLgRTY91rJuVaLF211zH4fGmuNokQiiqASYSxlJFAK"
    "EMZMS2ultGl7U6VVqtI0TfJWWCutNGlNSpNpmLXWBEaZx/FjDVKKqWTik0d8boMFG2mlSle/OdQnmuNf+c3nFKVCyGyFoMX18yXv"
    "LiJPZUZZ02oB0QK3ErlDOkeUiIplZojpWsDZ45paKU0Ak+d/wsBBbo5GVnefKTK4iW6RJo1s1wMRF0ATV3QjBo5WFojXRJGIJpuT"
    "37rwS4lKKPBmcQyIgizcPjd0KJ3VBBekJ05mlz8whlbZ769JCyHOv+Znjy55qKfWa+ZWxGJuvOx4xMAvn7hJJzvRzJEshGx1Whf9"
    "67dhxJJvb0Okv3DMVxett35HtRAEgM6Ve/nI12mF2bg/QAnQzjw8XRhP5BWRbKvWZutvefwBn9Zas4PHe8zmmZx9yelPPP9wb61P"
    "k8rbWo5OeLPifIE5DWB+MxJDeXgsCIN7yA5S8nYrwywRLRsyT+r0E+/AQ2OoeL9ng7RMZIREjHVrfxjy8hZYb05UGAqSm9TlrLB8"
    "BRMGXQfYYQgFgyClVX/P4K2P3HThtb8wyXzefiWfoUyB7BVDiByNESGyoQRHmXMTBft3syEIGuQnevqlJ86/8tz+Rn+q03xIYpVT"
    "xBs6sngZFRT2GCIGBKR02tsYuO7uq5568TEppDL+gYwTapinWtP8WQtPOexrU+0pKaUmAg2gmSGHE49my1yAmZNlWKcUkZRRJKJI"
    "xJGMYhlJlJGIDJ+f8u8mpWh2mice9Jn+nn4zxy1OG8wzueWR//7hul8P9QylOrFjRja+4kGRmCcXs5MVGReDY9B8BE/g2R4B815g"
    "thLIs6S81e8GXuSLzHN9CVo1HOYMdseyxDCSuWjMR+hNo5wHV678yYetprxBLq1lK8IRxLnJgtKqvzHwk8vOXPrKs44vDe5+4hBw"
    "PhB3/sNODYsZBZhxT0LbC7/vdfahZG6A7/3xG1PTE5GMNI8ysRUj+nIK5Mctonf0cyw4y8SSMpqYmvj5ZWfbLxMFwFEIobQ6dt+P"
    "7LXVu0bHRiKMnX4NLUMAIxEhgVZKaZWknXan3WxPTzenppoTE9NjYxMjI2NrRkbXDI8PD0+sGR5fMzK5pp22M3YDEZKYmpjcYt1t"
    "jtzrWJ3T/cHXWWqtEXF8euw7F3xVgiB/WaDT2VAwawLwkyXzmavD9ti5S9y5H4mKtXb+YDNnr4xmCF18FMhqspziN1/ABgJ1kq28"
    "tmXSIPSVv+RXO77a3jN3I7DrgbgQISCDZOBL+C+IKlFlfHrsrL98x3M1yP9Xg3bujYXJf/DxRd2Y66wGQM6Yzvo8GV1112W3PPjv"
    "gd7BRCfsmGWjdW4ETYT+SJoVsjyZMTcSB9Ra9fcM3vC/a+978q6dttxdKWXcgRDDhCMh5Xc/+sP3PrZX0klkRdotRxoEiiRJxppj"
    "A30DjXpfFMW1uFqJq9WoWq3Ueuq9PdWeeqVRq9YbtUZfva+n3luNK6NTozf879pnlj/RU+tTqRIE01Otj73vJBNMH1A+rZOcFPJX"
    "V5379IuPzRmal6SdgmevZ8/I7ETLCJ72gZTYSLiUVfQ0SNmlzR0xfIsfzyiOGUWQJfUEzMeMGcc6cHBBDTlTwDOFw8IAPjvpGUhF"
    "XuYi8zZmYrMcACAmxGXZtorUQM/QdXdfdf3/rn7vLocqpQya5+bWefvqZUS4Yoyd0mFMqtc+E5Q5+gPAZGvifV/d+6WVz1fiGlkv"
    "c0b48UKqZ9Ak+2qOQKshUE41x3ffeq9LTv+n27uFcYbSOpLyousv/OIvPz1n1pwkTcyxJ0C0O+1ZA7NOOPhz79zuXf09A5WoUq3W"
    "K1ElklJgJGVXavTE1PgnfvzBWx+6qb+nf6o1ueHCTa//0e31ah0IEEXQuWvQUsjnVjzzvq/sRVoBCtLadTwseeCt0FUDDTFj3WIY"
    "OwCFk8bWRlSK37GtwiQ1BTGB76lIDkAlZCaJPmOUZzV5SYXIkahs/Qfuc+ExAOSFQkAhOQ8AhWi1pzdae+Przr6jUW04obDHM6Fy"
    "eQaTQJbc5lYgG7q65O9VCPG76371zLLH69WG1pozdgxBKZQSYGhnUQjlzS8P1o+azdlbH7jt0ZuuuOVi1/EU3oxAkabph9/78eMP"
    "OXFkenUcxagRCZM07evp//O3rzrp8C9sscHW68xbb96sBQM9A/VqPZaxEEJn6JAyHNJUpUmaJGnaSdt9Pf1fOuabEURA0Elb/7f/"
    "cT21Hq1UiUoz/4Nz//aD8clRmVWllunm2QeWX7NYuHSJ0bty+QhTC5Rc0U71TpZcTWV1P5GX8+VV7iWGz9lMinwjWEAg4+XqgRx+"
    "Aj04XjlhboFQ4BiDz3tlnSyhTzd28KvR+jVqvU+88PjvrjvPcYR8z4PyE8ffDyLg2nC/OY8Gh271v/zGS7+77ld9jQGlVXnamgcF"
    "IBG3eWKFH4IDB9FTS9hXo0nXK41zLvnuqpHXTAlUMthAY0OrTzvunJ033XN0dCRCKVBMdyaOfvcHt97wbZ2kY1Z5BnNmIKcGBCEy"
    "xoQUUholpJSRiLXW84bmD80eaCbTa89ZdNiex5AJDvN9yQjIAP/3PHnnP+65qr93QPFNkjGDydFmPf5aGWYHfuUfGm0Uw4wCejmy"
    "ANhiU8qiPh1AnRflUCp5sR5YzuwYnU0IQKGq94pqFwPCPbUDxyEMojRZq0pO6gU89NQ4a6S9tb4Lrv3FK6uWCyEz0I8Kixy7xPll"
    "OajAEtOJUQwpdL2w3OBfXf2TVaOvx3GFglwUZ1TCEyytdUWZ2RaFVSN52BUR6Gqlvvy1l3/w528Z26kuxygSUL3aOO/zv1traJ3J"
    "5mRciSCGTdfbXCmFAhGEIS3nEiUMwBZ+knWSjhDi6ZeeHJsaS7Fz4O6Hzh2clyWeFyYSZjD6iyt/kiYJgrBYIYELqfTNSyhUQhTX"
    "XU7rQdvOQ1mEfS4JJ0u1ZP5rwBdeTsLxxV/Ium7ehgC3S2PvqDxgxlLiyHok8X4YGeSJ5DIEA2E+cFtTzqLPBA281c7FAlStVF8f"
    "fe38q881FmOFGFCY0bwr1wMEoB6VBU7kva98YtmjV9z0l4HGQKpS4MQ2xs/lhQE6F1W2XSzU7W8wZihgFTmolBrsnXXJjRddddvf"
    "IhkFQhliozGl1HoL1v/jty8bHBps61aM8WPPPyClJEVKpwRKZ4TQlHFCzSBM5wQFFIC1am10cvjsP35PTVNN9rx/7//zib3AiYBC"
    "iPuevuv2h//b39Nv7JwIPQGtHW4geOwDOwkiPmilAAG1SAgv00MAGT0xfJDZ52apiFzry6FGTwmJ7ALCQE/v28FaExyvabTca3KG"
    "imh1vcQDlKhAZfUXrJOMomPgOcq5GQsMXXLTn55a9piUbiwQ2pjgDMZYrHLyYVRE/oxyMOwXV/54cnoCUQJ5IzhmPUtheHO21gn9"
    "MXdBZumGkn4yAxFQvdL45gVffvHV5y30S5Y2b9+MkGmabr3xdr//2mUxxDKJLvvnX6+85ZJKpRJHcQb5yyiOYiOZj2QkZSSFNBoa"
    "BEhVMjyx+qaH/n3s9w9+ZsUTKtVbr7vdlhtsY4dfBB6Dyfz431//myRJEIXFdnIPACp8COgnDGS9QjbZsMNhZMNbi49b7mPOkA5x"
    "UMpkYi5dgBixOd9glvBAbsSA7ObwXcW8QSIGqgxWrBMwVNeMzpwXUR58T8z5xWcEo9fPePmrhIFdG3oQUiSiianxn19xTuhPFfi1"
    "lfrs5mxQDA/jQq9gaD/3P3vvkd/cvyor2iNzEzG8L0uEcs7tjp3CWhrk/NJQYF0UOANGIhqbGt1tq90v+e71UkhewATtjuGi3f34"
    "7R89/eixsdFqtbbf7u/Zc4d9e+u9zWar2Wo2O81W0pxuTkw1p7OZQGvSUEc71B5vjryyagUS9FT6X3/jjXM+e97HD/5UmiZSyMCd"
    "wTyTZ5Y/ddBX9iatgcXK5+va1wkUYkswV/HaeGvPrMRZJ3DVO8u4RNdtsuifgFqFPGnaCY3JT+Eg4HaRFn102bycnlR0nGA5YDar"
    "nfGHQtkYsZgyZPJlP4WOPQduCMF9z/P4i07SvuoHN759s52UUlJID89yqnlrwcMuz4wNSgWntxCV00LIj/7gqBvuuWagbyiXvHCJ"
    "EoO+yqRHLmyJQn9u4J4BvkI8G+QDaqBYRGvGV3/qsJO/+/EfMe9BsmiXFRybPfDQs/edcOb/vbJqhayIVCaChGprIgKRf0ja3NgC"
    "MSeFSpSRrMgqKuy0k1q18d+f373uvPUM7wr9+bz5Kd//87d+dtlZs/rnpCrtEqpM4S1XalJQiFSzxjs5/xLJpnvkSRo5QuiYi9wt"
    "lv1w8rn+4C1u6+bDwkbYCwi7yTAzkQoz/cL1RF5eGITs6ACC5+JEi95jkP/iCkIp5OjUyPt2PfR3p15mJpKB6tC+OIIyPQB7RBR+"
    "HWN9PvjsvYd8dd96pU5IbDzmnxMAxQmv93FQYU3kLCykcA5AHhZGhpk8NjF89qd+cdx7T0zTROZGDEXMS6k0iuJlrz5/wvc/8NgL"
    "D8+bNy9NFClCkYEzoCmYsRtsyCwPCXJ8auxdOx148bf+rrQWGIIK5idONicO+OIeK1a9HEfZrcj5nIhhIgRgYeyBJf0vC35keXns"
    "IecByb5EmFcP6KbqmEM86G0+HmkTtrxYpClyLyh+XbhJWXa9sFhVDjkRzxYsz2Nk9kPBCrYVBpPp8lsEEEQrmbry+//dcfFuqVKZ"
    "+Q3b3+TFErNYFoCCmLgMjgSAC677VSdtoxAU2l1bN0KG8xbwJwrhOnZQeYAMEvNoJvSWgFaqr9H/7Qu/dOMD10dRnOb+qUXEV8oo"
    "VekGCze64sx/HfaOo1eNvqYgRQEpqVSnSqVKp6lOU5WmKtdJktagSWvQAAIp0vvssJ8RtnNpKeWQACLe+vCNS1c8W42rzNEboTC0"
    "4V0YUUGR7DWWduYUACm+SaUjHVsat7NbyJ2js9+Yeh9ZfHjRQ8G2gBnyFB6WPIHC88d2ZC/EkK3kG4BbNAS7pNpz/nZw3DDyNWct"
    "sp0rsN1JLrj2F5BlzYEz76XQzMGbA1BITwk/NQN0PLnssf/c98/+ngGtlY+kMg5gt2EvUuARHxqPUSg1QDfrD1yJSQgphTzpxx+5"
    "/5m7Y7YHir+kkEqpwb6h35z65+8e/2MgmGxOxlGMgAb1sR9ermZyZYRSaV+tf+fFu2UT6LLxAwDccO813OwY/BwV8KNruO8bURdR"
    "u+95b+dQwZlJmW0DEpcPo9OiFoifGLBzcg848sFZtHMcLqlhNbxJMLMSKSiksDCbO+ZAgYU62141loNEYWdh93R+nXnkXu8G0Vr1"
    "9/T/98EbHn/hYSGk1hr92BjEkm2d0aG5ARJ2EdFc9K8LJyfHpYio2M5hCdDqA3R+hpSFaNFd4D7y4wbR9oO1EmOldRxXWp32R884"
    "6pElD3TbA5mEQggz/PrEwZ+/7LTrt97gbWtGVpHSkYjMzhUCGT3THoLYTtprz1lvw7U2dhQGzv0lklKuHlt171P31Ks9WikGS3j+"
    "JVjML/Y+Yg67Ux5vwTFyuz7R+lM5NxMMilbLWQuojxDGYTB+f1EvScRJ+CHkb0BevvF4wcSNt7g7EfCmnByoSO5qJShyqF00BwfH"
    "rUyBE35AimiqOfWH688vNhfIzDKCfyF431K8gMwn/cqaFdf/75reRq/S2k42zHsgLMgvAkI4ebFFXLvBoWo0BxdjJIKL+iJwkWkE"
    "AEqpWqU+PjV23BlHlO4B5OcgZJr37Tbd8aqzbvz6h74rRTQ6PiKFjKLYQI+5WivbnkLKTtLecsO35camGFSmpuB58Jl7V65aEUex"
    "ZtEUxKd5GA5V/RF5XgmG2HyBnsJgGuQgXo7qOKKoX+ESZ5oxtil5Ccfk5vdsEROGYnA/TI/z69Frb9xixyDsKweqqIyiVDQAzp4O"
    "5aASMWkEeRaQaAyjeut9/7r3Hy+++kJGjigOwvxPhEJBTKF0N5/032+75PWRlZW4QqB9RQdhyEsBPupGdJi/q2R9M90Cl5xy2W7m"
    "KuQLkrMjJ9VpvdYzOjHywe8devcTtxnDqeD457iZKYdqldrJR3/t2h/ecvCeRzRbUxNTY8b7BBFNAoR5lgIRBLx90+09wBvDU/zu"
    "J25PVco9TohBXsA0XoW5p2fqTdqQJTJ/cWKcCcytcgggrLJsMc5D7EqCXZH88YUjjQHD6d39Y5+/n27HzFCAoNRwmZiUyAz+EEO5"
    "TbcqgwJpgp1ekDOD99kDVrHnhh2xrKweXXXZzX+2CGsXGxhHlhZhce4ZsYMUstmZvvq2S+txXYN2VCGHhPG1zCozZNiQU4tm2Qhc"
    "fUC+rAMQhU1LQccksI8yD1GBVKtatT4xPfHB7x525a2XmGwpolJbLLDlUKrSzdbd/IJT//qXb1+919ve1WpOjzfHUGAsY4HSxoA3"
    "6o3NF23l1TDkZuRCSE3qkeceiEVszgh3bREj3xBTp6O/9DnbH8k7Kp11IDAJI1obOXLkD8qtQXIvNkaMACeL4AcQgVO0sxA+N7VH"
    "9G14kUJGtYffB5+hr/D06w4qJEmXsI+QHSJZxQQBkYnYvMs5MxKCJlWv9l51x6VjkyPGPwrsM/QQIPdNhcfNY9mPxiXNZFM//dKT"
    "tWrDEKc97klpvKrr0dim1WHOcJgcEaS3MS8adgE7wyUSCECpVtWoCho+85Pjf37ZWcYa0Ron+d7x2TzEXAVKqXe8bd9LvnPtX79z"
    "zSF7HCkjuWZ8TavVEiiEQAXpUN+sDeZvFOBL3AXk1TUrn1++tBpX+HLNB94c7nDtnMdq94SRwHz4uLqLg5pEIZmI8xwC5BqZAoUD"
    "SjaYxXmomS93zhvoNEnFg4R9KBAkn2FZGUH+62QdHTdMoGKnTs6igMnocqqF1ZM53iUREmqiWqX2/CtL/nv/P3PmmAknt8kvFBSb"
    "ohs8ZN/tpTf92SDkYTxnMEhkeAGxPNkwwaY7D76olgHO2fZyMd3QVGkto6in3ve9P37j8z89frI9KWWklNIFoM2+BGOJpVSqtd5j"
    "23eef8pFV59x08lHf3XdOYvaSVtEMlWdeQMLZg/MIwIUgtEFKN/MsHTFc8MTw3FcYaxJRGGlsr5vJnnRp9nEg8APekFy5AQX08xs"
    "aMH6J7s6hBzh1q5OYpvF0Vh4W4asYDXfTDumBAFQYUUTB8xy1bU9DqEAH7rinPzsMH61UAHqdHPRbKTgpO4eoGDRWi8BzWx3KeTf"
    "b7uESAtkYbuYq438ukt41CPfuNMwIu95/Paeep8mRVzgE7ZrGCBj+dyXWQcXVmQJnZawjDrv5Xi6/NqchaRJE9CswbmX3Pznw7/+"
    "roeW3BdFERCVUkeZnFIiorkNFq+3xTc+/L3//vJ/3/jI6a1WK2mqhUPrVCsVrRU/METGHtcAsOy1pVomQgpbaCAn9lIuq8j/vvsk"
    "HHsMwbOhtLV3Pk3IxbnOmpmQA4Nkd1K+7q1Js1dgWmcD28+6WwGzhB8kj8FIHk0ZQ1iD5cY6X2sbXpbFVyJb8miHDK6OJfBkWkzE"
    "j24O6FMDA8gqr/9Yjqgm3VPrvfepu59Y9pgQImc+EbkRuAf7i9D9zaqcSAPAtXdcMTo5GomIghKmwMtFl2jrU36d7Ke0H2HHP3lp"
    "uIwOFSxha2FpVapGR6+Gemc/9cLjR33zgF9d9WNAMtUOhTHXHqtXCIEoUqXanXaj2lgwe+12u5N21PrzNrQYAMNT3NN76Y1lGYas"
    "CRSAhsyJSEjLrRCIAswfCTO+kCiM6kCakDIjQkAps7+IYS+GhRERF71TPkNk+shC+pctGIj8JCU+XkPWhlqNAF8anHnq7EdDRiP3"
    "tg4HunZUhyw5jBkx5CUSj/fiJQGiN3VApyfxR3AgZTTRnLjmzivAmEFjbiTKVS75943468z7OSSASEbtpPXfB66vV+tKa/AylHz0"
    "mLGfiyiqV+piiUkzAIt7KlIDfAdNdvsgm7KhqR9SlTbqvSpNv33BV2975KZvffTMrdbfFgDSNBVCInar9ECiQJBE9NxLTwNpEeF6"
    "C9aDLqNQQwtduXqFAKk7WiekUaNApTOIFhSRt1izZ6vz/+eJECasIMt1EjLzq0BhMre01i4Vt2TUFMa4oXG4yjsVR0FkcjIq4JDO"
    "WdDSqtwm8U3gEPmQIUgsyi12newSPLaW4y1AYF3u2hLkCACbltpX42uALBbANp0RUf33/uu+cNSpjWpPOFL2qWiR9SAAz7FDCyEf"
    "fPZ/z738VL1a51pEZOxXSz4kr9hjTDwsSRb2ylGei+NcbPxSDZ2IvIxqhgwxA6UVCjF7YM4dj956+Nff/ZEDTvzkISfP6p8DANbl"
    "vEhvyU0xcOUbKxAx7o3mzZofsv+Z8zkAjE6MdJqdjTff9Eef+5XSqSadqI7WlKo0E5xZxjGC0qqTdJRSrU6zlbQ6nfZ0a2qqOTUx"
    "PT42NTY2OTw2OTo6OTY+PTrZMo6lbWM/EUeVSlSJhEQhzP1unfaBsWLyUBUkz4LRY/qT5c35bJ9c+I0sYQ6dtpYbzFnyqfXbACwc"
    "ZCxjzK1qPjGwkn+eB2WTm/Ic28wXzE0h8/cY2EBk6ibitDnQtUr9+RXP3vvU3fu8fT+lUoEyRGzy9xkhQoEWkv2Mfz9wQydt12s9"
    "BJqpWqCcyQdlo09iKE4QXmL5gkEOpJNtEwZTukLym1Ucojcvh1SnvY0+lapz/3bWdXdd+clDTz5qnw/VKvVsG5gj1hVQgICgAQDG"
    "pkZFhDKO+nsHSyFr0xq1k9bI6DABLViwcIv1t4L/v38RUbvTmmpPjUyMvD6y8pVVK156bdnzryxZtnLpyjeWj02OtNKWlFG1WqtE"
    "FYSMwpgFBFG+qtAqvjDg5FoXMip4rOehiwXSaI55oS+xRCx+0p6NrRvxeVU75rOssAggp1LLRfMgCip4JCzLNwjsFGxZJSBR6T/v"
    "uXKft+9H4Ftc+g17VEbOpUhGE9PjNz/4r3qtoUnbhPLc/gspFPYDX7n8BuQcLnJ8vrIBE7edKQWOQr8C9GQzjBVPBFppRJw1MPeV"
    "Na98+Vef+fMNFx7/vk8fsudRjVrDMJmNmh6YeTUQTCUTsioFimpcC3cz6/2nW9MTUxNCYr1e16S1VhkUiCHjxaNNMsI4MemciTCu"
    "Veu1an12/5yN197E/swkTVaNvv78K8899vwjDy6577HnH175xoqk06lVavVaXcjIpBwQMVazJ7hzxw3kPAsib3DpPDWJ+bIxgAjC"
    "a5eCQQGwDPB83QswxsLEh4PO8xm5sYRPjUNfjEP+TLvIM0evXczqdE26Vm3c9fhtIxMjQ31DjiNNzDYagIgiCLwFCBTpSIr7nr57"
    "2StLe+t9qVK2HhQu2Qq9yBLyiOHc85RRD339QEl3WyjSqTzjMjibs4WfpVEhsRo41Wklqtbi+jPLnz75vE9c+M9fHL3Phw/a/YiF"
    "s9c2f9eMcpGEea0KEyKiFCIZe2JR7wgBpVMtFUqUEAkUxhERRaD1LlAmRTkGkE9d3EDALD2BIpbxWnPWWWvOOntuuw8ArBlb/djS"
    "h2976Ma7nrh1ycpnm9PNeqVRrdQoU3U6+g1CIDD0XGY5jm9deRnaSO4+DVn05PzkyEf7nfog34lesYIe4cb5DWoKKXzkzHyJASGI"
    "pbIp9xedEQ+SpmpUffG1Zfc+fcd7djpYk5YoeaC1LbkiBgznik1AALjpwRuSNKHAozTknWf9XVECYyPOLTyK3EubjeVzY2rwbc8C"
    "aScxNUaeAe1wRocMs8Y4eySaSFNaqzYatZ4lK5795m9P+cXlP373Tgce8c5jdly8SxxVzPLQpJQmrUh3tOZRliic2pqViyIShJB0"
    "Eg/KYBLysM0stx4LR2zCjCa96XMW3CRQzB6Y887t3/3O7d/dSdsPPXv/P++55j//++eLrz5fqcSNei+RYW47e+IgBBNYr0Dom1c5"
    "kITRc7HAJPWrH/J53Q6z8IJwA+jQZWUTZ/bzHovXxgQF6RpygnXmn8W8l/OnJ5RSNz/0r/fsdDAyBjt/HzpPiuf2iySlbLan73rs"
    "9mpcMyFF5IE9GLxx8OMzwVNkIBvh8qO/xDs+3OVFTRU/1lxWmz9LJkakIydwMuBLNarX+hrjzfGL/n3BFXdcvOUGW++91bvfse0+"
    "22yyfaPeAwCUoGpTJDXrT0rWrhkYg4Jmq8Xl1B7mwdb9W7HEgi6KDQPiCSkMxGTeYSWq7rLlHrtsuccXjjr1P/de97dbLr7/mXtA"
    "U2+9L++VXXAGkktfzEVA3ForLwp8bov16CvMqtihzO8Yrz/1uPr+nIsCBSRHdp0VlgsQLioTyGOcoUtRIC+CQNfi2r1P3jPdmmrU"
    "euzq9+zdrTUiePkI8PiyR1587YVqXDVPEtEX95TIhYIqkcjyhTxSCXkxbsXpWOGvcCk4EtMfs5BhTpEiv8kHDRYMJMMPJIiiaKh/"
    "NhA9vuSx+x6797zLf7zx+ptuv/nOe2y518TEhBSRQtVO2tDNygpAykiKCADanaYmjQDGmImZd5BHny7LlCtthRGwBK8gN7XIL10i"
    "rYn0rP7Zx7z7uKPf9eH/3Hf9b6762f+eukNGUU+9V6lUA89zB7YlODzkoZkcIypo90pUfhR0wA7woBKDdPQpRe6MzCdm5MUrkbOk"
    "swU3cUEQZ78iH0gjAIFGXavUl61c+viyR3fefDcTKFG8nKNS/Oeex+9stpv1aiNN0xnLcbsXsWAvl8OZZPse7pbJZ8N8Jxd8tXNx"
    "OZJHwnJyAs9lrvA6hGOVIDPz01oBQa1Sb9R7NannX13y5EuP//mfv+up9jR6G1PtyenWFPDTh13HQFCNqxWoAFGzPSlQiOhN0uf5"
    "pISyskZ7tSCBl8/tAYvkgznZykAhAITJvZFC7r/zgfvt9N5/3nPVeX//0UPP3N/X6I/jWKmUCt4zLJ6Bsi1nC9W8+iYqa7bIu+2J"
    "u7DZ4UHgC+r5IUAYvYVETvXhUuY5HYz5liOySs2xs91sz2ofc/G+FM1287aH/7vz5rt1+2iCDUDmjLnniduZ+0jRwdVnQWHRfAuB"
    "qypyeM6F85Cf+OkCBvNcCOKTcPLNqUPKrw8bgBvHFD5DRD6tyzLjAaEa12tRAwgIlJCgQY1OjJiDVmRXM/EyqxrX6nEvpPjqqlf/"
    "98SdQ32zpJAyM1yJpZBxlF0RZjxs4jillBKlEBIFAoju94DOhNqIWJaPgyx0GhElSjvleN9uh79rhwMuuOYX51/909Gpkf7GQKoS"
    "bvJN5BXwga6JAMpHlUx2FojACDBMiSIAL6XdbwAC8SQW+Ouu3yO+05wJERKHJnzyTHbQ2jDxOKrc89SdWitjHRcGu5LtAfIbUgjx"
    "6poVT7/0RK1SI9Lcu5NJtv0RS7FMI/JSLICr9HwnaHLsaMiHKIjguwYiPyvBdXjkhSFzG3HwFVEFobX9msz4gLRjsRKmSTo8sbpI"
    "czdP0Fym8+cuiGrxyPSaD5xxcCwqggQCRlEURXGWrGqyxkzyGCAKEck4FlEcVxq1nkatp17rGeofHOyZNdg7OHdo/vzBhXMG584d"
    "mt/f6OeuvUopE4oohSi6FNuXKIU0s/Bapf7Z93/5PbsedNrvTv3Pff8c6BmQQipSrN8NjHGRK+apyNPHkA3tHqotYxhvOT/by/p9"
    "zhDmg1n0+m03xOZW127Q7KOIIUWAzchI1yq1Z19+asXqFevNW6QLcydEjIj5KJln+tjzj6wafaOv0ZcxIEreBgbZXkEzS4h8xIfg"
    "i1CZizpjj7igWtZ9MgiFlY7uMVgDgoBcTR6D3HMJJve30MmN3Y8gTaBo5eoVYVfqvWxYb8H6ipSUkghSpUCnBERp9mY0575pynXN"
    "2QvVWhG4GBtAFCCqlWpPo2/WwOx15qy3wfwNN11n8eJFW228zqZzh+YDyOzzJIUgeOgnsYEsIkqQWmul1SbrLL7421f98sqf/ORv"
    "30+SpF6tJyoFKmUiupIknBN5BywGI1ACBo9QKcvEDjS9dRoGGHkns8OIOD+AuFM/OQPmYK6b/7Ts6wkokvGasdWPLn1gvXmLDBga"
    "vMTID1MjAHhkyf1pmvDsnWx44KyaeAZwMGzEwGiPF2o8eZ7yhcjiMcnXVzPaF/PGsJE/ftHpQwG+x6ybRhD6HSr5Q8/sZUgZLX/j"
    "Rb/sDHOi1p+3IbVJ92QfH0qkIBjCQi8REqezsEWGBCCsUoXandbyV196YfmSm5N/a61rtfr8WfM3X3/L3bfZa8+t9918/a1MWWUY"
    "fubIBwAUGVvI7uFIRkopRPzMEV98+2Y7fv7cE15Ztbyvpz9VqWe8w7kISOg87bFMrsWmk+7kQEbkKk53HN7PT3ubUIDMNxUcJERc"
    "vQP8CEPmIxkmZ4VyTHvaKaXue+rug3Y7orSsizjuYhqAR55/OIpiW0s4m2tmI8O9kLwkwy6RtMRoRBYTAZ6Ild8rDE/wZgLuy4hI"
    "sOxFtP2BRQJY5COHy5wsjWFswEVs2SSkEsUrV61I0zSKIq2155GQP4RN1tmsVqvbXlY7Jz+OxhIA6oySmAczAAku4XWkZEDESlyp"
    "RjWoAxJppOHJNTc//J//PnTDQO/gZmttvvfW+x2456GLF22ZXQhaG8fsvJc1RyNmHyVBmqa7bfWOK3/wn4+feezDS+6fNTA7SZNQ"
    "GpuFJeZQjG/o4DHRXH9lk5dcVYFU1jtY2zpgzLyM5uL/DJ4J4E5IHx52XFPvlOT6AWKfvTl4KlHlkSUPeMkm7OYQ1t+fiASKVaOv"
    "L1n+dK1SA+1yTqmoDkIsm+CG1ZId2KEzunCjSE5kDeSR1pMbKHDjBOuQRla95DJZIXTZtsUoOTk3t+SmPMvbHnOadBTFK1aveH30"
    "teI0166wjdbedP6CeSkkCM4IlIfyAQtJcLmAZGxFCJBAIEp3/FOWTqk1KK2VIm1u8L76wFBjNqX08HMP/uDP3z7o1L1P+Mmxtzz0"
    "byIdyYiIlFbWmybgbkgpU5WuO3/R37533d5vf9fI5EgUVSyCxvIkyL52dFx+DlKEGCixeQtmMgnAUNOLzvLRcrHJ1bzFxcLtfIUt"
    "G8inR1OIgKA1EDa2pCwvUmtdq9ReWLlk5eoVNhbeH+mg1+o9/dKTrw+/FsuYq64QQyY1+rFqpTx/Qj+bxscTQrDBUiaCUDjfRI9L"
    "NKxQ0BOW5MRHszGYhUcYNReo9fi8JpLRmrHVS1c8W1oym+c4d3Dexuts2lZtI8kh7QRAaJR2WVKPAS0QiatMmMwJMyVAnnMgBQph"
    "cgqkNMk/gCAx6qv3L5izUKD4973/+OgP3n/Il/a95vbLpJCRjLRS3WIHTTk01Dfrj9+8Yt/t9xsZWxPJONesEdOr+fyTgi0Lk634"
    "1UZOMULP7aLAXAzs4pDKFkahxwCbd8RqbgeGM5/tvEwn39CFAOIoHp4ceeqlx0ubfOcObVJ9nn7piXannd1Pgh2njhUjgu1bYDIF"
    "i8acjuilGroCK3uXxNLSnfSJ+zaF+RIus42l79kKlUnR7awHA2cJ9mr8fGdEkaTtZ19+0lZ7hIywa4bqgNsu2qHd7KB1isyXNpmF"
    "LDPhS+ZELaJIGodqKYUUKBEEaCRNJqMmVWlHJa12a7o9PdmcGp8aH58cG50YHRkdXjOyes346jWTq94YfX14eKQzlWqtH1764Cd/"
    "/OEPfv+gR5bcH0WRMKFoQYWaZc4KpVSj2vObL1+80+LdxiZH4ijmJwwW4t1L7Bu83BMsCSqwf18XK3N0E4xsOI2AxXugeFK58LKM"
    "oOQASye5d/ZGtl739J8iSZNHlzwEZaPNiClvBQA88cKjQgjiFZfVwnItDYamZmWDfefsBcUKMXiLvggsGPUGIzK+dslzpaEg5Jyh"
    "DuRIwr5cjbhrZk4uQRSPLHkQmPwyiCgGgN233usXf/8xARmveHMDG9i/2W4maWJCaLLyhmVeSCHt9oijOI7iSlSNK5VKVKnElVql"
    "XqnUKrJar9Qa9UYtrtUq9Z5aT6PeU42qtbjWbE/fv+Su2x69KZL1Wx6+8c6Hbv/4QZ/98ge+Wa3UeKXLsS8Tp9nXGLjw65cc8fX9"
    "lq9+uV6tGzMl4jASeH512tLAcskMghPU+zqYgtrYcdgcDIKWfGs9fljmGPi+qwA+f5O88ACXdeNwTGI7JhdnExJRJORTLz5uAINw"
    "A2BOZhMoU5UsefmpWMbZYJ+DYg42QHImdQWml6O7cb42M4xx7tZW8297Vm6bSYX5tHGYJyrcmgHz0CeecC4vsaAVdk6I/NCyynNN"
    "1bj26NIHJ5uTvfVeS6a1EgXz39sv3mnjdTd98dUXKpWqQT4FiE7S1kSbrbvlogUb1Gv1erVRjWvVuN6o1XtqPfVqo1ap99Z7exu9"
    "jWpPo9pTq9QqlUo1rlWjqpQyiuJYxiXtWvjrlMtvvfhLPzupUenTWp976Zn3Pnv7Dz/1y8XrbZ2kaRTJXCDuJscSZarSBbMXXvC1"
    "v77/W+9ptppSSG1uX4IgiBv54IvYg6dirk0OvunsuifOJiK+QsiJmMkXxiOzQ8YS6Qrj+BQoQmwfkp/wlRfruhJXX3h1yWRzorfe"
    "l7cB2ZqNcqtXkgJfG371ldUr4ijSjlAYDmCLqbO+GzUvv4hNtyBsngp/FkYZeG7CCJ6aA9GLpXX1lW+uzh2qMJyl5kaXftRW1rVV"
    "K7Xlb7z89EuP77h4V00kwzBvVEr11Hp33WLPp55/olappQoEyGYyvcE6G532kbN332rvWrX6/1UdAy7mPh8noAA7WjAJ3u/f+4MP"
    "PXX/Bdf8YvasOXNnzXto6QOHfX2/cz/z2/fscpBSCpn2zaLAUsgkTbZYf+szT/zZiWd/sL93AEF5DBoH+lohsZ3Qd/M/JT5apIK8"
    "I/eRcB4tDM4kCM3skKiU8U/e6kcmtylxlAsfZiQrrw+vfG3NKxuvs5h448sikggAXnz1hdHJkUjGLseWirEcJWEaId3AYl/eI/N0"
    "2sQ6/9IBCviWUlwfERqykDVICmbRhe4E3QJjOAV4lFgAAJJCNDvT9zx5B+9y2CeWfa4H7Py+OIrMCZdSMjQw6/dfvXTf7fePo0qa"
    "pjZ80v5KzP8o9ksHv0w7TUYlI4SQWTtsAvyiSEZCRKYPeddO74mqMSGlOumr9k9Ptz72/WP+euMfpJRaK+e1xlxfIxmnKj14zyOP"
    "e++JI+PDUsbWL4fYnJGL7gCCOBMOswVzGzYZBuCfCk+/pAAwx8AFnlfFGKR5s7xhS0yhgB/tzUMxQ8MmpiaWvfoC5LRFV/jzHubF"
    "V59vJ22Rj9w99IQ88kco9CXfZNrbKRggAwTAmUI8k82hTAI9lm1xLskLM/aWct9lZHcBFbR8tuglCoKuLCioKY6rN97/L6XMuJcY"
    "fAUAJKQkol22ese2W2zX0tNxHE0lE+/b45BN1tmsk3ZQgImxMQtX5qr3yPwXSmF/5dY3aIP7Qt5oEEDqIJtWu0WpNi88SdJKHFer"
    "tS+e96k/XP+rKIp0nuGJzKwJgaQQWutvHPfdLTbYcro1KTO7gNyYkdG5EbBEiYhchl1GmLMSYvRZLWx2zAUeQSIWhoGJyHcIClHC"
    "X8rZY+TinjhRnxAwVakBgoLiWvBX/8KrS8MsKvAthHOiGlFRzOXKHYJguBu4XLjBue96iSzCltU0rCL1CiHHmrNsQshdeog1TUUJ"
    "D/JzKTD6Nr/XoBvVxqNLH3py2WMmgdAfOwsEUlpVK9XD9jy61WqhEFrpBUNrKa0yyEugJyhD786icvv9wj/k5RBjUpDSyqS73nD/"
    "P0AgEZJCIFBaCSl76/2n/urzl918cRRFqUpDjmL+CPt7Br99/FlKKZc4Bjaf2HHRPG2M/0aKHi7k3fjIkRqeyoXEnC49TkR+FHtu"
    "IPZ7ODyQ+4khXycuWpGHoWUYz/OvPFfUZgjzTM0M+IWVz0sZhcub2BblThsFfjsKtBNrxlMouq57KHLpNLsIKAceNGTdzhienLtI"
    "OXiICtmi5I8DMHA7Zq9FopxqTtz4wPUZ/x5D5Mikph60yxHzeha2W20BctmrL5hAbzMX4HLHzCkCSGutlNJ5TGXxPzbMmIiM9Ebk"
    "acaGXBfJSIro9//59XX3X9E/0K+0Ml08IhrTrp5a31d+edL/nrojkpFmoZr23BBCKJXuu/173r/3B0YnRiKU4HK40OuI/f4U/DK/"
    "wF3PZ8TEeYdceUxu2Mz57KyA5r1GzqFy0j+HwwWbjvEVsm1HTIpEFMto5epXiEgID8MV2dRNCKXT19a8EglpmZHkJdn7ehPmfYvB"
    "BMWz/IIZZCWBlT75CrMAgnZhh7m9gVdaBpgBEyWTdza58QPHI2w9xOMOtNa1av1f917b6jQjGRXdnhBRqXTh7LWPfNexE8no4MDg"
    "DXdd+8hzD9SqVSmlLWsE/4VCShmZX7L8P5lFkBDG47uTtKebU2OTo2+MvLbijZeWrVx688P//vjZH/jm+V+MqeJ5gOfRfZGMFdFn"
    "f3b8K6uWRzIKHPLyukgQ0ReP/eas/jmdpOOppYBY+FJJDDewLo4d3uj7+djZN/saT+2AAVrn34AcGPRTw4lnS1hEzx9g+5kcGiiK"
    "4hWrlo9OjmIe9py1RiY4XQgxPL5mv1N2XT36htGDAwQ8PhZaV7iy0bNAQD9I1MnnSr2qLeXJ0+uX8dJLTCHsEi86W1giFStYg2mb"
    "R6TjsVs5a08IMTk19sdvXrnfTu9N0xxi98MzhRArVr18wKm7t5rt9mR79uDsrx//3X23P6C/MaBIKZUmadLqtDpJu520O0knUZ1O"
    "0p5uTTc70+1Ou9VpTkyPT0yPNzvNVqc13Z5sd9qtpN1qN1vtVqs91Wy3OkknSZNO2jLmQhOdsaST9NX6rYo144NY9aaGKI7H2sP7"
    "bLffRV+7Cs1lVQx+1lpKecafvvXTS8+cMzTXhPw5cQXL+vbWJhaPAmISx1AaEg5wrMLX+ziJ0G8WiWV7+Wch88KBIIyHSfACeykE"
    "0oDi+nNu33TdxRmuTQAIkVWnrBp9fXxyTMrIBXRyWqrNgSoq4LkdrI+9oLGSRSybEvOQaPJ9+Cjbpg6BZqYCRL642uw/HcbWoDOh"
    "dPdHRi10HqfOvws9LZIFlVOlL735T/vt9F6GDbidKoTQpNedt+ioPT7007+dM3do/khz5PO/OGHd+YtmNea0W50k6aQ6aatOqpJE"
    "J0maKJWSJlPtUFbqaHKiaNcMo7HxNWwbFCiEeVm1qF6PGsq4dYjw4DCPJlXJYH3oX3dd9/vrfnXiIZ9L01RKGTx782w/duAnL7/l"
    "oonmhMSIOLnNDpPYoeeMR3ICpFfDuqFVyYvyAWu2PrmYEoOK2P0L832RuA89snzZooUHD9skgWKyPfX6yGubrruY86Ij+zxWj77R"
    "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"
    "ZAcEIquFhBUYZMw60zRHcRTFBlZCgdfcdcVJP/soCgJCpXS9UX9h9ZLXR15FRGXdQfz6cIsNtobc34KPO/MKBsleAm4xWMsQ3zK0"
    "aMKBAdMhn2Kxy8h6RZHLWeDDLsvnYlbCuR8b7+IsN57Y2YqIqUpbnSYwx+woH/doZaYk+SryowDRhVDl5WAZ04MHDjKnWs7lZOih"
    "tuSc/H4gLy+eDdJdpoMTw7GgGiyatXufEFqmtNfCI+edMytqe9zlEYmYZ/Ci0qq/MXjTA//61VU//tyRX80g0VDCikKINE133mr3"
    "H332V5/54fGxjL/228/+8pSLFi/6/xIolqo0STod1UnTNFVJK2k2W8120ukk7en25FR7utVptTqtZmtqqjnZbrfHmqOPLXvowSfv"
    "q8RxHFVUokBgFEUjk2ueXfH0wtlrE0uI5r8WL9qiUqnmWXr8JiQ7ri2YRnAaFzICJTF0xBc0ElsO4CnkKVB1IY+3LPDR2d8qSj8o"
    "207I+H1ImjqdjieKd1Cf8wxy/BLnH5U7wJShUR4YROTBsm72ZxoMGxFoUyQwxxCKjE4gDL3luIIYnQ2k/9Kh1LgltJOkYACJvsMu"
    "8swEy7cRmKq0vzH0k0t+sPvWe2+/2c7ZHiD+0RAiGlOqQ3c/6vVVr3/volNfWLns8K++59C9jtx163cM9Q01m81mezqlNNVpO2lN"
    "N6earelmu9lOOq2k3WpPtZN2O2k1O81Wu9lOmq1Ou520OipN004n6bSTdqpTrdKse4ZccKCAtMZI1Gq1RqOHlFZKgczO21az9dSy"
    "x/be9l1UICWaun/BrLV6G72dTlJwUQYoWfoc3LYTXUcLRT6mzFuELDQIHXcHXGtMedy35zzrxSWC7xkX9L0unjJny2XTIIuO61an"
    "xSv0yO4ZXXBhcIbQHGV1MDx18Ycjf0jkaLHI5ZXo+W2Ai5Fz2SrE/BUJS8y2WUACUgka0PUXB+uouOsK6Z/Bs5FCtpP2F3/5yau+"
    "f+Ng7yzHM/MBEoEiSdNPHPrZsenhc/54huyLL7zh/D/c+GtUIm0qQsIKoEStNCkiTXmcs6NFYC6DybtkYX6EECISEoVkGeyeNQGR"
    "Nm5ZmGUoIymgDi196TnfKoBREQCG+mYP9s5aufqVSlQhFitmJPl8kyOXqqNVwaBNfMQA1nA5JmhBC6bdIBsDSR50niOpOe5Kvj6/"
    "ZHxmSRnmdOQZHoiayESf2PM14mQU79wk7t/meIvMYfPNhzRBqoavgmZrK7CSB88lxsJo5O0oD0soSgxL2OnODJccIxRDq0nHpuP8"
    "C9/SWGnVU+t5etmTXzzvkxee+jfgFtY+shRJmar0Kx/4ztT09PnX/mTe7PlJO1WpwgoaayyT6JaDuDmx3nNnyC9SnxJimOY5ZgC6"
    "UOD5QkQAAknR8teXAxgvdUd0sEnzvfXegcbgy+lLlajKmjReywCP/c25oOh7+bBMN/IeOwKWiEXs+8KZ7B2C4rbrEeeHupJTuGUL"
    "3GW+ILhKEMnMALArzMjb0cLHjHz64P9swtBXooRmisBd3tDJ27sBkwUP78DvivxqEAJdUbBQyDtPMPz3+V7N+MaU1+VDvbOuuePv"
    "P/rr96SUSitgSY/8rQoUqUpP//jZJx/11deHX6OUiECDVpSmSqVJolSa6lTpNNW5i4ROtVYpZf2xNiGQmXc7V01zZMuL+XVpTUyy"
    "G8lozdiqdtKyjwuNZ6MJStY6kvFQ/2ytFVjAx1FwCHxDDbBYJwtIxZyim3tjAtlsXcRAWZaV6Tk1kXKCZKCMBSfSB+SC2PKtgX5N"
    "jn6ikSPRme8h7ATHWG+7AwPRayD8oqxYIrJwTn4lY/dBWgiYdhswlJ0DPtBLFJzQRblOCZ6MUI6oQ2i+7QV9sU2SajV7YO65l591"
    "9R2XxlGcapXnm+d3ezYCRIFCafX1D55xyuHfWLV6FRBgJAgRtGdQmNsNYN7YoRAChcg5e+i5DbrSx2+SymBW826iSjTVnpxuTQML"
    "77HsR4Pq9vX0WdZ+4UYjxl8JmJfoPOutiyEf+XAtvRcGhfbSA8qFyUhhVJyfvYn+CIIho940ILCXKDjDkqFCoLkTIzPc4UZtwQCK"
    "MDA8gYK6xRnq+zGC2I0lkScmOSFKmGlFDoMBZqKC2PUq9JklPKLS2mAVjHkp6LwDQRNDKjx0mIDqlZ4v/vKkuQPzd99m7yRNIhkV"
    "xzfmBFEq/dpxpw/0Dp3xl6/XoRZBnGb8c3+Qmp1HQmlFBJgnp1rz02AbI0Bm8Y3Ao7Ctt7ItCEUsOtRqd9olRUaOpPQ3BgkIhADP"
    "SNkd1rzb546ALJWB/JQ9bl2InOrAlhRi+LHlYUq5zp7nw5R/2IgF50tyNDA/38NuliirDoWMo9isM+3acGQiL0+zaZFzHlOM7MMm"
    "SygrjoZ52LJ79xTEkNjv6jAYCCeExVGou4eoEKDAqSEenY08tw/gJCC0fqve0BKs6kBHMmolrRPO/r8rzrhhiw22seJDDCzSEAFE"
    "qtJPH3HyugvXPeW8T05NTvXW+1JIUDpLbNPppkk6SdN9jf4IIqVSRVqpVKVKayOeNKGv7IXnNwRiZieNAgFQAAqQlbhCSCRIgOjo"
    "Tidt21YZg7cNEMsKUYBwElqJKpdk+9F3yPKtyGMzFwNYkW0MykF3fuygZU5wipMj/juNjg0+ZlYkRISGY+eNc81CMvFT6NigZh0L"
    "Ua1Ug8FGdgkw5ikWkjiK/C1r1lIKyFj4hVgAWKjyBSwgmIQlxX+QsJDzqdkJEtgHuOMAkZ9CvLorJaKyoQR5iWUISqt6pT4+PXbc"
    "94+8/Hs3rL9wozRNpIyzE0szthxmBp0H7XbEWrPW/fTZH3nhtaWzZ89ROiWls0UlRaI6ffXBrx/13T223bsW19qdllEItNrtdqfZ"
    "7Ew329PNdrPVaXeSdqvTnGpOTren20mrnXSarcl2p9VRidJpu5O0Os2x0ZHXR16t99WqUVWplJTWSjG2SZjIg66Uzz5HYTx8MEjX"
    "YwTz/GKwVpjc3slGwHkDcyY2Zn0dujlx6L6GPnvehs3xlBurneWMC+v7RICEAuMo9uYA5nsIxFqlrqkIsSCnzHN4QVOQfev51xW7"
    "XbRYP8vIoIK1jANTmFku8XlVlt7OjmcbwMIyMoF8c0ZgIm4MjLqKw28/ChNZ/URYzHZIteqp9b6yasWHvnf4X0+7Zt1566epklKi"
    "/xozjaWQaZpuv3ina3988zcuPOWf/7uyJ+6LRKRUaqQ1tUb1wq/+dZfuyZ5vhUeqiVKddjqdyebE3U/cfval33l9zWsVUUFAFDJk"
    "gYdlIBQ8iDlAnOU0+oF3bkTmyeQxG43lT9/DagJyDHqYNEtsReROvz4Wns39Mwt1fibnEhNyBAVAxFrF82wVtgdq1Hus6TQyxRVi"
    "Kfhj2zIgJ/8hlncfnqSUcxwIeSMUtp/IHYnAT77I5RlUCIe0N6Q1k7EpPAwI82CqcNYD5BhBLjXCVkxMIOuSZhxHLFVpX2Ng6fLn"
    "jv32QSveeDmKpFJpyXtjpNH5sxZc+JW/nn7cObqjJibH40ocSTk5Obn9xjvvsvluSZqkKs2kw8ZjNE2VVqlKlVK502hivixN09R8"
    "lVJaKwOJVmSlp9a7YPbCw/c6+qcnXQgKVaogRUvg84SE+XNpJ21EL0aBJa7ltQ15bsTAXMn4o0W7KB1OYsE9Cgb/LEzMuuA5ITEV"
    "EMqgGbIJFl6QJXGRHwKCQFGNa8AyioV93f2NASKNvJeCwsypwA7K5YiEnHwRerFR0HyTJ+cCj8VcWupg2RnB8VICKkTYE0uY90XJ"
    "LCs6kGYGFD5mBo7I0aTchy+zzAdESFXS3zvw/KtLP3D6wS+99kIURYlKiUp6FuMNYahsJx7yuSvP+s9OW+42PLla6bQaVyemx8kQ"
    "n7J60ZF/BAppCNAZKzqKZJTxglBmvGjMFPEaCJC01mma7LrlHtst3n6qORVhJRaxv468JzAxPYEeJTe4HpkbZ14KUh6AQ5y9gFDK"
    "K7Nq3+xWtWchY/Cjf/ViwKx3/LvgI/eRhNCRhLTWAkU1rvL1I2xpM9AzSBT0P95UOLcfpRJmfWjri0H+HBSxf/KaaQcmIIe5WW1Y"
    "zKV0QJ4vtCTyh8s8lKyEMxQCytmWRv/+YSKE/PwM9gMgpirt7xlYsuLZo7/zvmdefjKOovwewEB4kA3oUSRp8rbNdrjizH99+8Nn"
    "VqNqgq3Hlz14/d1XShkJEEprndtD6EzyrnVOm6M8jpUCc2KLL6IDT4UQIKi3t6dea5S3cAgAMDk9kQcP+nHNmf5IF+panhaW+fcg"
    "QUGiYm3imPkuWmsPi8x40jMf62STjxBadR+dn8XhszWIIhnXaz32xgBOiZo/e4GmnA5MPiLCXckxYMlSAMfbb+2wW+JjBO4U583u"
    "ee1jXe78RLKQII7hQJirnIrecrYYt571VDIDQAyHCshadGDEPXYB27zVNE37ewaWv/HykV/f/87Hbo7jOFUpcYKUPxs3GV6RiE46"
    "4kv//PEdx+zzkbSpP376B3537a+anak4ys74yHOZzhwUhZBCIPOZtix9nZmQakVEURTd9ehtDz7xQFytzJo9u7fRS8xxzJIlpZCp"
    "SkenRqwbdv6yPR4362AJ/GLQnxryTovs6Uaeu7hvc1WUnqBHKkanFLYFt+fe7vnms79s3ocmXYmrjVrDgYQIkQX8Z/XNZUQzHmlP"
    "rl/1PN4sUdL2pkChVQx3gEaPM+e1sOi1WbxPReauwoFFbT32uEQw2AkhWZGNbfIEAQKcmcvh+Its1wXcRE05OR4AIVVpT613fHL8"
    "uDPef/ZJPz9yr/8zZre5R0t+QOWflok6VVqtv3DDn5z0m2P2/uh5f/vhF3500nlX/OgdO+09b2Bhb723t95fldV6pdHXGOxt9PbW"
    "e2qVRq1Si2RciSvVSqUSVaMoNsbr/B100s4dj9x86m8+T0BJR80dmB/JSGmFnAye/2q2p8enxyIRERFocqYPDhxmnCxeJZcYJlvY"
    "3makInA7OGLhLr6pP3M8RncIg9s9bNCAJXNV5HlKTqepSdcq9Ua1wa/9yP6QWX1zJDqRB3hEkZLjl6BgNkIFmqDjE1lCcX4Ia0Io"
    "Iep4qKYP1DilcwHDsBAscrtJQP9wwVKFgxdCEg6Ji6cTP1zQC3RlJ4RSaa1eVzr97LnHv/TqC1885luGPeECvArWmlJIpRQA7LTl"
    "Ln/+3t/vfuyO317786vuvHSsOV2PoojiNNEq1VLKSrVSqcSRiCMRSxHFcVyr1mqVWiQrsYyqlVq9Vq/GlUjGiUpeWv7is8ufiOtR"
    "o6c+8cb4goG1IDOIxuLQcGxydHxyxLnfcYqOHQaF8nUqUg882THDt4nJZHK2hHNp9IBOpyLxwIgSGQzzcWODoJy+azXBiJq0Sabi"
    "pWJkd8O8oXnVSlWTBrbOEJjTY7DU0JuLMCcU5hxsnX7YQJf4FZIPQoj8R8sjGNBrLwoSIvDdcgkY9sjIc6ze8S3HiDuLFjmh5A2T"
    "qbC/bW43es8MFSmBsq8x+MO/fe/x5x/+wSd+ttacddM0EUIWfOUdOgQEZhvsts2eu22z59MvPXn57X/5x51XrVjxUiWuNPp7IimN"
    "kiZNk4Q6JtIvd0zMBDQ6O2EICSREtUoDNWhFIsZ15i+CEng9++c1Y6vGJ8fjKCa+uD3Nffia/Yg9AN+XksU7AHjmtXw8AmFJZv3p"
    "rMrAg5PZieZI+OR6QUuYYPvTaDkGegYbtR52/JGwPKq5g/N76r3a+Hhwzw8rg/Djtoj8URUwhxjiKoi8l6dgapWpVNAzwbCDOweg"
    "QRjcyWM1eE6Cd28xR0n2UNm374Y5eSI1z40LycOLXEYn+fyczByTjACOBhqzrr/7uoO/us8tD/07iuJiXLOvjAZT0But4+aLtvz2"
    "h35w3Zm3nP2Zn+2y1e5EevXEqonmhFY6ElEcVSITMhlXqnGtXmk0qj2NWm9fvb+v1j9QG+yrDzQaPSCyDjquVzZaZ5OyGToaNchL"
    "ry9rdZooRGEIUBpnykydw29ZJB4SEwtSKCQBF0ML4TnPTISKM4tc64Ru9Tl7NySPHq+VGuqfbUznLXdQ2ON1Vv/sgZ5BpZQdGjjP"
    "kfKxiDfby/OOyKef2eacuZ94DbbR4XpB60hY8DJhFHAKn6Ltai3waS19ys6swE+dSqQC5EiI5LNCmcrHj4NjCk/IA0WN1i5N06H+"
    "WW8Mv/7h7x7xg4u+1Unbhj3qQcnEnc8yhYrJNk3SZP6sBR8+4ITLz7jhmjNv/taHztxhk12AcHh8eHx6NKXE5GqbT0qD1qRTpVKt"
    "FSlFyqgfEVFT2lvv33CtjSFkBzoq4Qsrl2qthBAZk5MPWwJnIPLZ7OTDz84FIg8GcWHX7NOkELVEN31ibmcOkrLHsODOLhQwNMk3"
    "icuB7VTruYPzAYAPfCNbtPXUe9eas/ZLr79QxRqisqsfQ8ilmB7i4F9PDgfF05P8O5F3PU52bd8WsaEgsyoC5ugPoTyTqYtsuUXc"
    "48hNCB0bD5l0L1OlOJDUu72Jcc65BZcPdHoXs2mLq5U6oP7pFWfd88TtX/vQd3fbei/TFQgURpOOhMVKQ6BAgQb7FEJsvmjLzRdt"
    "edJhX3j25afuefLOu5+87bEXHlm55pVWsylBVqvVSlwxBqYZfTr7tIWU2EnTdfrXWTBrLbARbDbMJL8Yn391CUrJQFUM2zL0Fzoz"
    "5id30mGoJSLWLCNrC7hsnVgjZWFiYr5XAUeSh5FZxRVikaTpiiukteeta/++QQsiU4QorSIZrb9wozsfvwU9uZTz3CoKQf2Ds3Dc"
    "WhYeWTqeX8Oh9WMgf6DAfAqR7x5wxKoyb7ZcvWbRZCyUpMR9ELEAhCJy2aVX01I3aUIXwZkLFc62tQaAod45Dz33wDHfft8H9v/I"
    "KUd/fd7QQiIyOXxd1B3ZazfQpIkUkDJavGjLxYu2/Oh7PzE8PvzUS4899Mz9Dz93/5JXnnltZOXo+LjxX6nWqlEUCRSaNApsTrc2"
    "mLNJX71PaWXsEPknaDDQF15ZWolirQk88iAj3LAhhkV00Kn6sOAg6C1aloqN6I0SWPgk5XNa5KIgKKfvuoiUEq57Zp6fb1QhxKIF"
    "6zOSZUaHdlt0/QUb5pvUMf6cjSk7AgmDBA4Oz3hBOeBTOf0BLDLSqa8YIz/tHRhq4CEP+WGMjMTlX0a8mPT1w0Cl4gbGZ0EWLk/M"
    "DsdvwcF7psyqKb8/rQAWUpX01HuB6PfX/fqm+/71qcNP+cB+H6lV6vw24I8sIIcLFGAaZVLGX2NW/6w9tt57j633BoDxqbFlry59"
    "7PlHH1/6yNMvP/HyqmUjU8PtpBVjpQd7007ytk23K35Pexy+umblS6+9EMsKac1Gj/lBTlx6njMNNcuE81AFn0/unDp4vV9CTOcx"
    "5iyxBAJYiXwWArcqcsRqP8hAk67G1fXmbhBsFbMBsitxg4UbcVkMMJt1D6q3RUpmw4RURiNzGX4UTqvQltdo1xv6P8jLDykllTJL"
    "UAotONyrtRhDQRWGnGxbZCyg16wzByxE75LInAJzbokLtCFPnpJPE1CRQsBZ/bNXj6869Tef+9stfzrxfZ89aPcjqpWa2QbGG6tE"
    "MsSILQIFSAACG6eHAvt7BrbdePttN94e9get9evDry155Zlnlj/xr3v/8eBT9zf6GjtuuUsph9fUV0uWPzM8vqa33quUDkcExBN9"
    "goPY+oQw4I2YKMDS0230iI2ecadmntATKhGBFQk+183hSMSThTJSdeAxDqBU2t8zuHDO2sHsXyA73tdfuEFPvTfPlzUEaPJaTys4"
    "dr6zmI2BgEKuQp56mSfbkZsnI9cqYmAn4TIOiDzyJTOSdgCpyOccWEJwYR00ZQGdvlNlOMdh7hNB1gb4lC63SSizVOAoYHDKkvV+"
    "wYzMlyoVxZWh/tlPL3vy0z/66IFf2usv//n9xPS4mfsavkMgeqJiZhdmAnkppEBBmlSmG0iFEAvnrPWObfc58X2f227THSc6Ywvm"
    "rbXF+lsX9RT2/Ty85IFEJ6brAPSzOTFksTB9BL/onOUVBcdKvjf4eApc0ZQFxgT5Yo7KwZyA+PajXBJh1yc43J2VCYidpDN/aMG8"
    "wfnBoEcQa+fXm7f+3MG5SZpgoIEGnxTlz1UJyioJhxGj7crRxpNQqI5zOpos/RfzlCpOQ4LCKQGcvRRQev30eAwja/08Pj/iDHMj"
    "Co/Yjb4U1YsnBo/yh4ghA4yR8iyLVqWqFtcH+2Y9s/zJU8775Hu/uMdPLjnj+RXPGb6DibTIMVPC0gmaB6lBvhkkEaUqVakanx6/"
    "+d7/6DZsNH/xnP65WuviDjDp048+/0AkI0DMNoCH8/hB7z4pDH34klzWTpFDiaxPYOIN0yN528xLrkUsOrEQK9dZpQZhx4wIiCJR"
    "yaL5G9QqVZMCk3tvkcgVxqiJ+nsH152/qJN2EETp8JMcdG/5NLwnD5enBkLeHHhMCuLe7byusfVRzp3KOUhIXBnEvmExLTXDg3Ov"
    "MRvx4KcIuz+2lRSG+wK5TSL53BNObiwwURmtlArMSAuaGnV8vdoY7Bl66bUXz/zLaQd9be8Tzjn2n/dcNT41Zqg/xsstJ8CFm809"
    "7MK5ICO5dPmzL7/2YhxFu265OyIaEgR434SkkGvGVz259LGqqLU7bQIdyxhZF+XrWXN6OM8MIn+2hs7Cw8E1BqPORsvZMy+Xa0PJ"
    "sACJhWKSJUA7UDpPEyaPB5d/YFrrTdbdDDIKQn63kFGEEZkxgYiizRZtecvDN/KRKeXpZzaqjm9KnAEYspmOiNxniOmyGHZO/ODg"
    "WkvP5cgaNXsEceZ5GJqvBo23e0EeH8++waCmdEu/oE7I7zNPI1WsGXhin0U/8sGne8Faa41QrdXq9UYn6Vx399XX33P1onkb7r71"
    "O/ff+cAdN99loHcIIDNkV1rl4kl+I/mWarng/sGl/5tOJ/r6+/d4296QG2B5WBVpAPHEssdeXf1qNaoP9Qy1aGr18KpGtadWq5u0"
    "Mh4jq61bXi5kRUIsq5GII2gO/8OgW0Pyul/rfMk4mejFlrpvRexSzBFLZH6O+RVNBFKILdbfCgrsoYj7bwHA5uttyfEX3/E9iM/y"
    "1LY+PsWEX/lxQO6ApaCN4u4ybopEEDbzAMFgijuKlWvkc8jG0xdlqo4czOOmcECBjh68DVyMvS9LiwBnjVFidu+JD1m/jEBEihQA"
    "9tcHiGjFqhV//Pdv/nLT79efv8FOW+z+jm322W6THddfawOruzf8h5xZjaD5C0fDsbvryds7nWTzdTfectE2hrVRFOsBwP3P/k/H"
    "qRKd73/6J5utu/kl//nDtXdc+eKaF2rVej2qG2uWQGFI+Q1JWOwR3LiGEddcDDMfa+atAAMbCTTr7MhR6iy7gOPbhL5o0WKL9pVo"
    "req1nk3WXhxuSITIvnaDPGy+aKueam+qlCUlIU/hztXi6O3jEJ9iOLtdRE4rzIJ8MxqNF/lCTFvt0ntzayUim2vmI0+OVuULLR2R"
    "kGnvctzKzQssARZ4FisRt67k0mw21MSAeO7URkQUflo+Ym13FvkImwJFBJVKtV6rE9HKNa/87aY/XXrzRYM9szdZZ7O3bfL2HRbv"
    "suUG2647b1ElrvDtqLSiPGcwknLl8PLHlz0sKmLP7fauV+teWnD+giVKInrg2XsAab21F+2+zV79jf5vfPT7nzziC9fcecXF//nd"
    "k0sfrcX1eqORKsX8M9macPcmM0lwLS0Te3B+i5eASLzm4MRRZvvHqwPgPPx8uE/gXBQchwIRE9VZMGeddeaux4eANifY3vUCADZZ"
    "d/E6c9d76fUXq3GVXIGPQd1eIiP1osKYdSLZCsNxolh9jBTyOljoHxXV91gMgckTlgqeoIx2zFlF/EhmUVf+oI9C/iM6k1DHZAbn"
    "4kgWbfWybchJY4OGnP8mnzQwbS0iAKUqRcBKVK3GNSDoJK2Hnrv33qfvvPC6Xw32Da07b9Hm62615YbbbrnBNuvNXzRvaH4l9gSv"
    "tz9606vDK/r6Bt6zy0EA4QDRTgBeWPncE0se1QnssMmu/Y3+JE2IYHb/nOPf+8lj9vnQZTf/5ffXnf/My0/2NfqjKEqV8jIW0BuV"
    "uZI8MPj0jVaB/D3iBvusKeDguBt3YdCQkD+GsQ6ctjQQKNqd1qbrLh7oHdRaoT8EjHhFqJTuqfVuvv5WS195rl6tZ9N05OljUIK4"
    "8AQ7vixzH1RicfD+GkP0/4ysvYBHGcqtenzTaWaxqv3Rm4XwyNPulGUSFtBSj+MVEH6wJPLeap2Auz8U7oQSNpk3mnDkCwqIqESk"
    "tM5wWyEb9V6BgojanfZTLzz+yHMP0k1UrdQHGwML56y1yaLFm667xUYLN9lg4cbrzl/vutuumR5pb7nVdtttvAMRCZTBUMkwf+54"
    "5OY1E6sr1co7376fgcaFNK5E1Kj1fOS9Jx6+99G/u/ZXF1x73vDk8EDvoBGmsfdAFEr9rIc3lo6IrN+PJ/zzRj7OSMKZyVIBiHWJ"
    "lHyxAfdfNxDCthu93fBCpG9vGnkLhDSA2G6T7a++4/LQ7Z33FgzOISdsQ+dc7o9uQvcMKCclEz+x7UmbRZkChgMYwsBG1XvI1m/L"
    "zep5yBQFZ6+3/gE9Fy3kzyvj7RQbfw8XDufB7NhjpHleopGzoLf9K3GOd4ZBgNJaQ8Zvq1XqjVqPcTWcak4+9eITjy57WKMWqWjI"
    "3tkDc8cmxyjBd2yzb61SY4ZFXL4gAODWJ27UoBYMrLPrFntmjTKhwMwvUWvV3xj4wjFfO2j3w7//52/9695r67VGHFXTNDWiGSiS"
    "xKCQHOM9/ZyaU+JCUqJ19Pt7q4lBawHposIZMpl9hiILDN/RumyQh6uK4AcBwNs22aFWreVND9gG1fc6JPAhEPva0OMbkweglmDz"
    "gBhGggPTiSKjVjP7IAzyNzFgUZSYo+eYcTZhIIfiM38LrsFHD6HFMspHIAhhqQ7kJQT51O2SkZIrWfMjDDGf3LkxIpF1p8ifpjGJ"
    "ICIhZb3aGKgPDdVn9/cOYgSvj63UMhmc1X/gbodAYJaYvWAtpXxl9fIHn7lftfT2G+48b3B+mqrAeFPKSGudqmTjdTf7w9cv++nn"
    "Luit909MjUUyCq7VYuSaj6VhV7tk4Gg2FBzN7HzLkZTyFo0F+XmAuwNgE9WZP3uhgYAQBSDxCtaTz0kUALDFBtssmL12p9PyKWCM"
    "2e/VauwNcX884rS9LLgKrMEl49kTN8pgXtRsxZA3cLKDPl9oWfJUS/PD0UnHHI2Fe+gF3CN4E8Wk5TBnRuRMIe2L70NStpvOIblD"
    "jFsROKcRctlOyE89QuFcMcwyVWlKpESEjUajpaZ32mqXt228vdY6k4CxH21e9o33/2vlK68KEgfucTAAaSMOCiSOCFJEJrrm6H0/"
    "dPUPbtpj63euGV0lQFjWBrlSNUesETjyzbyqSkLkvEdFQMFohRhZkRyTiPwmg4IpKIBA0Wq3tlx/m7kDc5VWwtjwsrmdCFpGrfVA"
    "z+BW62/b6rSFpWYVc5gKE1n0LEW926zoqhlM8KCE54p+bYE2FgYsJcO4VrMz259N+81AMJ8mQtckIMvAYp0dFkJ+fCGgC1Tnjzw4"
    "5Sm8Q0LzOczmo8jnoj7okK8oY4jul1WWSIL5jDnzUkFA6LTSw/c81shrCt4CZHqJG+66NukkW2yw9UF7HJ5LJYmoZN5slnuqkg3X"
    "2uivp11z8lGnTrUmNWRUVrSaOM09wUt8PMI6KbRmQRcNZIFv9C8bYg74wCyCSwi6mKh0x8W7AY/qYMtQFHhRCgD23GZvbWWCCEFt"
    "4i818t3SqDigxxLfNaAipF429PaPLGZ0TuQV3F5bySaYBL4zG5EVaubqGXBFJHeUwED55u/53D+bx/F5KZehDplpZdjOorI+hNhO"
    "zbcFefwij09hueReE0+i1W5usPZG79nlYOP4ECw1TSSEWPLKMw+9cF9ffx/GdP/T9yCicalQWVa2T3UwxGwRGUL1Nz9yxi9O+X0k"
    "ok7ajmXMwnVCXBz99wWFwtUJaShkn+e2o95ZYm3egM35KdAF5mVerVrfeYvdgHfkzA1DFNafAIBdt95zsHdI6xTLZj9eRJqFZJGV"
    "LMRkQugR9zwqBwWySiirY6yTDB8rFqgHbBNwsb0zkfFY41gUQTqkmRiiZDNHS/y1ic+AiXzDE5qpbCKG9GKRlYYlF4nPNLIUq/zo"
    "y//ZyYsBm53p9+/7gaHeIaVVkf9jkgb+df+1Y+2Rvr7eF1574ZjvHHjUN997y4P/Mql7Sqlub0GgAMBUpYfvdcxfvnPtrL45E5Pj"
    "kYyoRCaSrwOiAGfgjmNQjHTEkN/j+V5hkSJcAKIIBIp20lp/wYZbLNoqk1wX30v4z0IQ0cZrb7bpeotbaZvz3h2uEhz27kDPY07R"
    "Kw0AArQzJ65it5CZgtcFB85z3gIPTGAEGLJ8qDLckUkp8n4ZkEt+yjYI8uirYP2TPwQuydDqIpgJBGQQttJOReo7d/FBBWMWksdE"
    "w45K5g0u/MA+H6EgL4WcAma6PXXNHX+vylpKaSziWFRve+TmD555xEk/O+75V56LosgAiFS2nRHRmJzusHjnv51+3XrzNxgbH41E"
    "jEFXZgeYXgWJhfAO+3RD1kg5iuyLMAkKO8I2AElr16326K33pcahrLC7BBReilGH7b71Xq2kJUAw0S1X1VqetqcKRvRZfG+2tCk4"
    "vkM0wH7ULkmBHM5MgeM7J3ihf6x6lvTkJyoTeX1OyH8O6zViDwGDBD7Ckka8oLTxBPzkl3g5dyVXChM3eMFcAcI3hm9bDAAoo2hi"
    "euKgPY5ce+662pL72bmrQaPA2x6++fFnHq1HdTKiAqCBvoGeWu9Vd1x28Kl7//TSszpJO5JSk8KSro+MTi1Jk03XXXzpGddtss5m"
    "k83xKIq8AzMvOdGLAfX7PDZ+h2DSU/AlwKBfxgLrjR93RBLFO7d7dzcxHxEvgdAL83jHtvtUZVVlFhuEnrdvSP3JP0i/8vXTDPyN"
    "4eBNCn1TySMQQg6tWpwE/RU3U46MZzdEloPIfORDPWdJMBXxuRwwwxju2lGYbRT5ocxqmMgzM7Ps4tx4jieieN76zrmYQeE+/1Qg"
    "JklnzuCcEw76NDH2dQFmgMtuusiwI4kABIAARVqlqq8y0Gy1v/e7bxzxjf2fevExEwnuUHG/IotknKTpevPX/9sZ12287maT0xOR"
    "jEgXDGXDAs9q4NG5tPmLM3/WzPaRudJnNC6HG4FNZ8pXh+gk7YWz1tlp8a552VYkFpgDnoIKD4HgbZvsuNHCjdvtFmR/k0qkJhSE"
    "OAXmm1lAIWfwhuZTPN3LriciX4BgBx4WKkRL9HaqXizIl0K0JY8gx4LDBnMC9eRlyPKtkO0bAl6Yom/DXmxsvFdE5JOivek199rA"
    "wnEB9tp12d7eHUgEAsX45NhR+3xww4UbK62sfS+xaZEU8qmXnrj9qZv7+vuUzhNCRbb+0lRJIefOnffQC/e9/7v7X37rn6WUJp2+"
    "dCNFMkpVuvbcdf/09cvnDy1stqelkNaZ1FPYo29xRbl4FRE9ZW/uNFrA87xTHC0k74ZItuCTKJrt5vaLd54zOE9r3a3cFgAQsukR"
    "U5X21Hresd27mq3pHD9GO14l/0MLPoM8SY2f32EVxA4xNimiHNZ0baubhSGzifG8UbHgMxGoDjzrREPfzf7bT/FECue6heGChVvQ"
    "+YmgZ5JL6MtnmC6Tc12Rc4NzBJ2YwIYVRkzJQxkBzfb2Nk3MPmrR7rQXzF77Ewd/zpkxZveAV1H88d+/nWhOxFFk4uEkRpGII/Ob"
    "SIpYkKCBnsE00Z/78Qnf/u2XUAiReYYif+rm5Zs9sMFaG//6y3+OMEpVaqtrz6mBneDlPUIgt/BAtEBBTxw25S4inIL07h0PgDz+"
    "rNQKR1DuxFO8/Q/Y5aBqtZbteyyuBntlW21yvunDmNOSyRG6B2GLBvSPuqA18GssF+HHA9OCUS54mhs3x0MOOjAOj5/8Azyb0Izz"
    "3MYmti0Zao68FQlfF3qJA7lW1LPDwIKLo33PyFLtrCrLr/8xEnJiauKEg05ae846oftDfjgJIZa9+vzVN13eVxkghJaaGpscGRkb"
    "Hh5ZMzwxPDI9PDI1Mjo9OplMTnUmx0fHMRHn//2nn/7x/023JvORQsj7yrx102SnLXY969M/n25POcilWAeTNdLLMmXIjLwx8Ddn"
    "JinWqAHz5iufBWcMZZ83iQI6aXuteWvvvd27Hf5T5lURZeWFJ0DPTDh22Gznxetv+cyLT9WrNct2DLgAVlUBnCyAWGx0A/mU89dg"
    "DbTHb2bR3FgySCEnX3RkGmQkWz/bnLzIPHugMjlnGGqQT4vtvIuR1z0doG9JUBx7oRV+5FZLCF7wGXsHFAr6HS/SGVR6KVoOABeI"
    "0+2pLTba6vgDP5lldxfEq1pTJMVF/75gZGK4r9bfUs3tN99l7cF1qnG9vzHQW++t1Xt66z199f7eRm+9Wl89vOqy/15839K7/nH/"
    "30fOHv71KX8d6puVKiWEyFFiVwVIGaVpetS+H3zixUfPv/qns/pmJ2kC5QpJYkpbqz0tOu+gNT3Mqbi8WfXJJOT0WwLlRGvsoN2P"
    "mDc4X6n8IGDZyIwNiiXtIwIolVbi2r5vf88jSx7sqfemOuGLmzzeKlp/JfQDah1IzlNugvOjBF5A1uIwz1HyfbU8KY83BfIDRcGP"
    "tvWAM+IGMBRQ5pwvOjFlQRnYDd5VUvh3ZJ0hM4UDeXmPtqYh8KQ3BdCBkDzJI3mTSkTRSTqnHP2NvkZ/qlKJsjDlJCnEyuEV1953"
    "RaNR64l7fvj53++/y/twRtz28L2POeX8E6+849LbH7r1uNMP+/PpVw/0DCmlTEUUENrMFfH1D333gafvfWTJA731PkUKC28FXRAG"
    "Z6uHom5AZ5nGGZXINYM2h8wSiQGBSIjooN0OLxyb4Wcnymbe5iMRAHDQ7of1NwaUSt0guGCTyIeZQVw5EXntvZe9HuJIHnBJIZTo"
    "egMbD5KjJOhbboMvUeXdu1dnIBvgQteJNPPAR+hChw7A/RKei03NDQz3yRuhu7uNGCGKPClJIIuwVZHEaHR85D27vO/QPY9USmWr"
    "v+j4i/jbf/zi1TWv9PX1/uYbF71nl4O00korpWzwUpqFdatUKdVJOgBwwns/Ax3srw7e9egdH/vBsVOtSc/h1EMNkABqlfo5nz6v"
    "XmlorYS9JhyRDPz0XGQzVR9T7Ka5ohIBVs4YQURsdZobr7PZblvvRUbyj13BeFH8CNFtZb3lhtvutHiX6daUyI4T9IwVSrw3yXe4"
    "cqF9SHZeb8FODGls6EzDbIPno7ToM/u7nMYIAS2WDeDZhtAOnSCPchKCXuQL/7uNdwsZ2yHPhYLaDXmUNA/kIv9dF3uDwCVBJEln"
    "sGfwW8d931lsOhoCAoAmLYVYsvzZv/7zT6pN+27/3t22ekcn6RjJiMvdEFKiFOY/UkRRpLXua/Q3aj1t1Z4zNO+me//9hZ+fiCLU"
    "PPDxU6rSrTbc9pRjvjY+NSZdcDKBZ9LNcl7JRr9hjsN5sl829MxpcKGLittAAsV0e/o9Ox/UU+tROvVcDQvFtEAsozabcQlpBDx4"
    "z/cnOvVaW2TTOTd4IiwQ8jz7TAyZdMGAz8pKGMGY5y34yTThrNf/GeG+8FM6clK2F89u6yMfWseu1EWO7HJGHI+Oo6B/9Wbnjj1E"
    "NiPH0Zm8HYghtErEDY6kkOPT41/+0Gkbr72ZUsqML3PAlrjB23l///H41CglsM7sdbVxgMvdMbQ2kketQWuttE61UkmaCCEeWfLw"
    "8OiaWMapTubPnX/lnZeedfG3pZBaKyqcO4aDoLT6+EEn7bB4l8npcekQoSDoMZv6kceGQAj46cVHT1njiyw+Nds2BEqnQ31DR+59"
    "bL6++fANAjhIMHFVQJrMZgf77/y+9eav307amPX1JpmqiMoGiR1MF4ClbwJzvyOvcSRL8AQKeXgQAOvop4QTcLocefQ5d9KDDW8u"
    "ShO88gJzliXxYQH6U4KQ10Dutkd75jMnO2dejWDTvdgYOx9755YZ3CHcniY+DzyS8ejkyHt3P+RjB35CKZUhHr4SX5OWQj7w7P+u"
    "vefywcGBSIp7H79LCFGJKkIIAQIRhZBZGp+QURRFUSxlVI2rT7742GkXfj0WMQGBxJTUrN4551565tV3XCZlpIuUIcxoJtVK7WvH"
    "fTd3igtwBiQGSCC/CgGEABejUcBuiPVl5A1iCQCFkJPNid232XuzdTfXWmeeZYhQyhVFiPgniN72RERQSs3un3PArgf99przZsWz"
    "UsOq5SiL5/lPwaGLXuif/w8+XoTFTcI1mBSSPphTnhtOoQ2FtAcp7/p50AMWUX7rWeAohhaiZD1+uRKaO3wHT4LsleVE98QDz/nN"
    "yoyQvSBpZ6qGHrgkULTbrYWz1z7zxHNZmhvDjfO2S5H64d++m6SdalTp6em5+9E7vvGbL5x0xBcHe4c6aSdJk07aaXearU671WmO"
    "T4+OTYy12u0nX3r077dfMjY+Wm/UFSkDA2qlG5W+r51/8tYbbrvR2ptZvJV9JmS83ffc5p2HvuOoy2/962DPYGqqEatusZQxq2V0"
    "slLf8pNK0ixZBgz4sCkQ4fv3/mC27S0SgOVoD/ozAmSWTSgQNSkh5CNLHzj86/tJEeU0dG77yUwxQtfmwNGUiioCTtsInG26qWnJ"
    "9bEUVt080YhK8BjWCnDFvoNayyxv3XPNh7VURHw4PuQ0wV7bz2cAjmDMwWWHNQGD2lgTzDJUrDWUmJwa//03Lz1g54OzCCbyNwmg"
    "0iqS8m+3XHTyz04c6BnMaGEKJ1vjC+YtnNs/f6o53e60U0oS1Wkn7SRN0k6SJgoQRASNWo8UEZHOaaTZ9Hd0fOSdb9/3r9/9R37o"
    "h1CeJi2EeO6lp973lb1Lst+ZP6KzHPCPJ6Zy5Pau/gdMeRw3gRCi2Wlutu4W151zay2ucVecYBjiJJHkTXYJXACS+Tqhtd52o7fv"
    "tPluU61JgdJGcBac+REKMWHeSCk4z7kxCRH34vWJUyVSKgQIbYMC20byKqis8MBwJoVh1FoYAe/FzKAjZmC5EJbhAARAZURvNhcl"
    "6+NOTgdnhWp5P27goxzvQuQdUSSjkfHhzxz1Jbf6gXGc0BY/4tXhlWdddHo9qqlUZStZUH9jcGxi/NnlT7+y6qU1o2+MT4522m3U"
    "WJWV3nrv0MDgrMGhwd4hAZKUthoX80kqnQ4NDt36+H//8t/fSyGp8D5NM6qV3mzRlgftccTE9JiUEjgAhEi+AUTRcjsgqoer3+vq"
    "sv67nbQ+uP9H65W60QB5O4psyDnjAiGEZN7gcDOF1Af3O15rhQEb1RMckvU0dIyFnH7GrEj4EiTEQvqWF1xEvjsfka8KspbDZToy"
    "tCM3ooIMz1Oj2IsUw648YHFT4HLgDeu7sl69mw6cq5INF8+LOvS1IcEGc1NfzFhoI+PD++1y4Ff/79tu6Fty7xIinnnxaSteeymO"
    "qqS1dWBSpKSQtbhWiWtRJc72TxazqhXpVKepUlmnkhtjZS9AAIHu6en7yRXfX7lquWHRAwugZX0dfeqwkwf7hlKtcn8DvyslYlS/"
    "ArffGMR6IXVUmLJmF2+z3dpw4aaH7nk0UZnrfUhmRQQUJU4/xPVf2VT43Tu+d+sNt5tuTwsh+UJjibs+6I5lHCY/IAf8BDsgbqiR"
    "zYaCgYDwAFhbXxBiQHmn4hApSDEq52r7mTdl0JknvCYKNKszygAolBg4Z67cQZ+YhZ71NytYLoEUcrI5scHam5z7uV9LjIHxjZG1"
    "5qb4+c8DN1x+y8WzBmcppZwBRdaVGG857YKaEbgul+yVl39bykEapXVFVFe+8cq5l51pqz7rnWH1VUrrTdZZfMAuh0xMjZk9Rs7W"
    "3lLRPZGXvxYdXZ14A1B4tEbhcNQ+/zfUO6RUagSfJXii511BIrx2sGRVaK2qldoH9z++nbREXjURawDy+4xT3MiFRVFODcdukyJn"
    "O4GeZCwoLa0FP7pSmjgPlaAYrUGl41oqj5jtAnj6rD4GGnipuRBYI5Z8y1xHQTknD4DHkFIQPp9ntjtIiACklJ2009vo++1X/jx/"
    "cIHSSghZLBQNP2d4cs3pf/pKJY4t6sgpiYHQCNGH1ih7vCgAhbWURUAAgYCYJGlfZfCy/17y2JKHpJCadMkQnAAAPvyej9fiutmB"
    "9mpAKASDhvxYHupM4FeW5KoMQMBO0l5r7trHvus4IjKIpSeE8I5EtwZEyafP3dVN/yAkER36jqM2XGuTZnta5N89N7sn4DqYol7e"
    "1dPBAItKQH7fqzSQHhoUErlSDinIiaTSUNMQT8VS0RgEDQ0fuVHAtws2Rzf5dxkvggfvsjYmiw7hOL/PYwIAAUIppbT+5Sl/2GbD"
    "7VKVSiEYJTT7zLQmTRoRv3PBV5e89Gy9WtNAIjJQpxAoBEopIimiWEaRiKSIIhlFIopEHIsoknEsIykiA4wKIZENOkyNJAAlCHPu"
    "/vzys6EgASE2Ud1ukx333HbvqdaklAL5ljc2Dd7J6F37FinGYEBIXB8FUsiJ5sT79/7AWrPX1loLED5n20ORuHe0gDJiQ+FzI5Ox"
    "euy7PzzVmsqgJfIGsuQpSDiJkdMrkD8jLOjG/Skse1Ac/rFACdcBYGAjwfXxnorBaW6oKwEi27HOSKQ4LGAhyYF7NM54kYTzfb4l"
    "qFSU5mEJQhDRVHPyRyf94l3bvydNEwP7hAMjRK1VJKO//OcPf73+j0P1WToFQkoo6UCnnbam29PTzanJ5vhEc2xsanR0cmR0cnhk"
    "fHhkYs3I+JrhieGRieHRiZGxydHxqfHJ1kSzNd1pd1SqBIhKVKlX6rW4NjUxNTY2Rpr6+wb+8/D1Dz53nxSCeORy/qqMLOHod32Y"
    "yBqakcsFRY8qXQhw6lJN+qB1J23PHZz/kQM+QZ7jTiF3CwLDHuMMF7LkuB+ipcALIjp23+P+9M8LRibWRFFsEg2sxYGlNNkxDjmw"
    "0uPQIM8ft5QngSWSaiIvNNmZ96CPyBMQlgz6ED2LUQekoOe/HkRBMNd5K8JxPn4coXO7kvwRSnkGa5aDEnwCBV2fCwgDBzQjkjFv"
    "HZ0e/8GJPz56nw+lKs1YBhhuLk0qiqInlz122h++Ojg02EmTicmJefPmzZ41N4riqqwIjCpRpRpXa3GjGlerlUq9Uq9VatVKrVap"
    "Vyv1WqVaq9RrlVqtUq9VqtW4Vs3+pF6JYiFkqtKXXn3hj9dd+K/7/zE0ODQ5PfbHf/16+013yoFdT/BvYKK9tn3X+gs3Wrl6RSWq"
    "6IJEkSzbD/kchiGh/Kki12eQFNHo5PAJB3923XmLMhYguvkKX4AubtiOKYPQ5gLe7jhYRiv8w0vOOPvi0+b0z0l06oW2c1mbCwV2"
    "sLzfyENYqYdsUA/NIAvLOs95LG1nbTAlV5NAYNTqISRFOgsG7s++oysf9ImuBhDIlNToWaNwhrbTohECFtJf0Ov+kXB0YuS0E876"
    "9KGnlPg8M+QKBUy2Jg/5yrueffmpem+tGtU/8b6TD9r9sNn9syMZR1EkTY7MWzlq3+zXZ3720avu+ltvtV9Det33b99svS3dNNrD"
    "EpWU0bcu/NL5V/109sDs1OR0FIMELN2Y+PXvO1O5EJisEEjSZFb/rH//5J65A/M0aYGitP1yUlO3qk1EUqmHEYbZNwIEEX3kgBPX"
    "m7eolbSFkzMREaEmBlOSN8aymhsIQuSp4FrFmC7IiyWH+iMjJvh3LaG1/eIGSMVgyJlKPpbo4ngyQeVCmcGEV6h3oZQKLzQXs5eY"
    "0aKZNLjciypf/QIIxqZGv3X8Dz596CmpSqWUXTYeatKI4msXnvzEskcacQ8o8fuvX/q5939xg7U27O8ZaNQblbgihCBNhgdquZ8J"
    "o4Im5h9VmoWOKZUawmjGG1WWKPr5I7/SqPcgiNHx0ctvubh7GY0AcPAeRzZqDc2mfo7kbStZngBcREoK9lNCRJPNyRMO+ty8wfkZ"
    "HBzGejNLegxjpUUe1R6QSz3uVl4To9Z67uC8jx74yanmRHYCOVOaTD9LOcrBhIREXrw6eQr5cILIUBR0jYZLYyGvSeSh8n6gl3+5"
    "FqoRpmbs2r2SPwzOZUn8KrW0Pl0OJlF2TOQ2Z0glidOh+JXZLYEprCebkz/45E8/e8SX0jSVgnHdPLYTKK2iKPrZledcduvFc+fM"
    "G20Nv2/Xw3bebLdO0jHR2cbYOX/vKEzEnuH/mDwmIQ0hNJJSCpmFjgkhhTCsIYHmdyKOIiLqrw/2VvvbnXZd9lx35zVjU6NSRqHr"
    "LYEBJbfd+O1brL91q9NEIawhPLg0BsvkouKEka9fu4SFkNOtyW03fvtxB3xca80VcDxznY0nw18CymOgAUgDy7ewikdN+kPv+fjG"
    "a2/WbDcF2jlaHhVuXy4xZSxZASPxZI9gpIoB7E7+jYFQYiZWGNz5oL9377DpXVZVQ5kcoNi9Aw9ZCm3PM69VDDxKWZ2vLe+LSIPG"
    "cI7tu5dapBQQECKMkiRJVOenJ//2Ywd+Kjv7eXICkySbGvUf91x5zsWnD9YHVaogpq023VppBYDCDD2Z6BqxpG/3o6XK+CT530qV"
    "QsRlr76wemSVFLIa1V5YsfSOR24x42e+rsx3U0pVoso7ttu3lbSEuxkxdz7LRRcs+iKfKdgzhAKXTgTspJ2Tjz61p9abKd8pjHZk"
    "tazn2ZoNj8uqn8JUyIL0iFrrwd6hTxz2uen2tJDSqZeoBADP6bicTYA8RoECjWM4eQVPb+mdlOTGeeTh9OQzDgMCdiCpKJ+LOcYN"
    "5T8/TPoDBxz4ktSSBcOBMAgYUVggOloediTj6fZUHFV++5W/HP3OD5okbV7X8VdlVv99T991ys8/WZU1UqA1QQovvvKCFJJIa9La"
    "tpW5MZ/9pc0XZP+ltcnkM/+nldZK5f9tsvo06TiKNalzLz0zaaUCECVCDNffdw2DLKmYc/zO7d5drVRNWGUOxhGGpz6iX484uzLW"
    "DkgZjU+PvfPt+71vl0ONQu2tTHTQTwkV5aOv7lNNRKFJH73Ph9+2yfZTzcl8MFwywrYDMgzUiSUn9gyWalj4xsREXcDjxShwXycs"
    "K/Ftnq2ZonlFp69wJ89zwpaq3Ly/qytRKaLpBhgE3vCfyOOhIGIsKqPjI+st2ODS7123/47vM4inJ6RkW9KAns+8/NTHz/pgknRk"
    "FGmltVaNuO+KG//29IuPVytVZg9h5lqC/8oqH/YrMpxoGUkZZVoZGUXZbySiePz5R4/7wdG3P3FTX6MvTZUm3VPvuffpu1aNviGl"
    "5GkStpUBgC032GbduYs6Sdvb/oxnlY2ovdvAse25QbrWqhJXvvKBbwkhiPsR+N0Zd/QvfiiRTd0KyJpEPiM5r6MEotKqXq1/4ahT"
    "jz/zmAYnc9lklYDuR+SqaCwDOt+EPIC+Mt5nsNm+GFk4CXrqTytEZ6udq3eweBZ7dvNWgUOMEO0i2TxWC/z/WvvyeLuKKt21au9z"
    "zp3vzQQkhMEwIyAgCIogCiKTwIMAKipKKygKKiqgSKAFpLsFUXB4igNOgAMgYsNP6QbFAdFGJQqNTIFACEnuPN9zdtXqP/beVWtV"
    "1T732u9F/yGQm33OrmGtb30Dp1BbrVOeYePad17P+fSvJEkQcGhsy+H7v/GGC76xzaLlmW4xXRVISzzKAZb1W55911WnDo8Ndnd1"
    "t3SGBMZQrV6bak2+48qTP/6WNW981bGL+5bkD5n3vq2sWXa6rZZuNVvN2dbsXHMu062Z1szs3Mz03Mzs3Mzs3MxMc3p2bmZ6bmqm"
    "OdvKmrPN6edfWv/I3/88NTs5MNCn8+PcQK1W3zj8wh8ff/DYg08sqMiSxqy17uvq32fV/s9seKpR78gje2VV4qIe7aXgZSrl7yhN"
    "ksHRwfeceO4rd32V1jpRCiEHAISQFILD1hOXpogsDE8sSvITgUpiQqISrfUxB5/4xgOO/Y//uru/ZyDTWakHJv53kHULxtIQQWar"
    "ICds8DY8fyiF3IiBc5jBI9AxTTtIhp1oKVzZX9rZojeW9akTDrFkLnLEiIxgU/GIUcQRXRBBcfKKPoYxatGlqRGmSTI7N6tN9qHT"
    "L7zojMtrSS3TWaLSKqadMSZJko1DG868cvX6l9b19vRlWYaAoIrEgI5a5+Dolg/dcPYOK162cul2CpJmq9lqNVumOddqNlvNTLe0"
    "yVpZK/9/pjMoemVtwJSecFDAHEkuE6EE0q7OriWdi7MsSzAxaAhAodKZfuAv9x178IkhLc+u9FfucdDtD9xaBJMURPA8NtYSTske"
    "p4Q2dt5Nf1Cp2bnZHZavuuC0T3KjBr/O8eR7kqdUbIBoU+lVvaIYcVMm/PgZl/72r7/UZEApzP0l0Y4B8u2APM0PmdGJHGwF1RdZ"
    "rhmy/CD3E+2fsqGntj7xht9sDIKMUEGIPMzVfTZPuORJhIn7kDkLAi++jY2mqLwXyfN5ReJWUERKqSRJRsaHVi7d4cpzrj3m4BPI"
    "kOX5hIvJct02jb70zitXP77u0b6e/izLbBmQH0fa6Fqt3mh0bBp8cf2LzwKCQpVbHihUiEU7ao3c6vV6USnapFdm61J8QkOGCBAm"
    "mmOkgYzp7enL+4Ya1R/+2+/nmnONesOjBtmL6xU77ddZ7yx4qeWQiMSdiOhNyDhUg6BQzTSnLzzj0q0XbZPpLJdcemYeNgXXfjzf"
    "WCCSD0CCXe/JT3lHr5TSWu+z035nHvPe0fHhFJPA20AU6gVe56eCEwRWYGSp/DY2QmbYEYtlk2nCHozkAaVIgrmN5E+eHAZT1QvJ"
    "wr34IQJ7pVLLz3YVekbG7ge7F52otNlsjowPn3DoqXdd88tjDj4hyzJCsmhJuPqNMWmSbBzacMblJ/3t6Uf6ewbK1W+nywWXjYi0"
    "NrWk3tPZ29PR29XR09XobNQ7a7VGrVZLkjRHQnNCI5UMUU35lEBnWmdGt8qZQJZpDRoUTM1Mr37tGd+95LZ3HXf2xNQEakWa6mnj"
    "+cHnNw5tiLFACjLvy7bZacnAspZuyUsWSxaNJH4yEV+OoSQqHZ8aP2zfI049/G3aaKUUIBujFFmRwr6Bx5V7ezKNlUkVDbS8VRCR"
    "yJx36oV3//6ulwY31BuNfKiMXogqNzL2KZiIIHI+PJPReP0GAoNnxznxI6D4LvnCc6mDNs8vuPf8RwvvRhewBMy4y8+g8fnewO4W"
    "IW5M0hSIRsaHtl283UXnXvbWI84EgFaWpUkS+cYk5vPClufOvOrUx575a3/vQKvVsnUroZjgkK2p88cx7OWU96O9TeVNhSwJmFyk"
    "Gqimae6w3cs+c87nuzq63vDKo/Us3fiTLy1bukzVcNKMPb3xiR2XryIy6HsTKSJa0rd0xZKVm4c3pR01y6W0krcyypyioxwE1EZ3"
    "dnStOesziUq0zhBUWfRSqNEGeS2w+POizVbzoEbM0Nh32kE0xizuXXLx2y+bbc4oUBjTQBE7kSpMuILdScwt3qfyleFIzv+HMDxe"
    "uZEjUjzBwy/tBFWKLIuBlzMurU0MI6VlMAlKbJn3g8LKJf8hKsFkamZizsy+89j33n3dA2894sxM67ywgaAG8Fb/0y8+8ZbLTnj8"
    "2UcHegdaWcsKipnoTzikRIS1yNgynncA+xTI+iZkWLJKVUs3W61WK2td9f5rj3n18SMTQ42ORqabz2x8MgoxI4IhkyTpyqXbt7Km"
    "ncWQ9L/xyT/s5E2TdGxq5H0nnb/Py/bLdKZUYrl0zOySSDLLwrPfvllVaWbPFU9+TmnxHyhUWuuTDj3tqIOOH50cSdLUMVjLhcmU"
    "7ShdeZg5iZUIIbv1UNwkxCdGUgBMwiiouErZ6yuORGbZLTNfief48HfPMkG4hX2Z6sVsKz1TCQFgEzgT3MLtASlJ02Zrbmx89MDd"
    "X/O9S35y7blfXrFkZZZlCSph5C12YK5xydIkXfv0n9+y5s3Pbnymt7uvmbXk6I5YyhmKcAUXdVilX0Oe8YSeA1SZB05gGvWOp579"
    "+x2/urVWq+V8hC9e+K0D9nrV2OxwaurPrF9XxQjPndO3WbLcgEFO8nVcMucHz00hcvRlcmZyn51e+cGTL9C55ZacXzJUVWbvhAUI"
    "i9KBWKeCQVHseeKiQ8kR17z76oHeRa2shS50FJHEXe+8QNwVgmIQgVG/rfjqApec6Fa1i7+wdzlrYNAb03mgMHFbckCR/SHiYqTy"
    "Rl56HmUoRhIigESl07OTO67Y8eufvPn2K39+6N6vzzKd4zlxl7rS2coYnSa1B9be99bLT3hpaGNPZ09x9ruRn02DBW/SIJSBLrnN"
    "b3cI+OYWN5DzSAUkQx21ji/efu2m4Y1pmma6tahv8Y0fv2XnrffY8OL01PQEVFju5T9tcd8SSxxD5plbut4SuuvPGSLn6twrz76m"
    "u6OnhO8DDhlRyEqkcNpryuEEv8wZ1FBGr5VhjIyaJqAOpZTW2S4rd7vg9IvHp0bT0jlChAW7PG0SsR+lKs8nkIdRAz4TQU5XyQKV"
    "yNw2UXjyoWdnVzh9+e4oJP7nFA62+UZib0RcM+FpCkRB8C0olUxOT+627Z4/+vQ9bz7k5JwjkMtEKBarURqXGSBIkvTm//zWmVeu"
    "npqd7Orssn0kr3etRhTtW3TxMoJ0yiIyrWeSi54i5KJQlN5SZIzpqHeu3/z8v33/cgRUqIzR2y7b/uY1d332Q1e/f/X5AKAUVpDi"
    "IC3HGlh23j6nV2RsEQGkSTo8NnzOiee95uWHiqE4N9l3jgrCM8LN6skrtEjxPSZMB5CniJaDLmdZ7L7CHBE66/gPvG7fI8emRvP8"
    "ZOIWnoX2gbgBEScP8ekwv40JhL0CIlY2KmTHJ+Tk5TKYjVjF5bPuGDOmKJ1RYmrkskqBWTPbLVAcZhQbo1trVwJENKAHBvq+9PGb"
    "tl60otlq5qsnYK3ycW/uL5Iopa7+/poLbjg3Uapeq2c6Qy98AUHoCCz9jlyNKUnGNmqj+KYISzTG2WrI+gGdOZfWWV9n380//85v"
    "1t6fpikBGDJbL9nmY2+7eI8d9i5OaIpXQlmWgWEZByL6KN+lyPN8U5VOTI8fsOerP/KWTxjmAOC0UFghPEIuM5A6ekQgUQJRQElA"
    "FgcGzmrBs/4mBIBaUrvqnOt6O/u0zpKCd2XlbeQmwijU9J6+CrwZdmml6nkNCCKXfUKUxTywfBdW08rgllK6wd0Lya5vB0rbQ1SC"
    "p0gy9Iqs77lf/lCJA6qJqfG3HnHm7tvt2cpataTGlzunqdi/JDNZopKJmbFzPveOz916dV9nP4IyxvDY8ChxhEeP2L1dGs+AyBXl"
    "nAMsFa7E3NeQGYyUPkVkSAGCgn/93j+3slb+XRtjsiwTZlMI4TaYmB7nNBNr25GfPkbGKyAqrbNGveOzH7yhu9FdmI1jhNlW7lyf"
    "AoyEwdFQrC5F4Dk3WUYq62KZlTg3QWZJRSrT2e7b73nJO68YmxpNVAI2Hgd5n2ltMAVXmzGaMfTIjmTMIPAFJwomEgANRfKIPb9F"
    "rlmQJgWyPSCee8giO3nOBR8QcPKJTYAjMsqkh+19JBGp3LW4mJWFwa9ARDng8+iza1evOebOX/94Sf9SQ9pQjCnsMHCHWxAL6bFf"
    "MHocYW4abGexZS6tNXO00TU8PF2D7uvuf+jvv/vR/d/NCwEASJLEtdJkfbXFK9wyukkhejiPZFi5B0sSNTY9ftEZl71i1X5alw4A"
    "nqcN+ilrUiIk0EvWkqISQSzclVsmd3Iw3YPm8ws+USrT2ZnHnn3a4WcMjQ/W0pqjNRAq5JLeMD0aICpFIMFpE1wNuUfY1S5TBY1/"
    "THJHUBGwzbcDEvi5T97tKCooiggJROCzBXUNmq6u7q0XLy9BSn5JibzUXNmUJumPf3XzyZcc/ehTawe6ClM3dCCal6FFHk/OGh9z"
    "IplzZc/NxMll5iDPuSRxq4oSwVkygSbT1dH9+dv+ZcvIJqUSz9aekDwTz9w8Yv3mdbVancjzvBQBQfn/0iQdGR85/pCTzjnxPK21"
    "7CsInNFQ6SzmN/covNnBP9lUJMHG4/QTcqsZXwTPNxAqQ+bq931h71X7jk+OpcJBSDjKkRU0YrWvlPv6nF8u53igmwcH/rfEEFI+"
    "e2KhJh6rP+jBfa2ZAPqp9HNF8MbkKA5Y10/ZKz5pcFjMoy0W+kttdKKSuWz20q9/7PzPn51lc90dPa2sBVwUhK6CR16msbw/F6SN"
    "6MdUe2WztWCyJEgevU1CuCuGOGQaSce6Dc986fbPKYVQFRHu+kDcPLLp2U3rOhoNQ4ZFw8oWFRAQEpXMNGe322r7q8++rrwtlNRG"
    "up61nXCTIhVA/nkV+BKzyqGw2HgU8WtWiETU1zPw5Y/eNNA90GzOJSopVg+FfOHS9y7GyI/Q60nCR9ItSpiNCISjXCLoET4q0FYh"
    "0iRm1CHOdofSOcAO+AXklHcoiEQJJNNzU89terrgGwjlMQCAJoOIaZI+8tSfTrnk6K/ecX1vZ0+SpBll6BzVECQ7XrqoleA0Yniu"
    "SZt6Hp1GHtuF5BAuSKxwlYLWWW9H/3fv+foT6x9TieIqc+nLALkq7U9P/GHT8EtpWrPcZ0Cx/vnclMhce95Xli/ZtpQ7kiyYAxp9"
    "paV62RRJcpGSlaR02MCKAzEo6aw4Iw9H2G2Hl1/34a81m808aYaZIzIvDAegYPSByWUOsEzZyKOhd8iA6KiDKbEc7oo/j1KwCeg5"
    "kiJ5ABj4KjVWZTJyimNxEZGeMz+877slDqYtA8oe/JqyG26/5pRPvelPT/xhcf8SrbXRBoNT0gkIsUiai2CuNB+9xU8izOHhnF3i"
    "lBw2nAwh6jBJaVIbnRr52k9vyDnJJOdqfCyFgD954EekjT33iQkV2FweU5WMTAxfeMalh+97ZJa1XOYNVp/T5NcU0ioHA5YDqGDk"
    "EioSOTzt+VATuwpdkdfKWkcdeNwl77pqZHw44YoZbsPrSMjEWky+GouBq3Acwph+t6ShYpkqw51yhc88ihKTIO4jLaozhmOQ+B7L"
    "fIPIihAtCzlRMBnSvZ199zz4s2/87Eu5xKSY7wAqpdIkffiJh0755Js+/Y1PgKGert6WbskTzPMVciyOkiNCPpk1ggpG5nPOFi7G"
    "HUC39LkG0/1N2mT9vQN3/f62dRufqqU1nRtQ80kUFgyO/17/t/se/nl3Z48x2gu4p9IKgghqaTo0PnjSYaedd/LHI4Rw8h1aIVrJ"
    "UMiz91W1yWWXXU7kyfw8Go5MYmP0KGc87ZzvoaBIGH3Qy1+zcfDFhx77bW9XrzZaVLpoyRUEDGCQAGtpIOTMV2RzGj/GyH9ox+bl"
    "cIkIkWdIZDnHDuI9kaUbFJ/A5+lQPM6ZzRTyT5smtXsfunt0YnjVtrv0dfUnSTLTnH7ihf/+4h2fW/P1jz23cd1A76K8RvI4eQqd"
    "phdjz1gc4AWhvVpfhMwqWX5Vfgttu0yU8Yy+1pOSJJ2YmXjmxSdPOGR1La1nWvPrymiTpmmzNffBL5y17sWnG7WO3ELLo9gUjW+a"
    "TkxN7LHjXt+4+JZGrQMJrfi8mMkyoZ7VoKDr5H1/W5m1KklRWutclRJV97nZsBcEUUVUZGw7QJhrzp225riHH/99b3d/lrXYZQRE"
    "gmLJzV0ZQk1RXrKJkeNJ+J/7O8KmKSALl+dJi/bYtgRVqzrgqQMQ5OgEGanERhvoDU7KvxGBYGxidNmirVdtt0tno2NodPC5jesm"
    "psd7u/qSNDFaF3VxSZl2ilkmbii3F7lBIWeylHVFaNhFTslWbmQm9UTPJYnL+8SwGkmmNSQqGRsfff0r33jteV/edtn23gt6YfC5"
    "C798/i//fG9vZ29mc2VQWsERqERlJuuodd521c/33GEvP/TAlWT2CGKfr2rbY/CPFoD1QrRJRG+6YAIX1SCXuR8P7FZb7oWUPL/5"
    "uZMufuOW0U0djQ6ttfNXy/+04kBbRDgGoTal8uMFj2cREisJcC2X573pTFit1BGZfZGvv3GlN8VptxC2ZyxrAAkQE0yazWZLN0FR"
    "Akm91kgw0UZ7ccLCBQy5fQ6VIhWUySTkrQbLFXayYKG9828yFDAAiiwSRjon36wKACDBZGxidOXylWceffah+xy5tH9ZotTgxJb7"
    "//SLm35+40ubN/R192c6Y8RI9oSUR3rhTHPmW5fcetQBx7nV70v0AKW5VFwxVH4LVI0UodZGMgx4PDVC+G+s2FKceViqosTlkMtV"
    "H3rsd2+97M2AmEiIQA627Or3Bs3e+YoonDjJ3+K+WxBxuiZ3N3TGbNyinMBvDDDmdOgLLsuvnxt5R3eCcYlHBapXjopMTldlrRH7"
    "+eKBEITxoGBwiSd2h5fHqo2cH+wkZodU3l0HXtckpW0A/ERJVNLUczOt6a60pzPtSVI1rScnpyc66931tKa15nx/EuEwmKhkZHzo"
    "M++77r1v/mAra1nKkK1yxGIQB76vfvXGAabCLwKNNvPgnvJgRe6/jywrOwztAAVA2mRpUrvjgR+ce+2Zvd39uR09hkE3IDn0hByO"
    "8dOV3NwkCqFWoMIha4ZZh1p7CIpfnBVCBqdMYBYSxBzsrADNysSsTx754BFvwNmCDIU78cMOZKmqrEkBj0hjeTwL+oxt/g3/dO46"
    "KqTeSiltMp1pQEzSJMFEkyHmf+wf60RpWhsa3XLe6o9efta/ucAbvt9EbxPEX0cyuMhPaAx2sqoIHKgoAq3jJQq8iLxcW3aKJipt"
    "Zdn/Oez0T77zipGxYYUpO8KJ5OQP+boRxZU/mbIkzXC+h21em3/V85uChItKmAEWHi9MLYWIMpSJ8+5cd42FVRybySL7MnOWi3Kq"
    "BQw287yrvxzilponIcKJSYECz7xKhBFE+gKxmQgnLxGgAZMZDYRJmuYujlkewUqRcOEy8KY+NDq4+g1vu+ysf9G6oLshd9jEyKi2"
    "5C2z47SCXEAx12nOQ/TNOcg3ZIGIPFlOOqgCnE2U0lqfd8rHzznx/C2jm9MkxbbeQ8gOfoIYCRQxHML6y3x+41d/sgnSlT62dWRV"
    "JkBF9J8TPSYH3zkl5Z0RqZEl/7HfRznlRqj04wWmC5dDwsiUP/KBZe0FEadMay1PnnyDEOSoHMQghiIos7vuiKCW1kbGh1633xHX"
    "nfdVe15wnJ2w8o16hEcS8wmkYFF5bCsVOfCYpTGil1rCWiAApxngiG+w/XLcUxt9xXuvecdRZw2Obim0YyR5cp5rVNG9sKAeCugi"
    "GF2gWHGE+fcvxeZ9pTLN84tm1u/gkx1Kwi1i5P5BR3JhIpXSlItXuACKxAjEz5V1JjPipaEvgxV+1GUae3U0FEaZWMSpUy5JkTC4"
    "KrCoWflW8zBmF/1OzrujkEelSW1sYnT/PQ688eLvdzW68mxJ368SkCKLS4YAuWmaeDK/cJA/RUXiMVBoxilKF0CU5UaltTFjuSAR"
    "XXv+V0449OThscFaUiMGX1sCLIEILENqV535IXYVkdweN99p5rGK6BHlR5GwGBNepegHTormtSgVuO1n6ZRClivEPC3IPmcREAcY"
    "Mq+8F1ohq8M4zwXk47mJu8XFmDNnIbFE6ZCKTKBaTIpt3I208Sc5SnVLxQCkSTo+MbrTyl1u+uQPl/Qu1QHXX0AAMoILbUobehGM"
    "GJvzexF1lguEzLwzyIokTtZvs8hLD2QvKggFioN5S/ClC246fL8jB0cH62mtBKo5LZYscYpL83gMMYKwBI2UsxgUKbzAoBKQdNeG"
    "HUcX5uXykwb0Y2RjUYvCkWRAWL9GcnB1vq7LpgZ5VcOeBfm2RxcSQiiJscRAG+c/LfkaFA5QhSIizO0seyMbYl8SAPkRzvYOHyFz"
    "WBBtAJqYnFrIk6CW1iZnJrbdauV3L719+eJtC50XVd1T6HZXmbiFAWLHcDkKhuYYVtrKDch4kgtxnrF1LxEzNkElR97+RrZJ2RQq"
    "Y3Rno+sbF9962CtePzw+XEvr/GMUy8hX7yIHtsvJDr+io67/hMxegkQCE7npbL7wjSlQWEeg5aU4WuKeezA+ArflirDXI75by2Qa"
    "V0yRXX4IvqiJuRFzqM0bQnLEgVdpLtXGs7xAZJYAbnkyAR0yI3DrtleKvP24elcJ+zRmJ8Iil4vqbMYxF0ZPTU9us2T5zf9856rl"
    "O2dZAftgJLOEk9VZ14Q+sZVEGYjEQZpoEktRdfJDhbxpOPq3ZXDX8lfn+Y+6Oqo84HLZRF93/02f+uHBe712ZHyoqIVYw8+zE8jl"
    "KDqiiGhHMXBzAymmdsnMQMgPr2I5I/PaRr9LK78mREeAEEhIecgHW56xlMgzzypH9lZnVnxEJPAVrUSizK8YD3oxaU5r4BS/Tt7K"
    "jjIi4bKN0m6SKPKbYtTjxC5oVybLeymsfvjpU1AXakltanpycd/S7625Y7eVe9rIDy4+psAQVy5T4nF6wuwG+e3G5qshkx9K4zGS"
    "ZDSSkb1uL2BYKUGsAiVOEQVp71HugYFvX/KjV+3xmuGxoXpSZ7R95LuQbXKUvWyp1BOLxrOe45GgBarI2AUFp4TYIUmejLgiMzB+"
    "QHjFBbkAS3cOoW0H/Lmm7Lo9VDeSukrs8isBj2JnsTwSu/iwXFmEzCuQuXgw4Yzt6knwpgRwzM4WEnC4LTXJCsGR96YEaVqbnJ5Y"
    "vmzbH175sz132JuPe4krNfk1GBudivVGQXckzW65rbRsgq2LkvC2ZFJw7yMCk/DF41dKPXsVdEXFHljUu/g7l9726r0OHR4fSpOa"
    "u8GcAwXG23Bw7mdkUSuMxQ/HAFS0swW0NT9a0wqoKqnKgWIUNS+rI28vomN2FGAQK7UgEmIpM7HAyQ6i8/N2kwq3KQsmHxMn2XVG"
    "sXa/ZJeCJzkJqwgRmwAKUU6w0XMlQSBIk9r45NiOK1b94Mq79th+ryLoEmQqTnv0mihUJ0BVooeDNBFC0w0g5SdpRDQAxA88u1wk"
    "5kcQd8QX5Tp/L3YPfG/NHYfvf+TQ6JZU1axTKkrtooQyKBKl6SbrTA/GCVsC6SSMUEnABb3Y8ZWECMUUgKiSrYTOO4GYPMCVqwQR"
    "YQGI2Sq6S0l8e+I3PXyrSGuT/GcPw3aCr1hJJTIKgA8a0f1hoopTIozy4gncAFCr1UcnhvbYcc8ffvrfd16+W5ZlKkm4ZiOO1UoM"
    "E8MsWqLwQiZrswi8f/AVBYpRI32tLiFGThp+6DPVdximQ2EjLFZYsQd6u/tuuuRHJ7z2lMHRwTStOWiEB6nZkswKAKP4pouY5/VH"
    "cNiKCD1X/ksxtS+J9uxUotkYrqEq7n5e5GCZ+VxeNRgEkvG+gxzfi+SQkVi4npjFO/OHiMcvlfcmeaMFDGA7VgtxyUcpm4wvC+7e"
    "LIDqsimopbXhscGDXn7IDz79s5XLttdaJ0lSGvCUVyhGHGyFkZ+QLUNsMsBFdggc/IjyY7TWEWox8/onj10eYV/65LvoS0WIJ6Ma"
    "MolKtMkuuP79N9970+L+pYYMkSGySRS2CnWE08jZS+54R4hr5RB5pnE8oRUqSUbVBFs5ECz9XbmpuywGWJcintZ7TEkzQ6ziyTrf"
    "WqKgjw6eV5IpwvQzx0eyubQoq9nSNzVCmOQTXtfIISYqGRrdcuwhJ33pI9/s7exzUaoECxvbe5wwYQMWuhvatyiZAyFPDJW/jdwl"
    "WLp5sMaoTHqQRCYCbwLv0DGUVGBvrlEALMoYozD5wodv/NBpF46MDxWUGDH7J2cVgxWXJFr3YutvYxEI4s8P3HnX84VlGoAKG67c"
    "L4a4U5ryrAgY5YcfLoTOf4H4VJT11mQLKOu95SAA4n4tiBGDDjt7wMrJIQUlCkQ2IAOosHTwYydyWCqj29bEoF0gpRABhkYH3338"
    "Od+8+FZ/9TtvEa8L9xBtxqLnXxb41skWDiOG5kauKftOjDGRhHbxNfnNNZERzoLelADKAxDR2RULf/SQ8pyfc0ap5MaffvHyb15c"
    "r9VL8zN7JGIJEmHUz7RIq0B3AgABKp+8KinGVPIrwT/4i+NUfEj2UUlGygQsICojVUKOGjsDyuqauBEz2RgldEHbJHTihCiAcX9W"
    "6h6RcVzRdU8kGRMYcGBj4S5ObUrgwkRsYV7wXaXxeqKSLGtNz01f/I7LP7z6onyxFdAjgUS+UfRZWF1WgGR3e6EYnnDG+e6iUPzY"
    "o7k6KV5MvixHilVHjA5ZbX1eijD8mTD5ZHfMkwzTJL37d3d+5PpzpmenurtK81efVROrbezNa4N3uBSEuwX6b12qVkqwjL9O1mij"
    "IBfyHxZpttnX7tYrSiILRF6PC3f2KwyxNC3mYvdwWVPZQCHn7wYxzSZ/AdK+zHWcTkaDbK7IjxJkyK3z7kjTdHpmsqPR9a/v/8LJ"
    "h71FG+0og5E1LeVqnkINg41IbsqM7co8d3DKlrhsTcrYZPGCKVSYCTEA/2tL2xsxZI3RUcBT7sQ3TR71/Nen/vKBz7377+sfW9S3"
    "ONOZzU+AoAdwMQ62FyenLGNj6Og3VVJdRGQweh6MRPPy4yNHlFsTyI4ifk2Ct8al5MhJwEqPPccM9yKfvBcP1sqYkEVfUPwjYMQ7"
    "hIPgXkIGryJDbqA74hKVjIwP77r9bjd89Jv77XRATnOwUmoSSSiIAZAvafhlWhaBg66Z0kIc8cQcgkNtFQV8CK2NF64kN4CiPPTY"
    "b4OE6Swi/kPdDFkyGIS2Cphb4A+Obbng+nPu/t1PF/cvoQJBIRtZJVYN8S6ZbTdHPSNeo/r0YfI8aQXlGIlVoAhs1tNuL5AYKSJr"
    "Acj/GiXJsH3jTXFbByZjZvQg8eGQ+xPH9V1O2WNJlGWEnet4BamAC3/cjClXsA+PDx110HGf/9BXtxrYJsuKQa9r5ctQo/k/dzjn"
    "IoiLljkShITtpTDAS6Cw0PW+Gd5n+Zk5wQ+PzqCi/xjXMAOUIShE9JnvrPnij6/paHQ2Go0sa7l0B+LKPzG2QRDqrGADlMIi/nl8"
    "RSHyfAsPsgx5F47ghRieopHaXRRtCCj+KhDRxcGrLpt0LLFWLuEHEOoo+eTE5aEek8I3/JPXNNkyKBgNgRTop0ltdnam2Zo7//SL"
    "LnzbpblPlFW1t1PzhFM24RdP4f1NEHHKcQVaBWAY2QARKZmod2PthSeKjztEeOvb3UsWN2xzBuS2GUqpu35z26e+9tEtY5sHehdl"
    "umXcjJ3tWuFFVdC5ra7QVgMs/FU+IoZ6Wr9AhnbAqPQ2KNt2IvDs24p5NUv2xAXFJAcPIymcvkiB2N2GPobFqb4S7SNHiUN2RnkT"
    "TYrXT6hUgsnY5MjKrXb49Hs+e+zBJ+TOLj7cWc1CF11lfOG1wa4XCFmLATuiLYEguAa9mocXVeQZ1YKXbs+x2MLM2gb3hFu2+kYj"
    "orwtXvfiUx/74gd+vfb+RX2LAUFr7X8WyREtsaOQEusbHJR3O3hHMS+g47p4b3VSlMyDbEhAPOKgnCh5Otd4RivTgzMijO12/Yo+"
    "9HkgL4mT20ig+JoiJxL3P7FFFTEbAQBIkrSVNSdnJk445JQr3nvNiiUrbYaFo9OhvwAiRXdoDhXaL82L2bTfaRDCoOzPiwXqmQAF"
    "531IoPZwNPFp3ZEs33bYasuyXGudpmkra133g6u/csd1mnR3R0+mW/x+j5SnkRvernO54kOelQdRoY24cFWHR1b7R+ZmAaYJXm5y"
    "1SiOSvqW5w0SHJPBqY9lkEMY++hm9BUmreKLJdEDEYBCTFQ6Njna3z3wiTMvf9fR5wCAkLQHHJJQS0leoHSbDVDZakbsRNhuD0a4"
    "Fkow2oAoJP3pphhzMo9LZuQI0lCEpC4tdg2RPNUwKJvlr/wmRcTfrv3VpTd+9G/r1g70LsplliATGFDWsh4AyjFHgWpH2CUe2mib"
    "BwphE+Xl3bdtb6rvd0LJTnOTb6uGWeAv+9Yi7Tmx8GbhZY0BcZJ/le6Dk6sF0jRttpqTMxOH7//GK95zze7b72mxTu4TNc+1Hx9R"
    "e2N1YiooH1Vn2xSrwFXpXIaxJljiAeIIJ4Zv2vcBAdgRa0d40rs/WKB51eu2moY8lypN04np8Wtuuerb93zNGN3d1ZvpFlchkDj+"
    "kdzCYsC2Nx5gh73xDMhY2SQ0yXwiK+F1oso9FaVwILmEG/9zY6iFZVC1CI3iH4eBAChd4IK371aN/Y4KgRmRP99gXA6C/Dwanxrd"
    "evE2Hz79oncf+34EzNEem/gbnptQTZ6QtY8EfWzsATOjIWHwwloWB5mTmPQxbxoxG64chMWxKPDbAJ/8U3HyxYAmv9xCQYSJLp4i"
    "ShHgwUd/fcVNn/yvxx/q6+pP07QwG4MITOkMeewr5h6jLKYH+UDBMxbm8co2aa6qOW7XMYtjHmO5VP7PQNYgUhwpd19bwd1ylyqx"
    "78AdoO6TAnNptT62SKXEk4i7CxUPlKp0amaylbVOOOyUS955xY7brDLGAASxjdxoQnaNfiw0VKBD0f5QrjHfYAm9zQo+miP9YKtQ"
    "oAUBroDz4xZhie9tAGG7WLUBRGFXdMZzrdn/+5Prv3L7dWNTI/09iwhIk/YpXSVFgkM5pRtmxH6UIvIyB7aUV4GAXWm+un+eNg4B"
    "IgR3QZaQOSrukAvHgbZOkTvYyVeFW2q4Fj3kXczsABATlWZZc2JqYtftdr/w7WtOfO3qsOKPtGPzlz8+B4IiNUy1MSfxcQhWGlQS"
    "k/Dbr8G7ATxiIwh3eKgkhMo5XDDgl7xQsEL4BfSNkk7sEoTKruCJ9Y9dc8tVdz90p0LV3dmjjTZGM7s/itrgcKyUWfRIxXekOS71"
    "i67iYFo74vMmrFzt4NkUshLAVj5CFsphSrEKZJlhSyLiPApvVgC+EbdPphEtE+tIEpUYovHJkcV9y/7pzeeefeJ5/V39ORynlBJA"
    "cIWzVvV0z67IIAqU2JgXYj0DePN9dtPJ68J/KRY7MYbi1tALhVdJWCe3uQ346c6vKopxisJrUWyJoh/NrwIAuPePd3/+h1c//MQf"
    "O+qdHfUOY4wBw/sCmV/quCTWCgqjXAa/dMGS4iazse1cW/pnSXYdkCFElLCjnboVLx9RJiuTU+tzsFIQDD1gkQWuWEYweiky5Foj"
    "DGz33I2UE9pQEcDE1Fij1jjx0NXnr75o55W75gd/zm6YZ4XLzcCmZi7fKVoZBbwBEA4xxD1MozaykavHnryCCxS3VY4u66rJF4mB"
    "lFexcGE1obg0HL2JSkY7GJtAJN8NqfIq4q/NFAnSSbM19/17v/W1O69/esOTPV29jXpHZrQxmqvt2S4k7pYrtbuA1SMwBUhe9255"
    "NiH4GOB9vGOm2Gg5YnssWUnoTzdiwIFV8BL5HNVA00jemNA9GymVIODk9DiCOuKVbzrv1I8duPur/aXP73k24abAKTnCbfAtJtBZ"
    "vJQ+F+Sb3gbDKAp0e22WKEapEHktBAYh4nXO1FjMB4uRFiumY+7pEPNFjYykCoghizjGLQ1rstgKM0YnKgWEkYnhb9/z1e/d+63n"
    "Nz3X3dFTr9WNMYaM3ANM5ccWnfB6C8oYZBzlsnwq2C/5ycRLFA+TqyqK5Ai2VCFKg0Jes9vZLBFUTJTazmzJsdgDP0m2tRATTIho"
    "cmZCoTr0FYeffeL5b9jvKADQOsvH8+0+V9mHhzwj8LgFFNB7Y21nFFjxR2binqmI/Ym+gigKRFEmHVTDV/Knh2AXr09RzuuDSi5K"
    "sSZfowuBTTaRNiZNUwDYPLrpO/fceMu9N72wZX1XR3ej3mEMGdKig+T+sYTu/Antl4Xvi0c+sPunhNpciDEnfjsyKMUYRt6s2ZGQ"
    "K6Y8BAGex0Y01Z03euUuygYbEROVapNNTk/U0tpr9z78n44/98gDjsnHkUW5H4f2K94axI8AMabAih6yeo95j11Zs3h4I3moWv6N"
    "GOL0VIRILyvXukeMAxmfIgs1ufO8n0Du1iOPFyzO5GCU6AceScKXKbfBppGXbv3P79zyH99+5vknG/VGV2cXAGqjSXQCDuNjOCgx"
    "v0tnFCnlJaxzC9fyPGPLEB8XRJiAmeAOEBvCzdw35YjGO/u9hUGcTO2qNqWUwqTZmpuamezu6HnDAW8669j3HbL36wDAGCLSysa0"
    "yC1b2fGLaB3idFjr9GS7FOamg7HQB39FI9N+gBeAQ/MbpLPGo7wBKjdAmx9HsiYSTbg3zeREQnBkcvbhffC+tH6EgEkCYXEZCJEJ"
    "XH88MjF82/23/uj+7/71mUcMmZ7O3iRJDWli6S4+NGjLl0KE4ANIkozalrgsNzNG0QH/aGMvkly7QmwE69wW0c5JA1JsqNERIFdx"
    "iiilyND03HQza267dLtjDnrz6W94+767HJBDbbliO7IC49pUHguFsTF0jItIMXcKguq0o2iDGlJu41Zw3lfBBTFBrSb2eqwI8evO"
    "iDAgRh6eZ3wWuYWgzeiTOPU5yqUDgGared+ffn7LL7794KMPjE2NdtY7OxqdBKDJ5LTT8ii1iCEyjjFVV5euIpH52zQfWUHw9qII"
    "huNsFxQMrFjaJCtI8my62DDeMTuUUipRzdbc1Oxkqhp7r9r3xMNOPenQ1csXr8iXPjsTYzxgaNOWxdZPlKAmPo/PFZeNaDAEqFqJ"
    "trjCdjzmBU2CBYCDlRAoS8+qSCMLwdC2GGvkFopSj7Fi0BaUTMboNKnlv/P4c4/+7He3//Q3tz/5wuOGqLvRVas1CMgUOZbyWhHv"
    "pmpZy4GMl98ToZmInyOgGu/lkfi8FOF5OkSLeLERG2nlN2eO3mQ6m5mb1kavWLby9fsfdcIhqw/Z67BaWssRnjyztR1v2UIAwKTP"
    "0ZiW6tdt1zk/WFFgsJwoUEnErH7xFVR2+SP9Eqgdgc1N2SOFuWCgozh6/CHmPFu4cvI3D9tAqMJDTTflWv78Qp+cnnxg7X33PHjn"
    "rx+5b+PgBqWSzo6uWloDBGPynGfDB+nt1r6raKVnUCiSY+yzaJkOEG1dCbjfgp+PhiLtMr/E5F+br/scw5mZm850tqh3yf67Hnjc"
    "a0464oCj8yMfALIsy2eL7RvQ6GUb+SMLowNWnOhV3bMgmodzfPLm+wD+EC3o0StugMAGgbv0+GddNYfNb5H5qomwFTn9Lo7aVhEG"
    "kSu+2jI48ss975IB4KWhF+//0y/+/cGfPvzEH4bHtyBAR72zUe9AxMKeaB5qD0UQd1XyhfwwDWGZGt2lWMra/XgGYj7dFrUHFH4X"
    "DEogAIVKoSKgZqs505who5f0L33FzvsfeeCxr9v3iF223a1Y91ojEMZGWgI8YZMZZ+Vd8UcogtY7GlUUIfDplCSMrSObKmoVQUDt"
    "w75CiZHYAAFXx6ZMx0aQIZAR/NuKjo947AMbCFlDptg8mIE2IrMX4zGpVQpP1iUjQFLmED6/+bnf/+3XDzxy31+efPi5l56dbc2k"
    "Sa2j3lGr1QpyGBjrmUkEvDEHH7V3FHu+9jFWArlEVuImtMhkJ3xm4BF1OJCD5S8FQJnO5pqzzayZJunyJSv23fnAw/c74rX7HL5q"
    "xc72IDDGKIWWwebk6oBVyIcbylbR17zjzdEQjZN6xLBsB4eQnDlFh2s+l0LY9FTNEaJOfGwDBB+dPRM7lBhMzpqEgGoA892A4d3n"
    "ikuIgPFRzIdiW46BsZzy7/Yce6PGGAJSiKrEOsanxtY+85cH1z7w0KO/fXz9o0PjW4wxtVq9UW+kaa2wAs/vETtzFUaEJXbFzTbk"
    "Mi5VciTyqb2Jjm2BeTIIs1ssE2RK/gRRS2fN1lwraylUA72Ld1qxywG7H/SavQ/db5cDtlq0jZ0YakNKoQhi8YxIAIO3D3GCN0G8"
    "xCcJUArs1VRa6LBaBgGpQmQrlx+0K0O4WkkOJe3SLa0RveGrjNrBecIJg/mKrBeN/CtDZlAlhwTlaMGhuBKeqWCqRt6PeC9C6GeK"
    "EE9ME4f6Pb/5ubVP//nhvz/016cfeXrDk0MTm2dbs4hYS2q1tJ6qNF9GROU0RZo7iwNc5htLqTA4Th4rHA3P+kFAzCv0onM0WmdG"
    "5yueiBr1xpL+rXZases+q16x/26vevmqfXbcelVSfhati/xthcof1Hgsb0AiUhWj3JBsKSh4yIzIYhNiaQIir+9I2nEb6o28guJ0"
    "siIRXZx6ARU/ZINC6TAZdXxgqe7s4WxGrE/UZRY3iBg/HsIxoXxov+ypKnIicJDPgGl3ZpArzIzRhihRiZ35A8CmkU3rNj716Lq1"
    "jz279skX/v7i5vXDE0MzczNaa5UkaZImKklU6v5URGNIEjkj6WccCBAK9YkxZLQxxujMZFpn+RfSUe9c0rd05dbb77zt7rut3H3P"
    "l+2904pdly9dwbqdfPgNgrJWwbAnlkDcrpKep6NlbnQxLMyj9C7EH0ROEEjk6HCpD8ZGCLy+qmDa5ZLIaPWPwK2IwpqsYkzmF05h"
    "8xqIkjDcwRSZkIsvM7Z/gpnDPD2TK1qqLEfLX0opvhmMMVtGN28YXP/UC0+u2/j0upee3rDl+S0jm8amxqZnp5qtucxk+WJQRWp0"
    "ooojXKFCz4eFyEDuBkxFZUVUQFEKMUnSRr3R29U/0DOw9aLlK7faYfttdthpxa47br3j8qUrlg1s7bX4OfHJVkftub3zo39OSIkx"
    "lnt8Ah4A8GKGg8x2oc283HI6ci4Z8Fli9KiNDhrlKDBkiMgNIBwgQQAyEG+i5xv3zKd7jGkXYpNUAaI73/DKVU3sxWE7UtRCyOvk"
    "/OapSP9OQgBnbHJ0dHJkZGL4pZEXt4xuHh4fGhkfnpyemJgZn5qZnG3NzDVn51pzzdZcprWdNCml0iSp1Wodtc7ORlcj7ejt7uvv"
    "Gejt6lvct3RZ/7IlA8uW9i9b3Lt0oGegu7PHfzRDmgqPDJ+c3B5opn/gPcYR/XZJ8jEDEZh/DZQmcJVe4pXuDSK9snSS43UEke+F"
    "mvPPFy6JrFq1YU3ZDgWKsjyrnLnanFfBxTfPVME38miHc1efiK6SJygKfyyzppQnCAx+GaNbOst0po0unZARFSQqSZO0ltYUJvP8"
    "BDLlX1ri+9h+LUqLwThtoz1uSHGp2j8ytJn3Em5TUMBCt9H/7he2nwQDtrGuC75PTlVfwCO3+0/aL80Fugz8f/uSGA/EFkvMPRp5"
    "YIy9JSx8bE14PKvu8NvOXSiNi3AU7j3AwmZygN9jpC/IP+d/8zYW+APam3zNP8pkqglc6Gfh/6mgXS70URZ8A1BlFfUPX3kLvGyr"
    "YNCFnFsLf7VtBETsBotoOP+hT0oEbe2DJJqMC1kuQQOz4APe7/fawnH/b7ul+gPN+6OKY8X66UroKHCiaksCrT7/8X8AnNHbEfKK"
    "IMEAAAAASUVORK5CYII="
)

# ============================================================
# SECTION 9: STREAMLIT WEBSITE
# ============================================================

def main():

    # Self-contained green DNA icon: no external PNG needed.
    green_icon = Image.open(BytesIO(base64.b64decode(GREEN_ICON_BASE64))).convert("RGBA")

    st.set_page_config(
        page_title="ORF Finder",
        page_icon=green_icon,
        layout="wide",
    )

    icon_col, title_col = st.columns([1, 11], vertical_alignment="center")
    with icon_col:
        st.image(green_icon, width=40)
    with title_col:
        st.title("ORF Finder")


    with st.expander("View Tool Features"):
        st.write(
            "Analyse DNA sequence length and GC content; identify complete and partial ORFs "
            "across six reading frames; explore ORF lengths, amino acid frequencies, "
            "intergenic distances and codon usage; and view protein translations."
        )

    st.write(
        "Find candidate open reading frames (ORFs) in DNA sequences."
    )

    with st.expander("Specs", expanded=False):
        st.info(
            "Complete ORFs begin with a selected start codon and end at an in-frame stop. "
            "Partial (missing-stop) candidates have a start but no downstream in-frame stop. "
            "These are candidates, not confirmed genes. "
            "Length is in nucleotides including the stop for complete ORFs; "
            "protein translation excludes the stop. Standard genetic code (table 1). "
            "Only A, T, G and C are accepted (ambiguous symbols such as N are rejected). "
            "Alternative start codons are optional assumptions whose biological validity "
            "depends on the organism and translation system."
        )

    # All analysis settings are chosen before the DNA is submitted.
    st.subheader("1. Codon settings")
    codon1, codon2 = st.columns(2)
    with codon1:
        start_codons = st.multiselect(
            "Start codons", AVAILABLE_START_CODONS, default=START_CODONS,
            help="ATG is the default. Optional GTG/TTG/CTG start assumptions depend "
                 "on organism and initiation context; translation uses standard table 1."
        )
    with codon2:
        stop_codons = st.multiselect("Stop codons", AVAILABLE_STOP_CODONS, default=STOP_CODONS)

    if any(codon != "ATG" for codon in start_codons):
        st.caption("Alternative starts are user-selected assumptions, not universally "
                   "valid initiators for every organism (standard code, table 1).")
    st.subheader("2. Overlapping ORFs")
    overlap_mode = st.radio(
        "Which ORFs should be reported for the same stop codon?",
        ["Longest ORF per stop (default)", "All candidates (including overlaps)"],
        help="Longest uses the earliest in-frame start before each stop; all includes later "
             "starts (warnings above 50,000 total bases or 5,000 ORFs per frame).",
    )
    st.subheader("3. Filters")
    orf_types = st.multiselect("ORF type", ["Complete", "Partial"], default=["Complete", "Partial"])
    c1, c2, c3 = st.columns(3)
    with c1:
        min_length = st.number_input("Minimum length (bases)", min_value=0, value=MIN_ORF_LENGTH, step=3)
    with c2:
        max_length = st.number_input("Maximum length (0 = none)", min_value=0, value=0, step=3)
    with c3:
        strand_label = st.selectbox("Strand", ["Both", "Forward (+)", "Reverse (-)"])
    frame_label = st.selectbox("Reading frame", ["All", "Frame 1", "Frame 2", "Frame 3"])
    sort_partial = st.checkbox("Order partial ORFs by length (longest first)", value=True)
    st.subheader("4. DNA input")
    input_method = st.radio("Choose how to provide DNA:",
                            ["Paste or type DNA", "Upload a DNA file"], horizontal=True)
    dna_text = ""
    uploaded_file = None
    if input_method == "Paste or type DNA":
        dna_text = st.text_area("DNA sequence", height=160, placeholder="A, T, G, C...")
    else:
        uploaded_file = st.file_uploader("Upload text or FASTA", type=["txt", "fasta", "fa", "fna"])
    submitted = st.button("Find ORFs", type="primary")

    if submitted:
        st.session_state.pop("orf_results", None)
        if not start_codons or not stop_codons or not orf_types:
            st.error("Select at least one start codon, stop codon and ORF type.")
        elif max_length and max_length < min_length:
            st.error("Maximum length must not be less than minimum length.")
        else:
            try:
                if input_method == "Paste or type DNA":
                    seq = validate_dna(dna_text)
                    records = [("Sequence 1", seq)] if seq is not False else []
                else:
                    records = parse_fasta(uploaded_file.getvalue().decode("utf-8")) if uploaded_file else []
            except (UnicodeError, ValueError) as exc:
                st.error(f"Could not read DNA: {exc}")
                records = []
            if not records:
                st.error("Provide DNA containing only A, T, G and C (no N or other "
                         "IUPAC ambiguity symbols), or a valid FASTA file.")
            else:
                st.session_state["orf_results"] = {
                    "records": records, "start_codons": tuple(start_codons),
                    "stop_codons": tuple(stop_codons),
                    "report_overlaps": overlap_mode.startswith("All"),
                    "overlap_mode": overlap_mode, "orf_types": list(orf_types),
                    "min_length": int(min_length), "max_length": int(max_length),
                    "strand_label": strand_label, "frame_label": frame_label,
                    "sort_partial": sort_partial, "allow_many_candidates": False,
                    "continue_large_input": False, "continue_overlap_large": False,
                }

    if "orf_results" not in st.session_state:
        st.caption("Set codons and filters, provide DNA, then press Find ORFs.")
        return

    settings = st.session_state["orf_results"]
    records = settings["records"]
    total_bases = sum(len(seq) for _, seq in records)
    if (settings["report_overlaps"] and total_bases > MAX_ALL_OVERLAP_BASES
            and not settings.get("continue_overlap_large", False)):
        st.warning(
            f"All-candidates mode exceeds the recommended {MAX_ALL_OVERLAP_BASES:,} "
            f"total-base threshold (input contains {total_bases:,}). "
            "Generating overlapping ORFs may be very slow or exhaust memory."
        )
        if st.button("Continue anyway — analyse large all-overlaps input"):
            settings["continue_overlap_large"] = True
            st.rerun()
        return
    large_records = [len(seq) for _, seq in records if len(seq) > MAX_INPUT_BASES]
    if (large_records or total_bases > MAX_INPUT_BASES) and not settings.get("continue_large_input"):
        st.warning(
            f"The recommended threshold is {MAX_INPUT_BASES:,} bases per record "
            f"and across all records; total input is {total_bases:,} bases. "
            "Analysis can be slow or exhaust memory."
        )
        if st.button("Continue anyway — analyse large input"):
            settings["continue_large_input"] = True
            st.rerun()
        return

    all_orfs = []
    try:
        with st.spinner("Scanning six reading frames..."):
            for name, sequence in records:
                record_orfs = cached_find_orfs(
                    sequence, settings["start_codons"], settings["stop_codons"],
                    settings["allow_many_candidates"], settings["report_overlaps"],
                )
                for orf in record_orfs:
                    orf["sequence_name"] = name
                all_orfs.extend(record_orfs)
    except CandidateLimitWarning as warning:
        st.warning(str(warning))
        if st.button("Continue anyway — generate more ORF candidates"):
            settings["allow_many_candidates"] = True
            st.rerun()
        return

    for number, orf in enumerate(all_orfs, start=1):
        orf["orf_id"] = f"ORF_{number}"

    strand_value = {"Both": "both", "Forward (+)": "+", "Reverse (-)": "-"}[settings["strand_label"]]
    frame_value = "all" if settings["frame_label"] == "All" else int(settings["frame_label"][-1])
    filtered_orfs = filter_orfs(
        all_orfs, min_length=settings["min_length"],
        max_length=settings["max_length"] or None,
        strand=strand_value, frame=frame_value,
    )
    filtered_orfs = [orf for orf in filtered_orfs if orf["orf_type"] in settings["orf_types"]]
    if settings["sort_partial"]:
        complete = [o for o in filtered_orfs if o["orf_type"] == "Complete"]
        partial = sorted((o for o in filtered_orfs if o["orf_type"] == "Partial"),
                         key=lambda o: o["length"], reverse=True)
        filtered_orfs = complete + partial

    st.success(f"Loaded {len(records)} sequence(s), containing {total_bases:,} bases in total.")
    st.subheader("Analysis settings used")
    st.caption(
        f"Starts: {', '.join(settings['start_codons'])} | Stops: {', '.join(settings['stop_codons'])} | "
        f"Minimum: {settings['min_length']} nt | "
        f"Maximum: {settings['max_length'] or 'none'} | "
        f"Strand: {settings['strand_label']} | Frame: {settings['frame_label']} | "
        f"ORF types: {', '.join(settings['orf_types'])} | {settings['overlap_mode']}"
    )
    st.subheader("Sequence statistics")
    overall_sequence = "".join(sequence for _, sequence in records)
    gc_col, length_col = st.columns(2)
    gc_col.metric("Overall GC content", f"{calculate_gc_content(overall_sequence):.2f}%")
    length_col.metric("Total sequence length", f"{total_bases:,} bases")
    if len(records) > 1:
        st.caption("FASTA records are analysed separately to prevent cross-record ORFs.")
    st.write(f"**Candidate ORFs detected before filtering:** {len(all_orfs)}")
    c1, c2 = st.columns(2)
    c1.metric("Complete ORFs found", sum(o["orf_type"] == "Complete" for o in all_orfs))
    c2.metric("Partial ORFs found", sum(o["orf_type"] == "Partial" for o in all_orfs))

    # --------------------------------------------------------
    # Summary statistics
    # --------------------------------------------------------

    stats = calculate_statistics(filtered_orfs)

    st.subheader("Summary statistics")

    metric1, metric2, metric3, metric4 = st.columns(4)

    metric1.metric("ORFs after filtering", stats["number_of_orfs"])

    metric2.metric(
        "Mean length",
        f'{stats["mean_length"]:.1f} bp',
    )

    metric3.metric(
        "Shortest length",
        f'{stats["shortest_length"]} bp',
    )

    metric4.metric(
        "Longest length",
        f'{stats["longest_length"]} bp',
    )

    # --------------------------------------------------------
    # Results table and CSV download
    # --------------------------------------------------------

    st.subheader("ORF results")
    st.caption("Nested ORFs and translated amino acid sequences are included in the table and CSV.")

    if filtered_orfs:

        show_sequences = st.checkbox("Show full DNA and protein sequences in table", value=False)
        table_orfs = filtered_orfs if show_sequences else [
            {k: v for k, v in o.items() if k not in ("dna_sequence", "amino_acid_sequence")}
            for o in filtered_orfs
        ]
        st.dataframe(
            table_orfs,
            use_container_width=True,
            hide_index=True,
        )

        # RECENTLY ADDED: Allow smaller metadata-only CSV downloads.
        export_sequences = st.checkbox(
            "Include full DNA and amino acid sequences in CSV", value=True
        )
        st.download_button(
            "Download results as CSV",
            data=results_to_csv(filtered_orfs, include_sequences=export_sequences),
            file_name=DEFAULT_OUTPUT_FILE,
            mime="text/csv",
        )

    else:
        st.info(
            "No ORFs match these filters. Try selecting Complete, "
            "Partial, or both; lowering the minimum length; or "
            "changing the codon, strand or frame settings."
        )

    # Graphs only: genuine Streamlit tabs. The map runs in an isolated HTML
    # frame: Plotly hover works, but mouse clicks are never sent to Streamlit.
    if filtered_orfs:
        st.subheader("ORF graphs")
        map_tab, length_tab, amino_tab, gap_tab, codon_tab = st.tabs([
            "Six-frame ORF map", "ORF lengths", "Amino acid frequencies",
            "Intergenic lengths", "Codon usage"
        ])

        with map_tab:
            sequence_names = list(dict.fromkeys(o["sequence_name"] for o in filtered_orfs))
            map_name = st.selectbox(
                "Sequence / FASTA record to map", sequence_names, key="map_record"
            )
            map_orfs = [o for o in filtered_orfs if o["sequence_name"] == map_name]
            record_length = next(len(seq) for name, seq in records if name == map_name)

            # Display limit changes only the GRAPH, never the biological results.
            MAP_DISPLAY_LIMIT = 60
            visible_orfs = sorted(
                map_orfs, key=lambda o: (-o["length"], o["start"])
            )[:MAP_DISPLAY_LIMIT]
            if len(map_orfs) > MAP_DISPLAY_LIMIT:
                st.caption(
                    f"Map shows {MAP_DISPLAY_LIMIT} of {len(map_orfs)} ORFs "
                    "to keep rendering quick; all remain in the selector."
                )
            st.caption(
                "Hover over a coloured block to see its ID, type, strand, frame, "
                "length and coordinates. Clicking the map never changes your selection."
            )
            components.html(
                make_hover_only_map_html(visible_orfs, record_length),
                height=360, scrolling=False
            )

            st.caption("Map is hover-only; use the ORF inspector below to view full details.")

        with length_tab:
            st.plotly_chart(make_length_histogram(filtered_orfs), use_container_width=True)

        with amino_tab:
            st.plotly_chart(
                make_amino_acid_chart(stats["amino_acid_counts"]),
                use_container_width=True
            )
            st.caption(
                "Amino acid frequencies count each selected ORF separately; "
                "nested or overlapping ORFs may count genomic regions more than once."
            )
        with gap_tab:
            gap_chart, gaps = make_intergenic_length_chart(filtered_orfs)
            st.plotly_chart(gap_chart, use_container_width=True)
            st.caption(
                f"{len(gaps)} uncovered gaps between candidate ORF intervals. "
                "Overlapping/nested ORFs are merged; intervals are evaluated "
                "within each FASTA record, including both strands. "
                "These are not confirmed intergenic regions without gene annotations. "
                "The graph changes with your ORF filters."
            )

        with codon_tab:
            st.plotly_chart(
                make_codon_usage_chart(stats["codon_counts"]),
                use_container_width=True
            )
            st.caption(
                "Heatmap of all 64 DNA codons: each row gives the first two "
                "bases and each column the third base. Uses existing ORF "
                "codon counts, excluding stop codons. Nested/overlapping "
                "candidates can count the same genomic sequence more than once."
            )
    else:
        st.info("No graphs: no ORFs match the current filters.")

    if filtered_orfs:
        st.subheader("Inspect an ORF")
        selected_id = st.selectbox(
            "Select an ORF to view its details",
            [orf["orf_id"] for orf in filtered_orfs],
            key="inspect_orf_id",
        )
        selected = next(orf for orf in filtered_orfs if orf["orf_id"] == selected_id)
        detail_cols = st.columns(4)
        detail_cols[0].metric("Strand", selected["strand"])
        detail_cols[1].metric("Frame", f'{selected["strand"]}{selected["frame"]}')
        detail_cols[2].metric("Length", f'{selected["length"]} nt')
        detail_cols[3].metric("Positions", f'{selected["start"]}–{selected["end"]}')
        st.write(f'**Type:** {selected["partial_type"]} | **Nested:** {selected["nested"]}')
        st.write("**DNA sequence**")
        st.code(selected["dna_sequence"], language="text")
        st.write("**Protein translation**")
        st.code(selected["amino_acid_sequence"] or "(empty)", language="text")

    # --------------------------------------------------------
    # Codon and amino-acid counts
    # --------------------------------------------------------

    with st.expander("Codon counts"):
        st.json(stats["codon_counts"])

    with st.expander("Amino acid counts"):
        st.json(stats["amino_acid_counts"])



if __name__ == "__main__":
    main()
