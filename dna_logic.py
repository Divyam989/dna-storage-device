"""
dna_logic.py — DNA Storage Encoder/Decoder
==========================================
2-bit mapping: 00=A  01=T  10=G  11=C
Error correction: even-parity base appended every 8 data bases (1 byte)
Oligo chunking: PRIMER_F + up to 180 data bases + PRIMER_R  (≤200 bases total)

Usage examples at the bottom of this file.
"""

import os
import math
from dataclasses import dataclass, field
from typing import Optional

# ── Constants ────────────────────────────────────────────────────────────────

BIT_TO_BASE = {"00": "A", "01": "T", "10": "G", "11": "C"}
BASE_TO_BIT = {"A": "00", "T": "01", "G": "10", "C": "11"}
COMPLEMENT  = {"A": "T",  "T": "A",  "G": "C",  "C": "G"}

PRIMER_F = "ATCGATCG"
PRIMER_R = "CGATCGAT"
OLIGO_MAX_LEN   = 200
OLIGO_DATA_LEN  = OLIGO_MAX_LEN - len(PRIMER_F) - len(PRIMER_R)  # 184


# ── Data classes ─────────────────────────────────────────────────────────────

@dataclass
class Oligo:
    index: int
    data_bases: str          # raw data portion (no primers)
    full_sequence: str       # PRIMER_F + data + PRIMER_R
    length: int


@dataclass
class EncodeResult:
    source_bytes:    int
    binary_string:   str
    dna_raw:         str     # without parity
    dna_with_parity: str     # with parity bases inserted
    oligos:          list    # list[Oligo]
    base_counts:     dict    # {A, T, G, C}
    gc_content_pct:  float
    gc_warning:      Optional[str]
    parity_positions: list   # index positions of parity bases in dna_with_parity


@dataclass
class DecodeResult:
    dna_input:       str
    binary_string:   str
    decoded_bytes:   bytes
    decoded_text:    Optional[str]
    is_text:         bool
    parity_errors:   list    # list of block indices with parity errors
    blocks_checked:  int


# ── Core encoding ─────────────────────────────────────────────────────────────

def bytes_to_binary(data: bytes) -> str:
    """Convert bytes to a binary string."""
    return "".join(format(byte, "08b") for byte in data)


def binary_to_dna(binary: str) -> str:
    """Map every 2 bits → one DNA base. Pads with a trailing '0' if length is odd."""
    if len(binary) % 2 != 0:
        binary += "0"
    return "".join(BIT_TO_BASE[binary[i:i+2]] for i in range(0, len(binary), 2))


def add_parity(dna: str) -> tuple[str, list[int]]:
    """
    Append a parity base after every 8-base block.
    Parity rule (even parity on G/C count):
        GC count even  → parity base = A
        GC count odd   → parity base = T
    Returns the new sequence and the list of parity base positions.
    """
    result = []
    parity_positions = []
    pos = 0
    for i in range(0, len(dna), 8):
        block = dna[i:i+8]
        result.append(block)
        pos += len(block)
        gc = sum(1 for b in block if b in ("G", "C"))
        parity_base = "A" if gc % 2 == 0 else "T"
        result.append(parity_base)
        parity_positions.append(pos)
        pos += 1
    return "".join(result), parity_positions


def chunk_to_oligos(dna: str) -> list[Oligo]:
    """Split a DNA sequence into oligos with forward/reverse primers."""
    oligos = []
    for i, start in enumerate(range(0, len(dna), OLIGO_DATA_LEN), start=1):
        chunk = dna[start:start + OLIGO_DATA_LEN]
        full  = PRIMER_F + chunk + PRIMER_R
        oligos.append(Oligo(
            index=i,
            data_bases=chunk,
            full_sequence=full,
            length=len(full)
        ))
    return oligos


def base_counts(dna: str) -> dict:
    counts = {"A": 0, "T": 0, "G": 0, "C": 0}
    for b in dna:
        if b in counts:
            counts[b] += 1
    return counts


def gc_content(dna: str) -> float:
    if not dna:
        return 0.0
    counts = base_counts(dna)
    return round((counts["G"] + counts["C"]) / len(dna) * 100, 2)


def encode_bytes(data: bytes) -> EncodeResult:
    """Full pipeline: bytes → EncodeResult."""
    binary      = bytes_to_binary(data)
    dna_raw     = binary_to_dna(binary)
    dna_parity, parity_pos = add_parity(dna_raw)
    oligos      = chunk_to_oligos(dna_parity)
    counts      = base_counts(dna_raw)
    gc_pct      = gc_content(dna_raw)

    if gc_pct > 70:
        gc_warn = f"High GC content ({gc_pct}%) — above 70% makes synthesis harder and risks secondary structures."
    elif gc_pct < 30:
        gc_warn = f"Low GC content ({gc_pct}%) — below 30% reduces stability and increases synthesis error rate."
    else:
        gc_warn = None

    return EncodeResult(
        source_bytes    = len(data),
        binary_string   = binary,
        dna_raw         = dna_raw,
        dna_with_parity = dna_parity,
        oligos          = oligos,
        base_counts     = counts,
        gc_content_pct  = gc_pct,
        gc_warning      = gc_warn,
        parity_positions= parity_pos,
    )


def encode_text(text: str, encoding: str = "utf-8") -> EncodeResult:
    """Encode a text string."""
    return encode_bytes(text.encode(encoding))


def encode_file(filepath: str) -> EncodeResult:
    """Encode any file from disk."""
    with open(filepath, "rb") as f:
        data = f.read()
    return encode_bytes(data)


# ── Core decoding ─────────────────────────────────────────────────────────────

def dna_to_binary(dna: str) -> str:
    """Map each DNA base → 2 bits."""
    return "".join(BASE_TO_BIT.get(b, "00") for b in dna)


def binary_to_bytes(binary: str) -> bytes:
    """Pack a binary string into bytes (truncates incomplete trailing bits)."""
    result = bytearray()
    for i in range(0, len(binary) - 7, 8):
        result.append(int(binary[i:i+8], 2))
    return bytes(result)


def strip_parity(dna: str) -> tuple[str, list[int], int]:
    """
    Remove parity bases (every 9th position) and verify them.
    Returns (data_dna, error_block_indices, total_blocks_checked).
    """
    data_bases = []
    errors     = []
    total      = 0
    for i in range(0, len(dna), 9):
        block = dna[i:i+8]
        pb    = dna[i+8] if i+8 < len(dna) else None
        data_bases.append(block)
        if pb is not None:
            total += 1
            gc = sum(1 for b in block if b in ("G", "C"))
            expected = "A" if gc % 2 == 0 else "T"
            if pb != expected:
                errors.append(total)
    return "".join(data_bases), errors, total


def decode_dna(dna: str, has_parity: bool = True) -> DecodeResult:
    """Full pipeline: DNA string → DecodeResult."""
    dna = dna.upper().strip()
    dna = "".join(c for c in dna if c in "ATGC")

    parity_errors = []
    blocks_checked = 0

    if has_parity:
        data_dna, parity_errors, blocks_checked = strip_parity(dna)
    else:
        data_dna = dna

    binary        = dna_to_binary(data_dna)
    decoded_bytes = binary_to_bytes(binary)

    try:
        decoded_text = decoded_bytes.decode("utf-8")
        is_text = True
    except UnicodeDecodeError:
        decoded_text = None
        is_text = False

    return DecodeResult(
        dna_input      = dna,
        binary_string  = binary,
        decoded_bytes  = decoded_bytes,
        decoded_text   = decoded_text,
        is_text        = is_text,
        parity_errors  = parity_errors,
        blocks_checked = blocks_checked,
    )


def decode_to_file(dna: str, output_path: str, has_parity: bool = True) -> DecodeResult:
    """Decode a DNA string and write raw bytes to a file."""
    result = decode_dna(dna, has_parity)
    with open(output_path, "wb") as f:
        f.write(result.decoded_bytes)
    return result


# ── Pretty-print helpers ──────────────────────────────────────────────────────

COLORS = {
    "A": "\033[92m",   # green
    "T": "\033[91m",   # red
    "G": "\033[94m",   # blue
    "C": "\033[95m",   # purple
    "P": "\033[93m",   # yellow (parity)
    "R": "\033[0m",    # reset
}


def colorize_dna(dna: str, parity_positions: set = None, max_len: int = 120) -> str:
    """Return a terminal-colorized DNA string (truncated to max_len)."""
    out = []
    for i, base in enumerate(dna[:max_len]):
        if parity_positions and i in parity_positions:
            out.append(f"{COLORS['P']}{base}{COLORS['R']}")
        else:
            out.append(f"{COLORS[base]}{base}{COLORS['R']}")
    if len(dna) > max_len:
        out.append(f"  \033[90m…+{len(dna)-max_len} more\033[0m")
    return "".join(out)


def print_encode_result(r: EncodeResult, show_oligos: int = 3) -> None:
    """Print a nicely formatted encode summary to the terminal."""
    W = 60
    line = "─" * W
    print(f"\n\033[92m{'═'*W}")
    print(f"  DNA STORAGE ENCODER — RESULT")
    print(f"{'═'*W}\033[0m")

    print(f"\n  Source bytes   : \033[97m{r.source_bytes:,}\033[0m")
    print(f"  Binary length  : \033[97m{len(r.binary_string):,}\033[0m bits")
    print(f"  Raw bases      : \033[97m{len(r.dna_raw):,}\033[0m")
    print(f"  Bases + parity : \033[97m{len(r.dna_with_parity):,}\033[0m")
    print(f"  Oligos         : \033[97m{len(r.oligos):,}\033[0m  (≤{OLIGO_MAX_LEN} bases each)")

    gc = r.gc_content_pct
    gc_color = "\033[91m" if (gc > 70 or gc < 30) else "\033[92m"
    c = r.base_counts
    total = sum(c.values()) or 1
    print(f"\n  GC content     : {gc_color}{gc}%\033[0m")
    print(f"  Base counts    : "
          f"\033[92mA={c['A']}({c['A']*100//total}%)\033[0m  "
          f"\033[91mT={c['T']}({c['T']*100//total}%)\033[0m  "
          f"\033[94mG={c['G']}({c['G']*100//total}%)\033[0m  "
          f"\033[95mC={c['C']}({c['C']*100//total}%)\033[0m")

    if r.gc_warning:
        print(f"\n  \033[93m⚠  {r.gc_warning}\033[0m")

    print(f"\n  {line}")
    print(f"  Binary stream (first 64 bits):")
    print(f"  \033[90m{r.binary_string[:64]}{'…' if len(r.binary_string)>64 else ''}\033[0m")

    print(f"\n  DNA sequence (first 120 bases, orange = parity):")
    pset = set(r.parity_positions)
    print(f"  {colorize_dna(r.dna_with_parity, pset, 120)}")

    print(f"\n  {line}")
    print(f"  Oligos (showing first {min(show_oligos, len(r.oligos))}):")
    for o in r.oligos[:show_oligos]:
        preview = o.data_bases[:40] + ("…" if len(o.data_bases) > 40 else "")
        print(f"  \033[90m#{o.index:04d}\033[0m  "
              f"\033[93m{PRIMER_F}\033[0m"
              f"\033[97m{preview}\033[0m"
              f"\033[93m{PRIMER_R}\033[0m"
              f"  \033[90m({o.length} bases)\033[0m")
    if len(r.oligos) > show_oligos:
        print(f"  \033[90m  … and {len(r.oligos)-show_oligos} more oligos\033[0m")

    print(f"\n\033[92m{'═'*W}\033[0m\n")


def print_decode_result(r: DecodeResult) -> None:
    """Print a nicely formatted decode summary to the terminal."""
    W = 60
    print(f"\n\033[94m{'═'*W}")
    print(f"  DNA STORAGE DECODER — RESULT")
    print(f"{'═'*W}\033[0m")

    print(f"\n  Input bases    : \033[97m{len(r.dna_input):,}\033[0m")
    print(f"  Decoded bytes  : \033[97m{len(r.decoded_bytes):,}\033[0m")

    if r.blocks_checked > 0:
        if r.parity_errors:
            print(f"  Parity check   : \033[91m✗ {len(r.parity_errors)} error(s) in blocks {r.parity_errors[:10]}\033[0m")
        else:
            print(f"  Parity check   : \033[92m✓ All {r.blocks_checked} blocks passed\033[0m")

    print(f"\n  Output:")
    if r.is_text:
        preview = r.decoded_text[:300]
        if len(r.decoded_text) > 300:
            preview += f"… (+{len(r.decoded_text)-300} more chars)"
        print(f"  \033[97m{preview}\033[0m")
    else:
        hex_preview = r.decoded_bytes[:32].hex(" ")
        print(f"  \033[90m[Binary file — hex preview]\033[0m")
        print(f"  \033[90m{hex_preview}…\033[0m")

    print(f"\n\033[94m{'═'*W}\033[0m\n")


# ── CLI / demo ────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    import sys

    if len(sys.argv) < 2:
        # ── Demo mode ──────────────────────────────────────────────────────
        print("\033[90mRunning demo (no arguments provided)\033[0m")
        print("\033[90mUsage:\033[0m")
        print("  python dna_logic.py encode-text  \"Hello World\"")
        print("  python dna_logic.py encode-file  path/to/file.png")
        print("  python dna_logic.py decode        ATGCATGC…  [--no-parity]")
        print("  python dna_logic.py roundtrip     \"Hello World\"\n")

        # Demo with a sample string
        sample = "DNA data storage: encoding digital information into nucleotide sequences."
        print(f"\033[92mDemo — encoding:\033[0m \"{sample}\"\n")
        result = encode_text(sample)
        print_encode_result(result)

        print("\033[92mDemo — decoding back:\033[0m")
        dec = decode_dna(result.dna_with_parity, has_parity=True)
        print_decode_result(dec)

    else:
        cmd = sys.argv[1].lower()

        if cmd == "encode-text":
            if len(sys.argv) < 3:
                print("Usage: python dna_logic.py encode-text \"your text here\"")
                sys.exit(1)
            r = encode_text(sys.argv[2])
            print_encode_result(r)
            # Save full sequence to file
            out = "dna_output.txt"
            with open(out, "w") as f:
                f.write(r.dna_with_parity)
            print(f"  Full DNA sequence saved to \033[97m{out}\033[0m\n")

        elif cmd == "encode-file":
            if len(sys.argv) < 3:
                print("Usage: python dna_logic.py encode-file path/to/file")
                sys.exit(1)
            path = sys.argv[2]
            if not os.path.exists(path):
                print(f"File not found: {path}")
                sys.exit(1)
            r = encode_file(path)
            print_encode_result(r)
            out = os.path.splitext(path)[0] + "_dna.txt"
            with open(out, "w") as f:
                f.write(r.dna_with_parity)
            print(f"  Full DNA sequence saved to \033[97m{out}\033[0m\n")

        elif cmd == "decode":
            if len(sys.argv) < 3:
                print("Usage: python dna_logic.py decode ATGC…  [--no-parity]")
                sys.exit(1)
            dna_str = sys.argv[2]
            has_p = "--no-parity" not in sys.argv
            r = decode_dna(dna_str, has_parity=has_p)
            print_decode_result(r)

        elif cmd == "roundtrip":
            if len(sys.argv) < 3:
                print("Usage: python dna_logic.py roundtrip \"text to test\"")
                sys.exit(1)
            text = sys.argv[2]
            print(f"\033[92mRoundtrip test:\033[0m \"{text}\"\n")
            enc = encode_text(text)
            print_encode_result(enc)
            dec = decode_dna(enc.dna_with_parity, has_parity=True)
            print_decode_result(dec)
            match = dec.decoded_text == text
            status = "\033[92m✓ PASS\033[0m" if match else "\033[91m✗ FAIL\033[0m"
            print(f"  Roundtrip result: {status}\n")

        else:
            print(f"Unknown command: {cmd}")
            print("Commands: encode-text, encode-file, decode, roundtrip")
            sys.exit(1)
