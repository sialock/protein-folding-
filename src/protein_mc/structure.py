"""Read protein C-alpha coordinates and export intentionally C-alpha-only PDBs."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import numpy as np
from Bio.PDB import PDBParser


@dataclass(frozen=True)
class CAProtein:
    coords: np.ndarray
    chain_id: str
    residues: tuple[tuple[str, int, str], ...]  # name, PDB sequence number, insertion code


def load_ca(path: str | Path, chain_id: str | None = None) -> CAProtein:
    file = Path(path)
    if not file.is_file():
        raise FileNotFoundError(f"PDB not found: {file}")
    structure = PDBParser(QUIET=True).get_structure("reference", str(file))
    model = next(iter(structure), None)
    if model is None:
        raise ValueError("PDB contains no model")
    chain = (model[chain_id] if chain_id is not None and chain_id in model else None)
    if chain_id is not None and chain is None:
        raise ValueError(f"Chain {chain_id!r} not present in first model")
    if chain is None:
        chain = next((c for c in model if sum(1 for r in c if r.id[0] == " " and "CA" in r) >= 5), None)
    if chain is None:
        raise ValueError("No chain with at least five standard residues containing CA atoms")
    residues = [r for r in chain if r.id[0] == " " and "CA" in r]
    if len(residues) > 300:
        raise ValueError("Demo limited to 300 C-alpha atoms to avoid excessive O(N^2) work")
    # Require no gaps in numbering for the simplified consecutive-chain move model.
    numbers = [int(r.id[1]) for r in residues]
    if any(b != a + 1 for a, b in zip(numbers, numbers[1:])):
        raise ValueError("Discontinuous residue numbers; use a complete single-chain PDB")
    coords = np.array([r["CA"].coord for r in residues], dtype=float)
    residue_info = tuple((r.resname, int(r.id[1]), str(r.id[2]).strip()) for r in residues)
    return CAProtein(coords=coords, chain_id=chain.id, residues=residue_info)


def write_ca_pdb(path: str | Path, protein: CAProtein, coords: np.ndarray, remark: str) -> None:
    xyz = np.asarray(coords, dtype=float)
    if xyz.shape != protein.coords.shape or not np.isfinite(xyz).all():
        raise ValueError("Coordinate shape mismatch / nonfinite coordinates")
    with Path(path).open("w", encoding="utf-8") as out:
        out.write(f"REMARK 900 C-ALPHA-ONLY DEMONSTRATION; {remark[:50]}\n")
        for serial, ((resname, seq, insertion), point) in enumerate(zip(protein.residues, xyz), 1):
            x, y, z = point
            out.write(
                f"ATOM  {serial:5d}  CA  {resname:3s} {protein.chain_id:1s}"
                f"{seq:4d}{(insertion or ' '):1s}   {x:8.3f}{y:8.3f}{z:8.3f}"
                f"{1.0:6.2f}{0.0:6.2f}           C\n"
            )
        out.write("TER\nEND\n")