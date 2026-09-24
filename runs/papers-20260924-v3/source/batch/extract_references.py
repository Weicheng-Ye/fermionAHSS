#!/usr/bin/env python3
"""Extract paper table facts from downloaded arXiv TeX sources (no OCR)."""
import argparse
import hashlib
import json
from math import prod
from pathlib import Path
import re


def group(text):
    text = text.replace("$", "").replace(" ", "").strip()
    if not text:
        return None
    answer = []
    for term in text.split(r"\times"):
        match = re.fullmatch(r"\\?Z(?:_(?:\{(\d+)\}|(\d)))?(?:\^(?:\{(\d+)\}|(\d)))?", term)
        if not match:
            raise ValueError("unsupported group expression: " + term)
        order = int(match[1] or match[2] or 0)
        power = int(match[3] or match[4] or 1)
        if order != 1:
            answer.extend([order]*power)
    return answer


def extract(space_source, finite_source, destination):
    groups = []
    tables = re.findall(r"\\begin\{tabular\}\{ccccccc\|ccccccc\|ccccccc\}(.*?)\\end\{tabular\}",
                        space_source.read_text(), re.S)
    assert len(tables) == 2
    for page, table in zip((6, 7), tables):
        for line in table.splitlines():
            if not re.match(r"^\d+\s*&", line):
                continue
            cells = line.split(r"\\")[0].split("&")
            for start in range(0, len(cells), 7):
                cell = [x.strip() for x in cells[start:start+7]]
                if not cell[0]:
                    continue
                # The paper omits the final ampersand on its blank SG219/228 rows.
                if len(cell) == 6 and int(cell[0]) in (219, 228):
                    cell.append("")
                assert len(cell) == 7, cell
                number = int(cell[0])
                layers = dict(zip(("p_ip", "kitaev", "complex", "bosonic"), map(group, cell[2:6])))
                total = cell[6].replace("$", "").replace(" ", "")
                match = re.fullmatch(r"(\d+)(?:\\times\\Z(?:\^(?:\{(\d+)\}|(\d)))?)?", total) if total else None
                if total and not match:
                    raise ValueError(total)
                printed = {"finite_multiplicity": int(match[1]), "free_rank":
                           (int(match[2] or match[3] or 1) if r"\Z" in total else 0)} if match else None
                flat = sum(layers.values(), []) if all(v is not None for v in layers.values()) else None
                actual = {"finite_multiplicity": prod(x for x in flat if x), "free_rank": flat.count(0)} if flat is not None else None
                groups.append(dict(number=number, group_id=f"SG{number:03d}", source_page=page,
                                   symbol_tex=cell[1], layers=layers, raw_cells=cell,
                                   printed_total=printed, product_of_layers=actual,
                                   printed_total_matches_layers=(printed == actual if printed and actual else None)))
    assert sorted(g["number"] for g in groups) == list(range(1, 231))
    space = dict(source_url="https://arxiv.org/pdf/2512.25069", tables="I and II, PDF pages 6–7",
                 source_tex_sha256=hashlib.sha256(space_source.read_bytes()).hexdigest(),
                 twist={"s": "w1", "omega": "0"},
                 ahss_mapping={"p_ip": [1, 0], "kitaev": [2, -1], "complex": [3, -2], "bosonic": [5, -4]},
                 limitations=["SG210, SG219, SG228 have unreported cells; blanks are not zero.",
                              "Printed totals are retained separately from products of printed layers.",
                              "Only the zero-omega electronic case is tabulated.",
                              "Extensions between layers are not determined."],
                 groups=sorted(groups, key=lambda g: g["number"]))
    # Extract Table III by hline blocks: its piecewise row has internal \\\\.
    text = finite_source.read_text()
    start = text.index(r"\subsubsection{Summary of classification examples}")
    table = text[start:text.index(r"\end{tabular}", start)]
    table = "\n".join(line.split("%", 1)[0] for line in table.splitlines())
    rows = []
    for block in table.split(r"\hline"):
        block = block.strip()
        if not block.startswith("$") or block.startswith("$G_f"):
            continue
        cells = [c.strip().removesuffix(r"\\").strip() for c in block.split("&")]
        assert len(cells) == 4, cells
        number = len(rows)+1
        row = dict(row=number, fermionic_group_tex=cells[0], raw_cells=cells)
        if number not in (2, 3, 12):
            row["groups_by_dimension"] = dict(zip(("1", "2", "3"), map(group, cells[1:])))
        else:
            row["parameterized"] = True
        rows.append(row)
    assert len(rows) == 21
    finite = dict(source_url="https://arxiv.org/pdf/1811.00536", table="III, PDF page 8",
                  source_tex_sha256=hashlib.sha256(finite_source.read_bytes()).hexdigest(),
                  limitations=["Only phase counts can be compared without solving stacking extensions.",
                               "Intrinsic p=0 Kitaev and p+ip layers are omitted in 1D and 2D."], rows=rows)
    destination.mkdir(parents=True, exist_ok=True)
    for name, value in (("spacegroups", space), ("table-iii", finite)):
        (destination/(name+".json")).write_text(json.dumps(value, indent=2)+"\n")
    print("Extracted", len(groups), "space-group rows and", len(rows), "Table III rows")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("space_source", type=Path)
    p.add_argument("finite_source", type=Path)
    p.add_argument("destination", type=Path)
    args = p.parse_args()
    extract(args.space_source, args.finite_source, args.destination)
