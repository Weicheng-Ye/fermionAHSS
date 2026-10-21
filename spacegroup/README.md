# Space groups in package degree 4: the two sheets

One row per space group (1 to 230) and per physical spin convention, package
degree $k=4$ (spatial dimension 3): the four E6 layers of the degree, the
full group of `koFull`, its status and certificate, the package commit and
wall time of the run that produced it, and the comparison with the
published tables of [FSPT_AHSS](https://github.com/Phy-Ren/FSPT_AHSS).

| sheet | twist | FSPT_AHSS table |
| --- | --- | --- |
| `spin_half.csv`, `spin_half.md` | crystalline spin-1/2: $s=w_1$, $\omega=0$ | `space_groups_spin_half` |
| `spin_zero.csv`, `spin_zero.md` | crystalline spinless: $s=w_1$, $\omega=w_2+w_1^2$ | `space_groups_spinless` |

The layer columns are the E6 cells $(1,0)$ p+ip, $(2,-1)$ Kitaev, $(3,-2)$
complex fermion and $(5,-4)$ bosonic, after the incoming identifications;
in FSPT_AHSS they are the final filtration quotients, the same notion.
Groups are written as elementary divisors ($\mathbb Z_2^3\times\mathbb Z_4\times\mathbb Z^2$);
`0` is the trivial group.

## Provenance

- 458 of the 460 rows are the full campaign of this package
  (`bench/spacegroups_full/`, one fresh GAP process per case on the LPTMC
  cluster, October 2026), `source = campaign`: 429 rows on commit 2ccde64
  and 29 reruns on 6ec182a (the transport rewrite, identical results).
- The two cases that campaign could not finish were computed on the branch
  `extension-rows-reduction` (commit fa11ecf), `source = branch`:
  sg160 with the spinless twist, which needed the transport budgets raised
  to `maxSupport` 8388608 and `maxTerms` 2147483648 (48 h 57 min, 31 GB; the
  campaign stopped at the default budgets after 51 h), and sg218 with the
  spin-1/2 twist (24 h 34 min). sg194 with the spinless twist was also
  recomputed on the branch (3 h 56 min against 50 h 41 min) and agreed with
  the campaign in every field; the sheet keeps the campaign row.
- `wall seconds` and `host` are those of the recorded run; nodes differ, so
  wall times are indicative only.

The sheets are written by `tools/spacegroup_sheets.py` from the result
JSON files; `comparison_summary.json` holds the counts and the list of
differences (empty).

## Comparison with FSPT_AHSS

Reference: FSPT_AHSS commit 5e8cb5d (2026-10-15),
`results/computed_examples/tables/space_groups_spin_half.csv` and
`space_groups_spinless.csv`. Each of the five groups of a case (four layers
and the full stacking group) was normalized to sorted elementary divisors
on both sides and compared.

| convention | cases | identical in all five groups | differing |
| --- | ---: | ---: | ---: |
| spin-1/2 | 230 | 230 | 0 |
| spinless | 230 | 230 | 0 |

Every space group and both conventions agree, including the three cases
computed only on the branch (sg160 spinless $\mathbb Z_2^5\times\mathbb Z_{16}$,
sg194 spinless $\mathbb Z_2^9$, sg218 spin-1/2 $\mathbb Z_2^3\times\mathbb Z_4$).
