Updated 2026-09-24T10:15:15.946800+00:00.

20 / 489 calculations have saved E6 tables.
Task states: {"PENDING": 128, "computed": 20, "not_submitted": 341}

These are five-row associated-graded calculations, certified_ko=false. Stacking extensions are unresolved.
Paper spatial dimension d uses p+q=d-2 (package cutoff k=d+1).
Table III omits intrinsic p=0 Kitaev/p+ip factors in 1D/2D; comparisons test phase counts only.

Table III comparison: {"not_computed": 63, "order_matches_extensions_unresolved": 27}
Space-group comparison: {"layer_discrepancy": 1, "layers_match": 6, "not_computed": 223, "twist_not_tabulated": 230}

Sources: [Table III](https://arxiv.org/pdf/1811.00536#page=8), [Tables I–II](https://arxiv.org/pdf/2512.25069#page=6).
The second space-group twist has no reference in these tables. SG210, SG219, SG228 contain unreported cells.
Paper rows whose printed total differs from their printed layers: 22, 27, 31, 37, 58, 63, 82, 83, 85, 127, 133, 142.

Measured wall time includes setup, resolution, twists, and pages; it excludes queue wait. GAP event runtime_ms is cumulative GAP CPU time per case, excluding external-worker CPU. GNU time RSS is the largest process high-water mark; sampled tree RSS sums live descendants every 2 s. Cgroup peak, when available, measures the job cgroup and includes process overhead.

Completed results and discrepancies:

- table3-C2-row1, row 1 {}, 1D: order_matches_extensions_unresolved; graded order 2, paper order 2.
- table3-C2-row1, row 1 {}, 2D: order_matches_extensions_unresolved; graded order 8, paper order 8.
- table3-C2-row1, row 1 {}, 3D: order_matches_extensions_unresolved; graded order 1, paper order 1.
- table3-C2-row1, row 3 {'k': 1}, 1D: order_matches_extensions_unresolved; graded order 2, paper order 2.
- table3-C2-row1, row 3 {'k': 1}, 2D: order_matches_extensions_unresolved; graded order 8, paper order 8.
- table3-C2-row1, row 3 {'k': 1}, 3D: order_matches_extensions_unresolved; graded order 1, paper order 1.
- table3-C2-row12, row 12 {'k': 2, 'n': 1}, 1D: order_matches_extensions_unresolved; graded order 1, paper order 1.
- table3-C2-row12, row 12 {'k': 2, 'n': 1}, 2D: order_matches_extensions_unresolved; graded order 1, paper order 1.
- table3-C2-row12, row 12 {'k': 2, 'n': 1}, 3D: order_matches_extensions_unresolved; graded order 1, paper order 1.
- table3-C2xC2-row21, row 21 {}, 1D: order_matches_extensions_unresolved; graded order 1, paper order 1.
- table3-C2xC2-row21, row 21 {}, 2D: order_matches_extensions_unresolved; graded order 2, paper order 2.
- table3-C2xC2-row21, row 21 {}, 3D: order_matches_extensions_unresolved; graded order 1, paper order 1.
- table3-C2xC2xC2-row16, row 16 {}, 1D: order_matches_extensions_unresolved; graded order 8, paper order 8.
- table3-C2xC2xC2-row16, row 16 {}, 2D: order_matches_extensions_unresolved; graded order 64, paper order 64.
- table3-C2xC2xC2-row16, row 16 {}, 3D: order_matches_extensions_unresolved; graded order 32, paper order 32.
- table3-C2xC4-row15, row 15 {}, 1D: order_matches_extensions_unresolved; graded order 2, paper order 2.
- table3-C2xC4-row15, row 15 {}, 2D: order_matches_extensions_unresolved; graded order 8, paper order 8.
- table3-C2xC4-row15, row 15 {}, 3D: order_matches_extensions_unresolved; graded order 4, paper order 4.
- table3-C2xC4xC4-row18, row 18 {}, 1D: order_matches_extensions_unresolved; graded order 64, paper order 64.
- table3-C2xC4xC4-row18, row 18 {}, 2D: order_matches_extensions_unresolved; graded order 4096, paper order 4096.
- table3-C2xC4xC4-row18, row 18 {}, 3D: order_matches_extensions_unresolved; graded order 256, paper order 256.
- table3-C3-row2, row 2 {'k': 1}, 1D: order_matches_extensions_unresolved; graded order 1, paper order 1.
- table3-C3-row2, row 2 {'k': 1}, 2D: order_matches_extensions_unresolved; graded order 3, paper order 3.
- table3-C3-row2, row 2 {'k': 1}, 3D: order_matches_extensions_unresolved; graded order 1, paper order 1.
- table3-C4-row3, row 3 {'k': 2}, 1D: order_matches_extensions_unresolved; graded order 2, paper order 2.
- table3-C4-row3, row 3 {'k': 2}, 2D: order_matches_extensions_unresolved; graded order 16, paper order 16.
- table3-C4-row3, row 3 {'k': 2}, 3D: order_matches_extensions_unresolved; graded order 1, paper order 1.
- SG001-w2_plus_w1_squared: twist_not_tabulated; p_ip=[0, 0, 0] (paper None, twist_not_tabulated); kitaev=[2, 2, 2] (paper None, twist_not_tabulated); complex=[2] (paper None, twist_not_tabulated); bosonic=[] (paper None, twist_not_tabulated)
- SG002-zero: layers_match; p_ip=[0, 0, 0] (paper [0, 0, 0], matches); kitaev=[] (paper [], matches); complex=[2, 2, 2, 2] (paper [2, 2, 2, 2], matches); bosonic=[2] (paper [2], matches)
- SG002-w2_plus_w1_squared: twist_not_tabulated; p_ip=[0, 0, 0] (paper None, twist_not_tabulated); kitaev=[] (paper None, twist_not_tabulated); complex=[2, 2, 2, 2] (paper None, twist_not_tabulated); bosonic=[2] (paper None, twist_not_tabulated)
- SG003-zero: layers_match; p_ip=[0] (paper [0], matches); kitaev=[2, 2, 2] (paper [2, 2, 2], matches); complex=[2, 2, 2, 2] (paper [2, 2, 2, 2], matches); bosonic=[2, 2, 2, 2] (paper [2, 2, 2, 2], matches)
- SG004-zero: layers_match; p_ip=[0] (paper [0], matches); kitaev=[2, 2, 2] (paper [2, 2, 2], matches); complex=[2] (paper [2], matches); bosonic=[] (paper [], matches)
- SG005-zero: layers_match; p_ip=[0] (paper [0], matches); kitaev=[2, 2] (paper [2, 2], matches); complex=[2, 2] (paper [2, 2], matches); bosonic=[2, 2] (paper [2, 2], matches)
- SG005-w2_plus_w1_squared: twist_not_tabulated; p_ip=[0] (paper None, twist_not_tabulated); kitaev=[2, 2] (paper None, twist_not_tabulated); complex=[2, 2] (paper None, twist_not_tabulated); bosonic=[2] (paper None, twist_not_tabulated)
- SG006-zero: layers_match; p_ip=[0] (paper [0], matches); kitaev=[2] (paper [2], matches); complex=[2, 2] (paper [2, 2], matches); bosonic=[2, 2] (paper [2, 2], matches)
- SG007-zero: layer_discrepancy; p_ip=[0, 2] (paper [0], differs); kitaev=[2, 2, 2] (paper [2, 2, 2], matches); complex=[2] (paper [2], matches); bosonic=[] (paper [], matches)
- SG008-w2_plus_w1_squared: twist_not_tabulated; p_ip=[0, 2] (paper None, twist_not_tabulated); kitaev=[2, 2, 2, 2] (paper None, twist_not_tabulated); complex=[2, 2, 2] (paper None, twist_not_tabulated); bosonic=[2, 2] (paper None, twist_not_tabulated)
- SG010-zero: layers_match; p_ip=[0] (paper [0], matches); kitaev=[] (paper [], matches); complex=[2, 2, 2, 2, 2] (paper [2, 2, 2, 2, 2], matches); bosonic=[2, 2, 2, 2, 2, 2] (paper [2, 2, 2, 2, 2, 2], matches)
- SG010-w2_plus_w1_squared: twist_not_tabulated; p_ip=[0] (paper None, twist_not_tabulated); kitaev=[2] (paper None, twist_not_tabulated); complex=[2, 2, 2, 2, 2, 2, 2, 2] (paper None, twist_not_tabulated); bosonic=[2, 2, 2, 2, 2, 2, 2] (paper None, twist_not_tabulated)
