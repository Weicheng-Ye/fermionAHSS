"""Group catalogue for the paper examples used by the AHSS batch runner.

Table III lists fermionic groups G_f, whereas the computation uses
G_b = G_f / C2^f.  Rows 2, 3, and 12 are infinite families of finite
groups.  This catalogue contains only their explicitly selected parameters;
it does not claim to enumerate every member of those families.

Primary sources:
  https://arxiv.org/pdf/1811.00536v3#page=8 (Table III)
  https://arxiv.org/html/1811.00536v3#S1.SS3.SSS2
  https://arxiv.org/pdf/1811.00536v3#page=64 (cyclic carry, Eq. E2)
  https://arxiv.org/pdf/1811.00536v3#page=69 (quaternion twist, Eq. E26)
  https://arxiv.org/pdf/2512.25069#page=24 (electronic space-group twists)

The space-group entries mean the full infinite discrete space groups,
including their translation subgroups, not their finite point groups.
"""

TABLE_III_SOURCE = "https://arxiv.org/pdf/1811.00536v3#page=8"
TABLE_III_HTML_SOURCE = "https://arxiv.org/html/1811.00536v3#S1.SS3.SSS2"
CYCLIC_TWIST_SOURCE = "https://arxiv.org/pdf/1811.00536v3#page=64"
QUATERNION_TWIST_SOURCE = "https://arxiv.org/pdf/1811.00536v3#page=69"
SPACE_GROUP_SOURCE = "https://arxiv.org/pdf/2512.25069"
SPACE_GROUP_TWIST_SOURCE = "https://arxiv.org/pdf/2512.25069#page=24"


def _positive_parameters(values, name):
    """Validate and deduplicate parameter selections without inventing bounds."""
    try:
        selected = tuple(values)
    except TypeError as exc:
        raise TypeError(f"{name} must be an iterable of positive integers") from exc
    if any(isinstance(value, bool) or not isinstance(value, int) or value <= 0
           for value in selected):
        raise ValueError(f"{name} must contain only positive integers")
    return tuple(sorted(set(selected)))


def _group_label(factors):
    return " x ".join(f"C{order}" for order in factors) or "1"


def build_catalogue(k_values=(1, 2, 3, 4), n_values=(1, 2, 3, 4)):
    """Return JSON-serializable entries for selected Table III groups and SGs.

    ``k_values`` selects rows 2 and 3, with quotients C_(2k+1) and C_(2k).
    ``n_values`` selects row 12, with G_f=C_(2k)^f and quotient C_k for
    k=2**n.  Both selections must contain positive integers.  Empty selections
    omit those family rows; the remaining fixed rows are always included.

    Equal bosonic quotient groups share an entry, but ``source_rows`` retains
    every applicable paper row and parameter choice.  ``factors`` are sorted
    with trivial factors removed.  Each source row additionally retains the
    quotient factors in paper order, identifying the distinguished factor of
    a nonsplit cyclic extension even after sorting.

    Source rows specify the only twists evaluated by the runner.  Their named
    cocycles are transported into the actual resolution before deduplication.
    Each space-group entry records both requested families (w1,0) and
    (w1,w2+w1^2), pulled back along its actual linear representation.
    """
    ks = _positive_parameters(k_values, "k_values")
    ns = _positive_parameters(n_values, "n_values")
    finite = {}

    def add(row, fermionic_group, factors, parameters=None, s="0", omega="0",
            omega_kind="zero", omega_cyclic_order=None):
        factors_in_paper_order = list(factors)
        canonical = tuple(sorted(order for order in factors if order != 1))
        if canonical not in finite:
            identifier = "x".join(f"C{order}" for order in canonical) or "1"
            finite[canonical] = {
                "id": f"table3-{identifier}",
                "kind": "abelian",
                "factors": list(canonical),
                "label": _group_label(canonical),
                "source": TABLE_III_SOURCE,
                "source_rows": [],
            }
        source_row = {
            "row": row,
            "fermionic_group": fermionic_group,
            "parameters": dict(parameters or {}),
            "quotient_factors_in_paper_order": factors_in_paper_order,
            "paper_s": s,
            "paper_omega": omega,
            "paper_omega_kind": omega_kind,
            "source": TABLE_III_SOURCE,
        }
        if omega_cyclic_order is not None:
            source_row["paper_omega_cyclic_order"] = omega_cyclic_order
            source_row["twist_source"] = CYCLIC_TWIST_SOURCE
        elif omega_kind == "quaternion":
            source_row["twist_source"] = QUATERNION_TWIST_SOURCE
        finite[canonical]["source_rows"].append(source_row)

    def add_carry(row, m, other_factors=(), parameters=None):
        fermionic_factors = [f"C{2 * m}^f"]
        fermionic_factors.extend(f"C{order}" for order in other_factors)
        add(
            row, " x ".join(fermionic_factors), [m, *other_factors],
            parameters=parameters,
            omega=(f"u_{m} on the first paper-order C{m} factor; "
                   f"u_{m}(a,b)=floor((a0+b0)/{m}) mod 2"),
            omega_kind="cyclic_carry", omega_cyclic_order=m,
        )

    add(1, "C2^f x C2", [2])
    for k in ks:
        add(2, f"C2^f x C{2 * k + 1}", [2 * k + 1], {"k": k})
    for k in ks:
        add(3, f"C2^f x C{2 * k}", [2 * k], {"k": k})

    split_products = (
        (4, (2, 2)),
        (5, (2, 4)),
        (6, (2, 8)),
        (7, (4, 4)),
        (8, (4, 8)),
        (9, (2, 2, 2)),
        (10, (2, 2, 4)),
        (11, (2, 4, 4)),
    )
    for row, factors in split_products:
        add(row, "C2^f x " + _group_label(factors), factors)

    for n in ns:
        k = 2 ** n
        add_carry(12, k, parameters={"n": n, "k": k})
    for row, m, other_factors in (
        (13, 2, (2,)),
        (14, 2, (4,)),
        (15, 4, (2,)),
        (16, 2, (2, 2)),
        (17, 2, (2, 4)),
        (18, 2, (4, 4)),
    ):
        add_carry(row, m, other_factors)

    add(19, "C2^f x C2^T", [2], s="x")
    add(20, "C4^(Tf)", [2], s="x", omega="x^2",
        omega_kind="cyclic_carry", omega_cyclic_order=2)
    add(21, "Q8^f", [2, 2], omega="x^2 + x*y + y^2",
        omega_kind="quaternion")

    result = [finite[factors] for factors in sorted(finite)]
    for number in range(1, 231):
        result.append({
            "id": f"SG{number:03d}",
            "kind": "spacegroup",
            "number": number,
            "label": f"Space group {number}",
            "source": SPACE_GROUP_SOURCE,
            "twist_source": SPACE_GROUP_TWIST_SOURCE,
            "spacegroup_twists": [
                {
                    "spacegroup_twist": "zero", "s": "w1", "omega": "0",
                    "description": "Orientation character with zero central-extension twist.",
                },
                {
                    "spacegroup_twist": "w2_plus_w1_squared",
                    "s": "w1", "omega": "w2+w1^2",
                    "description": (
                        "Pin-minus obstruction of the actual linear representation G to O(3); "
                        "pull back to the full space group, retaining translations."
                    ),
                },
            ],
            "electronic_twist_description": (
                "Effective internal twists for the spin-half electronic double "
                "space group: s=w1 (orientation character of the spatial linear "
                "part), omega=0. Use the full space group including translations."
            ),
        })
    return result
