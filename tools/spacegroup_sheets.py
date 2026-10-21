#!/usr/bin/env python3
"""The two space-group sheets of degree 4 and their comparison with FSPT_AHSS.

usage: python3 tools/spacegroup_sheets.py <campaign dir> <extra result JSON>... --fspt <FSPT_AHSS checkout> --out spacegroup

<campaign dir> is bench/spacegroups_full of the working copy (results/
sgNNN-<twist>.json with their .time logs, groups.csv); every extra JSON is a
result of bench/spacegroup_case.g that fills a case the campaign did not
finish (its layers are read from the E6 table: the rows are q = -4..0, the
cells p = 0..). Writes spin_half.csv, spin_zero.csv, the same as Markdown,
and README.md with the comparison summary. Groups are compared as sorted
elementary divisors (prime powers, 0 for Z).
"""
import argparse, collections, csv, json, os, re, sys

LAYERS = ('A', 'B', 'C', 'D')
COLUMNS = {'A': 'p+ip (1,0)', 'B': 'Kitaev (2,-1)', 'C': 'complex fermion (3,-2)', 'D': 'bosonic (5,-4)'}
E6_CELL = {'A': (4, 1), 'B': (3, 2), 'C': (2, 3), 'D': (0, 5)}   # (row, column) of the E6 table of degree 4
TWISTS = {'spin_half': ('s=w1, omega=0', 'space_groups_spin_half.csv', 'half'),
          'spin_zero': ('s=w1, omega=w2+w1^2', 'space_groups_spinless.csv', 'spinless')}


def prime_powers(n):
    out, p = [], 2
    while p * p <= n:
        if n % p == 0:
            q = 1
            while n % p == 0:
                n //= p; q *= p
            out.append(q)
        p += 1
    if n > 1:
        out.append(n)
    return out


def normal(inv):
    out = []
    for n in inv:
        out += [0] if n == 0 else prime_powers(n)
    return sorted(out)


def show(inv):
    if inv is None:
        return ''
    inv = normal(inv)
    if not inv:
        return '0'
    parts = []
    for n, c in sorted(collections.Counter(inv).items(), key=lambda x: (x[0] == 0, x[0])):
        parts.append(('Z' if n == 0 else f'Z{n}') + (f'^{c}' if c > 1 else ''))
    return ' x '.join(parts)


def wall_of(path):
    try:
        for line in open(path):
            if 'Elapsed (wall clock)' in line:
                parts = line.split()[-1].split(':')
                return round(sum(float(p) * 60 ** i for i, p in enumerate(reversed(parts))))
    except OSError:
        pass
    return ''


def layers_from_e6(e6):
    out = {}
    for name, (r, c) in E6_CELL.items():
        try:
            out[name] = e6[r][c]
        except (IndexError, TypeError):
            out[name] = None
    return out


def load_case(path, twist_default=None):
    d = json.load(open(path))
    layers = d.get('layers') or layers_from_e6(d.get('E6', []))
    time_file = os.path.splitext(path)[0] + '.time'
    if not os.path.exists(time_file):
        time_file = os.path.join(os.path.dirname(path), 'case.time')
    return dict(number=int(d['sg']), twist=d.get('twist', twist_default), status=d.get('status'),
                reason=d.get('reason'), invariants=d.get('invariants'), layers=layers,
                certificate=d.get('certificate', ''), commit=d.get('commit', '2ccde64'),
                host=d.get('host', ''), wall=wall_of(time_file))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('campaign'); ap.add_argument('extra', nargs='*')
    ap.add_argument('--fspt', required=True); ap.add_argument('--out', required=True)
    ap.add_argument('--extra-commit', default='')
    a = ap.parse_args()
    groups = {int(r['number']): (r['symbol'], int(r['point group order']))
              for r in csv.DictReader(open(os.path.join(a.campaign, 'groups.csv')))}
    ours = {}
    for twist in TWISTS:
        for n in range(1, 231):
            p = os.path.join(a.campaign, 'results', f'sg{n:03d}-{twist}.json')
            if os.path.exists(p):
                c = load_case(p); c['source'] = 'campaign'; ours[n, twist] = c
    for p in a.extra:
        c = load_case(p); c['source'] = 'branch'
        if a.extra_commit:
            c['commit'] = a.extra_commit
        ours[c['number'], c['twist']] = c
    fspt = {}
    for twist, (_, file, _) in TWISTS.items():
        for r in csv.DictReader(open(os.path.join(a.fspt, 'results', 'computed_examples', 'tables', file))):
            n = int(re.match(r'SG(\d+)', r['Gb']).group(1))
            fspt[n, twist] = {'A': json.loads(r['pip']), 'B': json.loads(r['majorana']),
                              'C': json.loads(r['complex_fermion']), 'D': json.loads(r['bosonic']),
                              'full': json.loads(r['stacking_group'])}
    os.makedirs(a.out, exist_ok=True)
    summary = {}
    differences = []
    for twist, (label, _, _) in TWISTS.items():
        rows = []
        counts = collections.Counter()
        for n in range(1, 231):
            sym, po = groups[n]
            c = ours.get((n, twist)); f = fspt.get((n, twist))
            row = collections.OrderedDict([('number', n), ('symbol', sym), ('point group order', po), ('twist', label)])
            if c:
                for L in LAYERS:
                    row[COLUMNS[L]] = show(c['layers'].get(L))
                row['full group'] = show(c['invariants']) if c['invariants'] is not None else ''
                row['status'] = c['status'] if c['status'] == 'computed' else f"{c['status']}: {str(c['reason'])[:120]}"
                row['certificate'] = c['certificate']; row['package commit'] = c['commit']
                row['wall seconds'] = c['wall']; row['host'] = c['host']; row['source'] = c['source']
            else:
                for L in LAYERS:
                    row[COLUMNS[L]] = ''
                row.update({'full group': '', 'status': 'not computed', 'certificate': '', 'package commit': '', 'wall seconds': '', 'host': '', 'source': ''})
            if f:
                for L in LAYERS:
                    row['FSPT_AHSS ' + COLUMNS[L]] = show(f[L])
                row['FSPT_AHSS full group'] = show(f['full'])
            else:
                for L in LAYERS:
                    row['FSPT_AHSS ' + COLUMNS[L]] = ''
                row['FSPT_AHSS full group'] = ''
            if c and c['status'] == 'computed' and f:
                bad = [L for L in LAYERS if normal(c['layers'][L] or []) != normal(f[L])]
                if normal(c['invariants']) != normal(f['full']):
                    bad.append('full group')
                verdict = 'identical' if not bad else 'differs: ' + ', '.join(COLUMNS.get(b, b) for b in bad)
                if bad:
                    differences.append((twist, n, sym, bad, row['full group'], row['FSPT_AHSS full group']))
            elif f and not c:
                verdict = 'not computed here'
            elif c and not f:
                verdict = 'no FSPT_AHSS entry'
            else:
                verdict = 'not computed here' if not c or c['status'] != 'computed' else ''
            row['agreement with FSPT_AHSS'] = verdict
            counts[verdict.split(':')[0]] += 1
            rows.append(row)
        with open(os.path.join(a.out, f'{twist}.csv'), 'w', newline='') as out:
            w = csv.DictWriter(out, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
        with open(os.path.join(a.out, f'{twist}.md'), 'w') as out:
            keys = ['number', 'symbol', 'point group order'] + [COLUMNS[L] for L in LAYERS] + ['full group', 'status', 'certificate', 'package commit', 'wall seconds', 'FSPT_AHSS full group', 'agreement with FSPT_AHSS']
            out.write(f'# Space groups, package degree 4, twist {label}\n\n| ' + ' | '.join(keys) + ' |\n|' + '|'.join(' --- ' if k not in ('number', 'point group order', 'wall seconds') else ' ---: ' for k in keys) + '|\n')
            for r in rows:
                out.write('| ' + ' | '.join(str(r[k]) for k in keys) + ' |\n')
        summary[twist] = counts
    with open(os.path.join(a.out, 'comparison_summary.json'), 'w') as out:
        json.dump({'summary': {t: dict(c) for t, c in summary.items()},
                   'differences': [dict(twist=t, number=n, symbol=s, fields=b, ours=o, fspt=f) for t, n, s, b, o, f in differences]}, out, indent=1)
    for t, c in summary.items():
        print(t, dict(c))
    for d in differences:
        print('DIFF', d)


if __name__ == '__main__':
    main()
