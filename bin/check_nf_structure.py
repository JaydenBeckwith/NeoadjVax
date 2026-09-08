#!/usr/bin/env python3
"""Static structural check for this repo's Nextflow (.nf) files.

Not a substitute for actually running Nextflow (see the CI stub-run job,
or `nextflow ... -stub-run` locally against each tests/*/smoke.config):
this only catches brace/paren/bracket imbalance, duplicate process/
workflow/def definitions, and include { X as Y } from '...' references
that do not resolve to a real definition in the target file. It exists
because Nextflow itself cannot always be installed where these checks
need to run (see docs/ARCHITECTURE.md); when Nextflow is available,
prefer -stub-run, which catches far more.

Usage: python3 bin/check_nf_structure.py [repo_root]
Exit status: 0 if clean, 1 if any issue was found.
"""
import re
import sys
import os
import glob


def strip_comments_strings(text):
    out = []
    i = 0
    n = len(text)
    in_sq = in_dq = in_tsq = in_tdq = False
    in_line_comment = in_block_comment = False
    while i < n:
        c = text[i]
        two = text[i:i + 2]
        three = text[i:i + 3]
        if in_line_comment:
            if c == '\n':
                in_line_comment = False
                out.append(c)
            i += 1
            continue
        if in_block_comment:
            if two == '*/':
                in_block_comment = False
                i += 2
            else:
                i += 1
            continue
        if in_tsq:
            if three == "'''":
                in_tsq = False
                i += 3
            else:
                i += 1
            continue
        if in_tdq:
            if three == '"""':
                in_tdq = False
                i += 3
            else:
                i += 1
            continue
        if in_sq:
            if c == '\\':
                i += 2
                continue
            if c == "'":
                in_sq = False
            i += 1
            continue
        if in_dq:
            if c == '\\':
                i += 2
                continue
            if c == '"':
                in_dq = False
            i += 1
            continue
        if two == '//':
            in_line_comment = True
            i += 2
            continue
        if two == '/*':
            in_block_comment = True
            i += 2
            continue
        if three == "'''":
            in_tsq = True
            i += 3
            continue
        if three == '"""':
            in_tdq = True
            i += 3
            continue
        if c == "'":
            in_sq = True
            i += 1
            continue
        if c == '"':
            in_dq = True
            i += 1
            continue
        out.append(c)
        i += 1
    return ''.join(out)


def main():
    root = sys.argv[1] if len(sys.argv) > 1 else os.path.normpath(
        os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir))
    files = glob.glob(os.path.join(root, "**", "*.nf"), recursive=True)
    files = [f for f in files if "/.nextflow/" not in f and "/work/" not in f
             and "/test-results/" not in f]

    errors = []
    defs = {}

    for f in files:
        text = open(f, encoding='utf-8').read()
        stripped = strip_comments_strings(text)
        bal = {'(': 0, ')': 0, '{': 0, '}': 0, '[': 0, ']': 0}
        for c in stripped:
            if c in bal:
                bal[c] += 1
        if bal['('] != bal[')'] or bal['{'] != bal['}'] or bal['['] != bal[']']:
            errors.append(f"{f}: unbalanced brackets (paren {bal['(']}/{bal[')']}, "
                           f"brace {bal['{']}/{bal['}']}, bracket {bal['[']}/{bal[']']})")
        for m in re.finditer(r'^\s*process\s+(\w+)\s*\{', text, re.M):
            name = m.group(1)
            if name in defs:
                errors.append(f"{f}: duplicate process definition {name} (also in {defs[name]})")
            defs[name] = f
        for m in re.finditer(r'^\s*workflow\s+(\w+)\s*\{', text, re.M):
            name = m.group(1)
            if name in defs:
                errors.append(f"{f}: duplicate workflow definition {name} (also in {defs[name]})")
            defs[name] = f
        for m in re.finditer(r'^\s*def\s+(\w+)\s*\(', text, re.M):
            name = m.group(1)
            if name in defs:
                errors.append(f"{f}: duplicate def {name} (also in {defs[name]})")
            defs[name] = f

    for f in files:
        text = open(f, encoding='utf-8').read()
        for m in re.finditer(r'include\s*\{([^}]*)\}\s*from\s*[\'"]([^\'"]+)[\'"]', text):
            names_part, relpath = m.groups()
            target = os.path.normpath(os.path.join(os.path.dirname(f), relpath))
            if not target.endswith('.nf'):
                target += '.nf'
            for piece in names_part.split(';'):
                piece = piece.strip()
                if not piece:
                    continue
                if ' as ' in piece:
                    orig, _alias = [p.strip() for p in piece.split(' as ')]
                else:
                    orig = piece.strip()
                if not os.path.exists(target):
                    errors.append(f"{f}: include target does not exist: {target}")
                    continue
                ttext = open(target, encoding='utf-8').read()
                if not (re.search(rf'^\s*process\s+{re.escape(orig)}\s*\{{', ttext, re.M)
                        or re.search(rf'^\s*workflow\s+{re.escape(orig)}\s*\{{', ttext, re.M)
                        or re.search(rf'^\s*def\s+{re.escape(orig)}\s*\(', ttext, re.M)):
                    errors.append(f"{f}: include {orig!r} not found as process/workflow/def in {target}")

    print(f"Checked {len(files)} .nf files, {len(defs)} process/workflow/def definitions found")
    if errors:
        print(f"\n{len(errors)} ISSUE(S):")
        for e in errors:
            print(" -", e)
        return 1
    print("0 errors")
    return 0


if __name__ == "__main__":
    sys.exit(main())
