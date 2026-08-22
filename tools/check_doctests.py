#!/usr/bin/env python3
"""Fail if the doctest suite shrinks.

These doctests are the only characterisation tests on the numeric core, and
the baseline Phase 3 validates the niwrap rewrite against. Two checks:

1. No docstring stubbed to `>>> pass`. A doctest that cannot run goes in the
   skip ledger in pyproject.toml with a reason.
2. Real doctest count stays at or above FLOOR.

Usage: python tools/check_doctests.py
"""

from __future__ import annotations

import ast
import doctest
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent / "src" / "mindboggle"

# Real (non-stub) doctests currently in the tree. Ratchet: raise when you add.
FLOOR = 145

_parser = doctest.DocTestParser()


def _docstrings(path: Path):
    """Yield (qualified_name, docstring) for every def/class in a file."""
    try:
        tree = ast.parse(path.read_text(encoding="utf-8"))
    except (SyntaxError, UnicodeDecodeError):
        return  # quarantined Python 2 remnants under x/; triaged separately

    def walk(node, prefix):
        for child in ast.iter_child_nodes(node):
            if isinstance(child, ast.FunctionDef | ast.AsyncFunctionDef | ast.ClassDef):
                qname = f"{prefix}{child.name}"
                docstring = ast.get_docstring(child)
                if docstring is not None:
                    yield qname, docstring
                yield from walk(child, f"{qname}.")

    yield from walk(tree, "")


def main() -> int:
    stubs: list[str] = []
    real = 0

    for path in sorted(ROOT.rglob("*.py")):
        rel = path.relative_to(ROOT.parent.parent).as_posix()
        for qname, docstring in _docstrings(path):
            try:
                examples = _parser.get_examples(docstring)
            except ValueError:
                continue
            if not examples:
                continue
            if len(examples) == 1 and examples[0].source.strip() == "pass":
                stubs.append(f"{rel}::{qname}")
            else:
                real += 1

    ok = True

    if stubs:
        ok = False
        print(f"ERROR: {len(stubs)} stubbed docstring(s) (`>>> pass`):\n")
        for stub in stubs:
            print(f"  {stub}")
        print(
            "\nA doctest that cannot run goes in the skip ledger in pyproject.toml\n"
            "with a reason and an issue, or gets `# doctest: +SKIP` on the line.\n"
            "Keep the expected output -- that is the part we cannot reconstruct.\n"
        )

    if real < FLOOR:
        ok = False
        print(
            f"ERROR: {real} real doctests, below the floor of {FLOOR}.\n\n"
            "Tests were removed. If that is genuinely intended, lower FLOOR in\n"
            "tools/check_doctests.py in its own commit, so it shows up in review.\n"
        )

    if ok:
        print(f"OK: {real} real doctests, no stubs (floor {FLOOR}).")
        if real > FLOOR:
            print(f"NOTE: {real} > floor {FLOOR}; consider raising FLOOR.")

    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
