#!/usr/bin/env python3
"""Readable reference implementation of the binary attack in newRSDAttack.pdf.

This is deliberately a *verification* implementation:

* equations are represented as Python sets (symmetric difference is addition);
* the names B1,...,B5 and F1,F2'' follow Sections 3.2, 3.3 and 4;
* no complexity-estimation, row-cut, bit-packing, or experiment code is mixed in;
* every accepted answer is checked against Hx=s and the regular-weight condition.

It is suitable for toy instances.  attack.c remains the optimized experiment
implementation.
"""

from __future__ import annotations

import argparse
import random
from dataclasses import dataclass
from itertools import combinations
from typing import Dict, Hashable, Iterable, List, Optional, Sequence, Set, Tuple


CONST = ("const",)
Monomial = Tuple[Hashable, ...]
Equation = Set[Monomial]


class TracePrinter:
    """Human-readable matrices for a small reference attack.

    Positions are printed 1-based to match mathematical notation.  Internally
    the implementation continues to use Python's 0-based indices.
    """

    def __init__(self, p: "Parameters") -> None:
        self.p = p
        self.slot_positions: List[int] = []

    def section(self, title: str) -> None:
        print(f"\n{'=' * 78}\n{title}\n{'=' * 78}")

    def step(self, number: str, title: str) -> None:
        self.section(f"STEP {number}. {title}")

    def note(self, text: str = "") -> None:
        print(text)

    @staticmethod
    def positions(values: Iterable[int]) -> str:
        return "{" + ", ".join(str(value + 1) for value in sorted(values)) + "}"

    @staticmethod
    def grouped_bits(bits: Iterable[object], group_size: int = 5) -> str:
        """Render a bit row with a visual separator after every group_size bits."""
        text = "".join(map(str, bits))
        return " ".join(
            text[start : start + group_size]
            for start in range(0, len(text), group_size)
        )

    def set_slots(self, positions: Sequence[int]) -> None:
        self.slot_positions = list(positions)
        print("\nThe other variables slot map (used by polynomial matrices):")
        for slot, position in enumerate(self.slot_positions):
            print(f"  v{slot} = x{position + 1}  (block {position // self.p.beta + 1})")

    def term_name(self, term: Monomial, global_terms: bool = False) -> str:
        if term == CONST:
            return "1"
        if global_terms:
            return "*".join(f"x{int(index) + 1}" for index in term)
        return "*".join(
            f"x{self.slot_positions[int(index)] + 1}"
            for index in term
        )

    def equation_text(self, row: Equation, global_terms: bool = False) -> str:
        terms = sorted(
            (term for term in row if term != CONST),
            key=lambda term: (len(term), tuple(int(index) for index in term)),
        )
        lhs = " + ".join(self.term_name(term, global_terms) for term in terms) or "0"
        return f"{lhs} = {int(CONST in row)}"

    def expression_text(self, expression: Equation, global_terms: bool = False) -> str:
        terms = sorted(
            (term for term in expression if term != CONST),
            key=lambda term: (len(term), tuple(int(index) for index in term)),
        )
        names = [self.term_name(term, global_terms) for term in terms]
        if CONST in expression:
            names.append("1")
        return " + ".join(names) or "0"

    def system(
        self,
        title: str,
        rows: Sequence[Equation],
        columns: Optional[Sequence[Monomial]] = None,
        *,
        global_terms: bool = False,
        row_labels: Optional[Sequence[str]] = None,
    ) -> None:
        """Print an augmented GF(2) matrix and the equation represented by each row."""
        if row_labels is not None and len(row_labels) != len(rows):
            raise ValueError("row_labels must have the same length as rows")
        if columns is None:
            columns = sorted(
                {term for row in rows for term in row if term != CONST},
                key=lambda term: (len(term), tuple(int(index) for index in term)),
            )
        columns = list(columns)
        print(f"\n[{title}]")
        print(f"shape: {len(rows)} x ({len(columns)} + 1 RHS)")
        if columns:
            names = ", ".join(
                self.term_name(column, global_terms) for column in columns
            )
            print(f"column legend: c0:c{len(columns) - 1} = {names}")
        else:
            print("column legend: (no coefficient columns)")
        description = "source constraint" if row_labels is not None else "equation"
        print(f"matrix rows: coefficient bits | RHS | {description}")
        for index, row in enumerate(rows):
            bits = "".join("1" if column in row else "0" for column in columns)
            equation = (
                row_labels[index]
                if row_labels is not None
                else self.equation_text(row, global_terms)
            )
            print(
                f"  r{index:02d}: {self.grouped_bits(bits) or '-'} | {int(CONST in row)}"
                f" | {equation}"
            )

    def instance(self, instance: "Instance") -> None:
        self.section("Input binary RSD instance")
        p = instance.p
        print(f"(n,k,h,beta)=({p.n},{p.k},{p.h},{p.beta})")
        print("H | s:")
        for index, (row, rhs) in enumerate(zip(instance.H, instance.syndrome)):
            print(f"  h{index:02d}: {self.grouped_bits(row)} | {rhs}")
        print(f"planted e: {''.join(map(str, instance.planted))}")
        for block in range(p.h):
            start = block * p.beta
            end = start + p.beta
            print(
                f"  block {block + 1}: positions {start + 1}..{end}, "
                f"e={''.join(map(str, instance.planted[start:end]))}"
            )

    def parameter_summary(self) -> None:
        """Print the binary-paper quantities and the planned matrix dimensions."""
        p = self.p
        q = pre_row_cut_target(p)
        b1, b2, b3, b4 = p.counts()
        self.section("Parameter derivation and paper-step matrix plan")
        print("Basic and balanced-distribution values:")
        print(f"  beta = n/h = {p.n}/{p.h} = {q['beta']}")
        print(f"  nlin = (n-k)+h = ({p.n}-{p.k})+{p.h} = {q['nlin']}")
        print(f"  eta = eta1+eta2 = {p.eta1}+{p.eta2} = {q['eta']}")
        print(
            f"  (omega,flin)=divmod(nlin,h)=({q['omega']},{q['flin']})"
        )
        print(
            f"  (mu1,f1)=divmod(eta1,h)=({q['mu1']},{q['f1']}), "
            f"(mu,f)=divmod(eta,h)=({q['mu']},{q['f']}), "
            f"(nu,fpre)=divmod(eta_pre,h)=({q['nu']},{q['fpre']})"
        )
        print(f"  per-block B1={b1}, B2={b2}, B3={b3}, B4={b4}")

        print("\nFirst and second systems:")
        print(f"  n1 = k-h-eta1 = {p.k}-{p.h}-{p.eta1} = {q['n1']}")
        print(f"  m1 = sum_i C(beta-|B_i,1|,2) = {q['m1']}")
        print(f"  N1 = C(n1,2)+n1 = {q['N1']}")
        print(f"  Delta1 = free equations/monomials at round 1 = {q['Delta1']}")
        print(f"  delta1 = Delta1-zeta = {q['Delta1']}-{q['zeta']} = {q['delta1']}")
        print(
            f"  reduced round-1 size = "
            f"(m1-Delta1) x (N1-Delta1) = "
            f"{q['m1']-q['Delta1']} x {q['N1']-q['Delta1']}"
        )
        print(f"  n2 = k-h-eta = {p.k}-{p.h}-{q['eta']} = {q['n2']}")
        print(f"  m2 = sum_i C(beta-|B_i,1|-|B_i,4|,2) = {q['m2']}")
        print(f"  N2 = C(n2,2)+n2 = {q['N2']}")
        print(f"  Delta2 = free equations/monomials at round 2 = {q['Delta2']}")
        print(f"  delta2 = Delta2-zeta = {q['Delta2']}-{q['zeta']} = {q['delta2']}")
        print(
            f"  reduced round-2 rows before PRE removal = "
            f"m2-Delta2 = {q['m2']-q['Delta2']}"
        )
        print(f"  zeta = within-block PRE-PRE free monomials = {q['zeta']}")
        print(
            f"  tau = (m1-Delta1)-(m2-Delta2) = {q['tau']}"
        )
        print(f"  gamma = C(eta_pre,2) = {q['gamma']}")
        print(
            f"  N2' = N2-gamma-delta2 = {q['N2p']}, "
            f"N2'' = N2'-eta_pre = {q['N2pp']} (=U)"
        )
        print(
            f"  m2' = m2-gamma-delta2 = {q['m2p']}, "
            f"m2'' = m2'-eta_pre = {q['m2pp']}"
        )

        print("\nCorollary 2 backward row budgets:")
        print(f"  xi(U) = floor(U/2) = {q['xi']}")
        print(f"  hm2'' = U+xi(U) = {q['hm2pp']}")
        print(f"  hm2' = hm2''+eta_pre = {q['hm2p']}")
        print(f"  hm1' = hm2'+tau = {q['hm1p']}")
        print(f"  hm1 = hm1'+gamma+delta1 = {q['hm1']}")
        print(f"  Step 1-5 rows = hm1-Delta1 = {q['target']}")

        print("\nPlanned coefficient-matrix sizes (RHS is one additional column):")
        print(f"  Step 1-1 H'                  : {p.n-p.k} x {p.n-p.eta1}")
        print(f"  Step 1-2 H''                 : {q['nlin']} x {p.n-p.eta1}")
        print(f"  Step 1-3 H''_B2              : {q['nlin']} x {q['nlin']}")
        print(
            f"  Step 1-4 reduced F1          : {q['m1']-q['Delta1']} "
            f"x {q['N1']-q['Delta1']}"
        )
        print(
            f"  Step 1-5 pre row cut         : {q['target']} "
            f"x {q['N1']-q['Delta1']}"
        )
        print(
            f"  Step 2 after PRE-PRE cut   : "
            f"{q['target']-q['gamma']} x "
            f"{q['N1']-q['Delta1']-q['gamma']}"
        )
        print(
            "  Step 3-1 after B4=0           : Step 2 size minus actual "
            "0=0 rows/B4 columns"
        )
        print(
            f"  Step 4-1 after x_pre removal : actual Step 3-1 rows-{p.eta_pre} "
            f"x N2''({q['N2pp']})"
        )
        print(
            f"  Step 4-2 final row cut       : hm2'' x N2'' = "
            f"{q['hm2pp']} x {q['N2pp']}"
        )
        print(
            f"  Step 4-3 rank/solve          : {q['hm2pp']} x {q['N2pp']}"
        )


def spread(total: int, blocks: int, offset: int) -> List[int]:
    """Section 3.2 balanced cyclic distribution."""
    base, extra = divmod(total, blocks)
    out = [base] * blocks
    for i in range(extra):
        out[(offset + i) % blocks] += 1
    return out


@dataclass(frozen=True)
class Parameters:
    n: int
    k: int
    h: int
    eta1: int
    eta2: int
    eta_pre: int

    @property
    def beta(self) -> int:
        return self.n // self.h

    @property
    def nlin(self) -> int:
        # Section 4: n-k syndrome equations + h binary block-sum equations.
        return self.n - self.k + self.h

    def counts(self) -> Tuple[List[int], List[int], List[int], List[int]]:
        flin = self.nlin % self.h
        f1 = self.eta1 % self.h
        f = (self.eta1 + self.eta2) % self.h
        b2 = spread(self.nlin, self.h, 0)
        b1 = spread(self.eta1, self.h, flin)
        b4 = spread(self.eta2, self.h, flin + f1)
        b3 = spread(self.eta_pre, self.h, flin + f)
        return b1, b2, b3, b4

    def validate(self) -> None:
        if not (0 < self.h <= self.n and 0 < self.k < self.n):
            raise ValueError("need 0 < h <= n and 0 < k < n")
        if self.n % self.h:
            raise ValueError("h must divide n")
        if min(self.eta1, self.eta2, self.eta_pre) < 0:
            raise ValueError("eta values must be non-negative")
        b1, b2, b3, b4 = self.counts()
        actual_totals = tuple(map(sum, (b1, b2, b3, b4)))
        expected_totals = (self.eta1, self.nlin, self.eta_pre, self.eta2)
        if actual_totals != expected_totals:
            raise ValueError(
                "invalid B-set totals: "
                f"got {actual_totals}, expected {expected_totals}"
            )
        for b in range(self.h):
            if b1[b] + b2[b] + b3[b] + b4[b] > self.beta:
                raise ValueError(f"block {b}: B1+B2+B3+B4 exceeds beta")


@dataclass
class Instance:
    p: Parameters
    H: List[List[int]]
    syndrome: List[int]
    planted: List[int]


@dataclass
class Layout:
    role: List[str]
    B1: List[int]
    B2: List[int]
    B3: List[int]
    candidates: List[int]
    b4_count: List[int]


@dataclass
class Trace:
    outer_trials: int = 0
    inner_trials: int = 0
    singular_B2: int = 0
    wrong_or_unsolved: int = 0
    f1_rows: int = 0
    row_cut_rows: int = 0
    preprocessed_rows: int = 0
    final_rows: int = 0
    final_columns: int = 0


@dataclass
class PreprocessedSystem:
    """Step 2 RREF after deleting quadratic, then linear, pivot blocks."""

    quadratic_residual: List[Equation]
    pre_linear_pivots: Dict[Monomial, Equation]
    linear_residual: List[Equation]
    step7_columns: List[Monomial]


def generate_instance(p: Parameters, rng: random.Random) -> Instance:
    """Generate a planted binary RSD instance."""
    H = [[rng.randrange(2) for _ in range(p.n)] for _ in range(p.n - p.k)]
    planted = [0] * p.n
    for b in range(p.h):
        planted[b * p.beta + rng.randrange(p.beta)] = 1
    syndrome = [
        sum(a * x for a, x in zip(row, planted)) & 1
        for row in H
    ]
    return Instance(p, H, syndrome, planted)


def sample_layout(p: Parameters, rng: random.Random) -> Layout:
    """Choose B1, B2, B3; B4 is chosen later for each inner trial."""
    b1, b2, b3, b4 = p.counts()
    role = ["B5"] * p.n
    sets: Dict[str, List[int]] = {"B1": [], "B2": [], "B3": []}
    for b in range(p.h):
        positions = list(range(b * p.beta, (b + 1) * p.beta))
        rng.shuffle(positions)
        start = 0
        for name, amount in (("B1", b1[b]), ("B2", b2[b]), ("B3", b3[b])):
            chosen = positions[start : start + amount]
            start += amount
            for pos in chosen:
                role[pos] = name
            sets[name].extend(chosen)
    candidates = [i for i, name in enumerate(role) if name == "B5"]
    return Layout(role, sets["B1"], sets["B2"], sets["B3"], candidates, b4)


def split_B2_by_linear_source(
    p: Parameters,
    layout: Layout,
) -> Tuple[List[int], List[int]]:
    """Split B2 into n-k syndrome pivots and one block-sum pivot per block."""
    block_sum: List[int] = []
    for block in range(p.h):
        in_block = sorted(
            pos for pos in layout.B2 if pos // p.beta == block
        )
        if not in_block:
            raise ValueError(
                f"block {block + 1} has no B2 variable for its block-sum equation"
            )
        block_sum.append(in_block[0])
    block_sum_set = set(block_sum)
    syndrome = sorted(pos for pos in layout.B2 if pos not in block_sum_set)
    if len(syndrome) != p.n - p.k or len(block_sum) != p.h:
        raise ValueError(
            "invalid B2 split: expected "
            f"{p.n - p.k} syndrome and {p.h} block-sum variables"
        )
    return syndrome, block_sum


def rref(
    equations: Iterable[Equation], columns: Sequence[Monomial]
) -> Tuple[Dict[Monomial, Equation], List[Equation]]:
    """RREF on selected columns; return pivot rows and the remaining rows."""
    rows = [set(row) for row in equations if row]
    pivot_rows: Dict[Monomial, Equation] = {}
    top = 0
    for column in columns:
        pivot = next((r for r in range(top, len(rows)) if column in rows[r]), None)
        if pivot is None:
            continue
        rows[top], rows[pivot] = rows[pivot], rows[top]
        for r in range(len(rows)):
            if r != top and column in rows[r]:
                rows[r] ^= rows[top]
        pivot_rows[column] = rows[top]
        top += 1
    return pivot_rows, [row for row in rows[top:] if row]


def rref_preserving_size(
    equations: Iterable[Equation],
    columns: Sequence[Monomial],
) -> Tuple[Dict[Monomial, Equation], List[Equation]]:
    """RREF selected columns while retaining every input row, including zero rows."""
    rows = [set(row) for row in equations]
    pivot_rows: Dict[Monomial, Equation] = {}
    top = 0
    for column in columns:
        pivot = next((r for r in range(top, len(rows)) if column in rows[r]), None)
        if pivot is None:
            continue
        rows[top], rows[pivot] = rows[pivot], rows[top]
        for r in range(len(rows)):
            if r != top and column in rows[r]:
                rows[r] ^= rows[top]
        pivot_rows[column] = rows[top]
        top += 1
    return pivot_rows, rows


def build_affine_forms(
    instance: Instance,
    layout: Layout,
    trace: Optional[TracePrinter] = None,
) -> Optional[Tuple[List[Equation], Dict[int, int]]]:
    """Paper Steps 1-1 through 1-3: set B1=0 and express B2 affinely."""
    p = instance.p
    syndrome_equations: List[Equation] = []
    for row, rhs in zip(instance.H, instance.syndrome):
        eq = {(j,) for j, value in enumerate(row) if value and j not in layout.B1}
        if rhs:
            eq.add(CONST)
        syndrome_equations.append(eq)
    remaining = [(j,) for j in range(p.n) if j not in layout.B1]
    if trace:
        trace.step("1-1", "First zero guessing and constructing H'")
        trace.note(f"B1 (first zero guess) = {trace.positions(layout.B1)}")
        trace.note(
            f"H' coefficients: (n-k) x (n-eta1) = "
            f"{p.n-p.k} x {p.n-p.eta1}"
        )
        trace.system(
            "H' after deleting the B1 columns",
            syndrome_equations,
            remaining,
            global_terms=True,
        )

    # Binary-only equations sum_{j in block} x_j = 1 (Section 4.1).
    block_sum_equations: List[Equation] = []
    for b in range(p.h):
        eq = {
            (j,)
            for j in range(b * p.beta, (b + 1) * p.beta)
            if j not in layout.B1
        }
        eq.add(CONST)
        block_sum_equations.append(eq)
    equations = syndrome_equations + block_sum_equations

    b2_columns = [(j,) for j in layout.B2]
    if trace:
        trace.step("1-2", "Adding linear equations and constructing H''")
        trace.note(
            f"H'' coefficients: ((n-k)+h) x (n-eta1) = "
            f"{p.nlin} x {p.n-p.eta1}"
        )
        trace.system(
            "Added binary block-sum equations",
            block_sum_equations,
            remaining,
            global_terms=True,
        )
        trace.system(
            "H'' after setting x_B1=0 (syndrome rows + binary block-sum rows)",
            equations,
            remaining,
            global_terms=True,
        )
        trace.step("1-3", "Gaussian elimination on H'' and eliminating B2 variables using linear equations")
        trace.note(
            f"H''_B2 = {p.nlin} x {p.nlin}; "
            f"remaining variables n1=k-h-eta1={p.k-p.h-p.eta1}"
        )
        trace.system(
            "Selected square matrix of coefficient block (for checking invertibility)",
            [{term for term in row if term in set(b2_columns)} for row in equations],
            b2_columns,
            global_terms=True,
        )
    pivots, _ = rref(equations, b2_columns)
    if len(pivots) != p.nlin:
        if trace:
            trace.note(
                f"\nabove matrix's rank = {len(pivots)} < {p.nlin}: "
                "singular, so this outer partition is rejected."
            )
        return None  # Remark 1: resample when H_B2 is singular.

    free_positions = layout.B3 + layout.candidates
    slot = {pos: i for i, pos in enumerate(free_positions)}
    if trace:
        trace.note(f"\nabove matrix's rank = {len(pivots)}: full rank.")
        trace.system(
            "RREF rows of H''_B2",
            [pivots[column] for column in b2_columns],
            b2_columns + [(j,) for j in free_positions],
            global_terms=True,
        )
        trace.set_slots(free_positions)
    forms: List[Equation] = [set() for _ in range(p.n)]
    for pos in free_positions:
        forms[pos] = {(slot[pos],)}
    for pos in layout.B2:
        row = pivots[(pos,)]
        form: Equation = set()
        for term in row:
            if term == CONST:
                form.add(CONST)
            elif term != (pos,):
                other = term[0]
                if other in slot:
                    form.add((slot[other],))
        forms[pos] = form
    if trace:
        print("\nAffine form for every original variable:")
        for pos, form in enumerate(forms):
            if layout.role[pos] == "B1":
                rhs = "0  (B1 guess)"
            else:
                rhs = trace.expression_text(form)
            print(f"  x{pos + 1} = {rhs}")
    return forms, slot


def product(
    left: Equation,
    right: Equation,
    slot_block: Sequence[int],
) -> Equation:
    """Multiply two affine forms over GF(2), using x^2=x."""
    out: Equation = set()
    for a in left:
        for b in right:
            if a == CONST and b == CONST:
                term = CONST
            elif a == CONST:
                term = b
            elif b == CONST:
                term = a
            else:
                i, j = sorted((a[0], b[0]))
                if i == j:
                    term = (i,)
                elif slot_block[i] == slot_block[j]:
                    # Free equation x_i*x_j=0, eliminated exactly as below (5).
                    continue
                else:
                    term = (i, j)
            if term in out:
                out.remove(term)
            else:
                out.add(term)
    return out


def build_F1(
    instance: Instance,
    layout: Layout,
    forms: List[Equation],
    trace: Optional[TracePrinter] = None,
) -> List[Equation]:
    """Paper equation (5), with free equations/monomials eliminated."""
    p = instance.p
    free_positions = layout.B3 + layout.candidates
    slot_block = [pos // p.beta for pos in free_positions]
    rows: List[Equation] = []
    free_equations = 0
    free_pairs: List[Tuple[int, int]] = []
    kept_pairs: List[Tuple[int, int]] = []
    for b in range(p.h):
        positions = [
            j for j in range(b * p.beta, (b + 1) * p.beta)
            if layout.role[j] != "B1"
        ]
        for left, right in combinations(positions, 2):
            if layout.role[left] != "B2" and layout.role[right] != "B2":
                free_equations += 1
                free_pairs.append((left, right))
                continue  # free equation; its monomial is substituted by zero
            eq = product(forms[left], forms[right], slot_block)
            if eq:
                rows.append(eq)
                kept_pairs.append((left, right))
    if trace:
        q = pre_row_cut_target(p)
        trace.step(
            "1-4",
            "Constructing the algebraic system with x_i*x_j=0 within each block, excluding free monomials/equations",
        )
        trace.note(
            f"planned reduced F1: (m1-Delta1) x (N1-Delta1) = "
            f"({q['m1']}-{q['Delta1']}) x "
            f"({q['N1']}-{q['Delta1']}) = "
            f"{q['m1']-q['Delta1']} x {q['N1']-q['Delta1']}"
        )
        trace.note(
            f"free equations removed (Delta1) = {free_equations}"
        )
        trace.note(
            "free equation/monomial pairs = "
            + (
                ", ".join(f"x{left + 1}*x{right + 1}=0" for left, right in free_pairs)
                if free_pairs
                else "(none)"
            )
        )
        trace.note(f"F1 rows kept = {len(rows)}")
        trace.note(
            "The coefficient bits are computed after affine substitution; "
            "the rightmost label keeps the original unexpanded constraint."
        )
        trace.system(
            "F1 reduced Macaulay system",
            rows,
            row_labels=[
                f"x{left + 1} * x{right + 1} = 0"
                for left, right in kept_pairs
            ],
        )
    return rows


def pre_row_cut_target(p: Parameters) -> Dict[str, int]:
    """Corollary-2 backward quantities for the pre row cut in Step 1-5."""
    choose2 = lambda value: value * (value - 1) // 2
    b1, b2, b3, b4 = p.counts()
    free1 = [p.beta - b1[i] - b2[i] for i in range(p.h)]
    free2 = [
        p.beta - b1[i] - b2[i] - b4[i]
        for i in range(p.h)
    ]
    b5 = [free2[i] - b3[i] for i in range(p.h)]
    m1 = sum(choose2(p.beta - b1[i]) for i in range(p.h))
    Delta1 = sum(choose2(value) for value in free1)
    m2 = sum(choose2(p.beta - b1[i] - b4[i]) for i in range(p.h))
    Delta2 = sum(choose2(value) for value in free2)
    zeta = sum(choose2(value) for value in b3)
    tau = (m1 - Delta1) - (m2 - Delta2)
    b5_total = sum(b5)
    core_columns = (
        p.eta_pre * b5_total
        - sum(b3[i] * b5[i] for i in range(p.h))
        + choose2(b5_total)
        - sum(choose2(value) for value in b5)
        + b5_total
    )
    xi = core_columns // 2
    gamma = choose2(p.eta_pre)
    eta = p.eta1 + p.eta2
    n1 = p.k - p.h - p.eta1
    n2 = p.k - p.h - eta
    N1 = choose2(n1) + n1
    N2 = choose2(n2) + n2
    delta1 = Delta1 - zeta
    delta2 = Delta2 - zeta
    N2p = N2 - gamma - delta2
    N2pp = N2p - p.eta_pre
    m2p = m2 - gamma - delta2
    m2pp = m2p - p.eta_pre
    hm2pp = core_columns + xi
    hm2p = hm2pp + p.eta_pre
    hm1p = hm2p + tau
    hm1 = hm1p + gamma + delta1
    return {
        "beta": p.beta,
        "nlin": p.nlin,
        "eta": eta,
        "omega": p.nlin // p.h,
        "flin": p.nlin % p.h,
        "mu1": p.eta1 // p.h,
        "f1": p.eta1 % p.h,
        "mu": eta // p.h,
        "f": eta % p.h,
        "nu": p.eta_pre // p.h,
        "fpre": p.eta_pre % p.h,
        "n1": n1,
        "m1": m1,
        "N1": N1,
        "Delta1": Delta1,
        "delta1": delta1,
        "n2": n2,
        "m2": m2,
        "N2": N2,
        "N2p": N2p,
        "N2pp": N2pp,
        "m2p": m2p,
        "m2pp": m2pp,
        "Delta2": Delta2,
        "delta2": delta2,
        "zeta": zeta,
        "tau": tau,
        "gamma": gamma,
        "U": core_columns,
        "xi": xi,
        "hm1": hm1,
        "hm1p": hm1p,
        "hm2p": hm2p,
        "hm2pp": hm2pp,
        "target": hm1 - Delta1,
    }


def apply_pre_row_cut(
    rows: List[Equation],
    p: Parameters,
    rng: random.Random,
    trace: Optional[TracePrinter] = None,
) -> List[Equation]:
    """Randomly keep the Corollary-2 backward budget hm1-Delta1."""
    quantities = pre_row_cut_target(p)
    requested = quantities["target"]
    kept = min(requested, len(rows))
    indices = rng.sample(range(len(rows)), kept)
    cut = [set(rows[index]) for index in indices]
    if trace:
        trace.step("1-5", "Pre row cut (Corollary 2, backward)")
        trace.note(
            "backward quantities: "
            f"Delta1={quantities['Delta1']}, zeta={quantities['zeta']}, "
            f"tau={quantities['tau']}, gamma={quantities['gamma']}, "
            f"U={quantities['U']}, xi(U)={quantities['xi']}"
        )
        trace.note(
            f"random row selection: {len(rows)} -> {kept} rows "
            f"(requested hm1-Delta1={requested})"
        )
        if kept < requested:
            trace.note("row shortage: all available rows are retained.")
        trace.system("F1 after pre row cut", cut)
    return cut


def preprocess(
    rows: List[Equation],
    eta_pre: int,
    trace: Optional[TracePrinter] = None,
) -> PreprocessedSystem:
    """Eliminate PRE-PRE, then precompute the following x_pre elimination."""
    present = {term for row in rows for term in row}
    pp = [(i, j) for i in range(eta_pre) for j in range(i + 1, eta_pre) if (i, j) in present]
    pre_linear = [(i,) for i in range(eta_pre)]
    leading = set(pp) | set(pre_linear)
    remaining_columns = sorted(
        {
            term
            for row in rows
            for term in row
            if term != CONST and term not in leading
        },
        key=lambda term: (len(term), tuple(int(index) for index in term)),
    )
    step6_columns = pp + pre_linear + remaining_columns
    full_pivots, full_rref = rref_preserving_size(
        rows,
        pp + pre_linear,
    )
    pp_pivots = {
        column: full_pivots[column]
        for column in pp
        if column in full_pivots
    }
    pre_linear_pivots = {
        column: full_pivots[column]
        for column in pre_linear
        if column in full_pivots
    }
    pp_rank = len(pp_pivots)
    pre_rank = len(pre_linear_pivots)
    unpivoted_pp = [column for column in pp if column not in pp_pivots]
    step7_columns = unpivoted_pp + pre_linear + remaining_columns
    quadratic_residual = [
        set(row) for row in full_rref[pp_rank:]
    ]
    linear_residual = [
        set(row) for row in full_rref[pp_rank + pre_rank:]
    ]
    if trace:
        trace.step(
            "2",
            "Preprocessing: Gaussian elimination on PRE-PRE quadratic and x_pre linear columns",
        )
        b3_positions = trace.slot_positions[:eta_pre]
        trace.note(f"B3 (eta_pre) = {trace.positions(b3_positions)}")
        trace.note(
            "x_pre slot order = "
            + (
                ", ".join(
                    f"v{slot}=x{position + 1}"
                    for slot, position in enumerate(b3_positions)
                )
                if b3_positions
                else "(none)"
            )
        )
        trace.note(
            "PRE-PRE quadratic columns moved to the front: "
            + (", ".join(trace.term_name(term) for term in pp) or "(none)")
        )
        group_notes: List[str] = []
        if pp:
            group_notes.append(f"c0:c{len(pp) - 1}=PRE-PRE quadratic (eliminate now)")
        if pre_linear:
            start = len(pp)
            group_notes.append(
                f"c{start}:c{start + len(pre_linear) - 1}=x_pre linear "
                "(keep until Step 4-1)"
            )
        if remaining_columns:
            start = len(pp) + len(pre_linear)
            group_notes.append(
                f"c{start}:c{len(step6_columns) - 1}=remaining monomials"
            )
        trace.note("column groups: " + "; ".join(group_notes))
        trace.system(
            "1) Step 2 reordered input matrix",
            rows,
            step6_columns,
        )
        trace.system(
            "2) Same-size RREF on B3 quadratic and linear columns",
            full_rref,
            step6_columns,
        )
        trace.note(
            "Each pivoted B3 target column has exactly one 1 in the "
            "same-size RREF above."
        )
        trace.system(
            "3) Reduced matrix after deleting PRE-PRE columns and pivot rows "
            "(passed to Step 3-1)",
            quadratic_residual,
            step7_columns,
        )
        trace.note(
            f"Step 2 ranks: PRE-PRE={len(pp_pivots)}/{len(pp)}, "
            f"x_pre={len(pre_linear_pivots)}/{len(pre_linear)}; "
            f"same-size={len(full_rref)}x{len(step6_columns)}, "
            f"after PRE-PRE deletion={len(quadratic_residual)}"
            f"x{len(step7_columns)}, "
            f"after later x_pre deletion={len(linear_residual)}"
            f"x{len(step7_columns) - pre_rank}"
        )
    return PreprocessedSystem(
        quadratic_residual,
        pre_linear_pivots,
        linear_residual,
        step7_columns,
    )


def choose_B4(p: Parameters, layout: Layout, rng: random.Random) -> Set[int]:
    chosen: Set[int] = set()
    for b, amount in enumerate(layout.b4_count):
        available = [
            pos for pos in layout.candidates
            if pos // p.beta == b
        ]
        chosen.update(rng.sample(available, amount))
    return chosen


def substitute_B4_with_zero_count(
    rows: Iterable[Equation],
    layout: Layout,
    slot: Dict[int, int],
    B4: Set[int],
) -> Tuple[List[Equation], int]:
    """Apply B4=0 and count rows that become the zero equation 0=0."""
    zero_slots = {slot[pos] for pos in B4}
    result: List[Equation] = []
    zero_rows = 0
    for row in rows:
        new_row = {
            term for term in row
            if term == CONST or not any(index in zero_slots for index in term)
        }
        if new_row:
            result.append(new_row)
        else:
            zero_rows += 1
    return result, zero_rows


def substitute_B4(
    rows: Iterable[Equation],
    layout: Layout,
    slot: Dict[int, int],
    B4: Set[int],
) -> List[Equation]:
    return substitute_B4_with_zero_count(rows, layout, slot, B4)[0]


def columns_after_B4(
    columns: Sequence[Monomial],
    slot: Dict[int, int],
    B4: Set[int],
) -> List[Monomial]:
    """Keep the existing column order while deleting B4-related monomials."""
    zero_slots = {slot[pos] for pos in B4}
    return [
        term
        for term in columns
        if not any(index in zero_slots for index in term)
    ]


def solve_F2pp(
    linear_residual: List[Equation],
    shared_pre_pivots: Dict[Monomial, Equation],
    step7_columns: Sequence[Monomial],
    layout: Layout,
    slot: Dict[int, int],
    B4: Set[int],
    rng: random.Random,
    trace: Optional[TracePrinter] = None,
) -> Optional[Tuple[Dict[int, int], int]]:
    """Project the Step 2 x_pre elimination by B4, solve, and recover x_pre."""
    eta_pre = len(layout.B3)
    pre_columns = [(i,) for i in range(eta_pre)]
    projected_columns = columns_after_B4(step7_columns, slot, B4)
    core_columns = [
        column for column in projected_columns if column not in pre_columns
    ]
    projected_pivot_rows = substitute_B4(
        [shared_pre_pivots[column] for column in pre_columns if column in shared_pre_pivots],
        layout,
        slot,
        B4,
    )
    pre_pivots = {
        column: row
        for column, row in zip(
            [column for column in pre_columns if column in shared_pre_pivots],
            projected_pivot_rows,
        )
    }
    F2pp = substitute_B4(
        linear_residual,
        layout,
        slot,
        B4,
    )
    if trace:
        trace.step(
            "4-1",
            "Constructing F2'' by applying B4=0 to the residual system after x_pre elimination in Step 2",
        )
        trace.note(
            "The B3 linear columns already have isolated pivots from Step 2. "
            "Delete those columns and their pivot rows; no new elimination is needed."
        )
        trace.note(
            f"x_pre deletion: rows {len(F2pp) + len(pre_pivots)}"
            f"->{len(F2pp)}, columns {len(projected_columns)}"
            f"->{len(core_columns)} (=N2'')"
        )
        trace.note(
            "\nFinal core columns z2'': "
            + (", ".join(trace.term_name(term) for term in core_columns) or "(none)")
        )
        trace.system(
            "2) Reduced matrix after deleting B3 linear columns and pivot rows",
            F2pp,
            core_columns,
        )
    if len(pre_pivots) != eta_pre:
        if trace:
            trace.note(
                f"x_pre rank {len(pre_pivots)} < eta_pre {eta_pre}: reject this guess."
            )
        return None

    if trace:
        trace.step("4-2", "Corollary 2 final row cut: U+xi(U) rows")
    requested_rows = len(core_columns) + len(core_columns) // 2
    available_rows = [
        row
        for row in F2pp
        if any(term == CONST or term in core_columns for term in row)
    ]
    if len(available_rows) < requested_rows:
        if trace:
            trace.note(
                f"final row shortage: available {len(available_rows)} "
                f"< requested U+xi(U)={requested_rows}; reject this guess."
            )
        return None
    selected_rows = rng.sample(available_rows, requested_rows)
    if trace:
        trace.note(
            f"final Corollary-2 row selection: {len(available_rows)} -> "
            f"{requested_rows} rows (U={len(core_columns)}, "
            f"xi(U)={len(core_columns) // 2})"
        )
        trace.system("Selected final system F2''", selected_rows, core_columns)

    pivots, full_rref = rref_preserving_size(selected_rows, core_columns)
    residual = full_rref[len(pivots):]
    zero_rows = sum(not row for row in full_rref)
    inconsistent_rows = sum(row == {CONST} for row in full_rref)
    if trace:
        trace.step("4-3", "Solving the final system, reconstructing the candidate, and checking the solution")
        trace.note(
            f"final system dimensions: hm2'' x N2'' = "
            f"{len(selected_rows)} x {len(core_columns)}"
        )
        trace.system(
            "F2'' full RREF after final Gaussian elimination (all rows)",
            full_rref,
            core_columns,
        )
        trace.note(
            f"RREF row classification: pivot={len(pivots)}, "
            f"0=0 rows={zero_rows}, 0=1 rows={inconsistent_rows}"
        )
        trace.note(
            "consistency: "
            + ("INCONSISTENT (a 0=1 row exists)" if inconsistent_rows else "consistent")
        )
    if inconsistent_rows:
        if trace:
            trace.note("Residual row 0=1 found: the final system is inconsistent.")
        return None
    if len(pivots) != len(core_columns):
        if trace:
            trace.note(
                f"Final rank {len(pivots)} < {len(core_columns)} columns: "
                "consistent but not unique."
            )
        return None
    if trace:
        trace.note(
            f"final rank = {len(pivots)}/{len(core_columns)}: "
            "full rank and unique."
        )

    linearized_solution = {
        column: int(CONST in pivots[column])
        for column in core_columns
    }
    if trace:
        trace.note("\nUnique linearized solution z2'':")
        for column in core_columns:
            trace.note(
                f"  {trace.term_name(column)} = {linearized_solution[column]}"
            )
    values: Dict[int, int] = {
        column[0]: int(CONST in pivots[column])
        for column in core_columns
        if len(column) == 1
    }
    # Each pre pivot is x_pre + known core terms = rhs.
    for column in reversed(pre_columns):
        row = pre_pivots[column]
        value = int(CONST in row)
        for term in row:
            if term in (CONST, column):
                continue
            if len(term) == 1:
                value ^= values.get(term[0], 0)
            else:
                value ^= values.get(term[0], 0) * values.get(term[1], 0)
        values[column[0]] = value
    if trace:
        trace.note("\nRecovered free-variable values:")
        for slot_index in sorted(values):
            trace.note(
                f"  v{slot_index}=x{trace.slot_positions[slot_index] + 1}"
                f" = {values[slot_index]}"
            )
    return values, len(selected_rows)


def reconstruct(forms: List[Equation], values: Dict[int, int]) -> List[int]:
    answer: List[int] = []
    for form in forms:
        value = int(CONST in form)
        for term in form:
            if term != CONST:
                value ^= values.get(term[0], 0)
        answer.append(value)
    return answer


def verify(instance: Instance, answer: Sequence[int]) -> bool:
    p = instance.p
    if any(sum(answer[b * p.beta : (b + 1) * p.beta]) != 1 for b in range(p.h)):
        return False
    syndrome = [
        sum(a * x for a, x in zip(row, answer)) & 1
        for row in instance.H
    ]
    return syndrome == instance.syndrome


def attack(
    instance: Instance,
    rng: random.Random,
    outer_budget: int,
    inner_budget: int,
    trace: Optional[TracePrinter] = None,
) -> Tuple[Optional[List[int]], Trace]:
    """Run the actual two-stage randomized attack (no planted-solution oracle)."""
    stats = Trace()
    if isinstance(trace, TracePrinter):
        trace.parameter_summary()
    for outer_index in range(outer_budget):
        stats.outer_trials += 1
        layout = sample_layout(instance.p, rng)
        if isinstance(trace, TracePrinter):
            trace.section(f"Outer trial {outer_index + 1}: balanced B1/B2/B3 partition")
            b1, b2, b3, b4 = instance.p.counts()
            print("Per-block target counts:")
            print("  block | B1 | B2 | B3 | B4 | B5")
            b5_counts = []
            for block in range(instance.p.h):
                b5 = instance.p.beta - b1[block] - b2[block] - b3[block] - b4[block]
                b5_counts.append(b5)
                print(
                    f"  {block + 1:5d} | {b1[block]:2d} | {b2[block]:2d}"
                    f" | {b3[block]:2d} | {b4[block]:2d} | {b5:2d}"
                )
            print(
                f"  total | {sum(b1):2d} | {sum(b2):2d} | {sum(b3):2d}"
                f" | {sum(b4):2d} | {sum(b5_counts):2d}"
            )
            trace.note(f"B1 = {trace.positions(layout.B1)}")
            b2_syndrome, b2_block_sum = split_B2_by_linear_source(
                instance.p,
                layout,
            )
            trace.note(
                "B2 = {syndrome-equation pivots} U {h block-sum equation pivots}"
            )
            trace.note(
                f"   = {trace.positions(b2_syndrome)}"
                f" U {trace.positions(b2_block_sum)}"
            )
            trace.note(f"B3 = {trace.positions(layout.B3)}")
            trace.note(f"B4 candidates = {trace.positions(layout.candidates)}")
            wrong_b1 = [pos + 1 for pos in layout.B1 if instance.planted[pos]]
            trace.note(
                "B1 guess is "
                + (
                    "WRONG; contains planted nonzero " + str(wrong_b1)
                    if wrong_b1
                    else "correct"
                )
            )

        built = build_affine_forms(instance, layout, trace if isinstance(trace, TracePrinter) else None)
        if built is None:
            stats.singular_B2 += 1
            continue
        forms, slot = built
        F1 = build_F1(instance, layout, forms, trace if isinstance(trace, TracePrinter) else None)
        cut_F1 = apply_pre_row_cut(
            F1,
            instance.p,
            rng,
            trace if isinstance(trace, TracePrinter) else None,
        )
        preprocessed = preprocess(
            cut_F1,
            instance.p.eta_pre,
            trace if isinstance(trace, TracePrinter) else None,
        )
        shared = preprocessed.quadratic_residual
        stats.f1_rows = len(F1)
        stats.row_cut_rows = len(cut_F1)
        stats.preprocessed_rows = len(shared)

        tries = 1 if instance.p.eta2 == 0 else inner_budget
        if isinstance(trace, TracePrinter):
            trace.step(
                "3-1",
                "Second zero guessing: choosing B4 in the inner loop",
            )
            trace.note(
                f"inner trials planned = {tries}; "
                f"B4 candidates = {trace.positions(layout.candidates)}"
            )
        for inner_index in range(tries):
            stats.inner_trials += 1
            B4 = choose_B4(instance.p, layout, rng)
            if isinstance(trace, TracePrinter):
                trace.section(
                    f"STEP 3-1 inner trial {inner_index + 1}/{tries}: "
                    "choose B4 and set x_B4=0"
                )
                B5 = set(layout.candidates) - B4
                trace.note(f"B4 (guessed zero) = {trace.positions(B4)}")
                trace.note(f"B5 (surviving)    = {trace.positions(B5)}")
                wrong = [
                    pos + 1
                    for pos in B4
                    if instance.planted[pos]
                ]
                trace.note(
                    "B4 guess is "
                    + ("WRONG; contains planted nonzero " + str(wrong) if wrong else "correct")
                )
            F2, zero_rows = substitute_B4_with_zero_count(
                shared,
                layout,
                slot,
                B4,
            )
            if isinstance(trace, TracePrinter):
                trace.note(
                    f"0=0 rows (less and equal to tau) created by this B4 guess: {zero_rows}; "
                    f"row count {len(shared)} -> {len(F2)}"
                )
                trace.system(
                    "1) Matrix reduced by the second zero guess "
                    "(Step 2 column order retained)",
                    F2,
                    columns_after_B4(
                        preprocessed.step7_columns,
                        slot,
                        B4,
                    ),
                )
            solved = solve_F2pp(
                preprocessed.linear_residual,
                preprocessed.pre_linear_pivots,
                preprocessed.step7_columns,
                layout,
                slot,
                B4,
                rng,
                trace if isinstance(trace, TracePrinter) else None,
            )
            if solved is None:
                stats.wrong_or_unsolved += 1
                continue
            values, final_row_count = solved
            answer = reconstruct(forms, values)
            if verify(instance, answer):
                stats.final_rows = final_row_count
                stats.final_columns = len(values)
                if isinstance(trace, TracePrinter):
                    trace.note("\nCandidate reconstruction and verification:")
                    trace.note(f"reconstructed e: {''.join(map(str, answer))}")
                    trace.note(f"planted e:       {''.join(map(str, instance.planted))}")
                    trace.note(f"block weight one: {all(sum(answer[b * instance.p.beta:(b + 1) * instance.p.beta]) == 1 for b in range(instance.p.h))}")
                    trace.note(f"H*e = syndrome:   {verify(instance, answer)}")
                return answer, stats
            stats.wrong_or_unsolved += 1
            if isinstance(trace, TracePrinter):
                trace.note("Linearized solution failed the final RSD verification.")
    return None, stats


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("n", type=int)
    parser.add_argument("k", type=int)
    parser.add_argument("h", type=int)
    parser.add_argument("--eta1", type=int, required=True)
    parser.add_argument("--eta2", type=int, required=True)
    parser.add_argument("--pre", dest="eta_pre", type=int, required=True)
    parser.add_argument("--seed", type=int, default=12345)
    parser.add_argument("--outer-budget", type=int, default=2000)
    parser.add_argument("--inner-budget", type=int, default=2000)
    parser.add_argument(
        "--trace",
        action="store_true",
        help="print every paper-stage matrix (educational mode, requires n<=20)",
    )
    args = parser.parse_args()

    p = Parameters(args.n, args.k, args.h, args.eta1, args.eta2, args.eta_pre)
    p.validate()
    if args.trace and p.n > 20:
        parser.error("--trace is limited to n<=20; use the normal summary for larger instances")
    instance_rng = random.Random(args.seed)
    attack_rng = random.Random(args.seed ^ 0xA5A5A5A5)
    instance = generate_instance(p, instance_rng)
    matrix_trace = TracePrinter(p) if args.trace else None
    if matrix_trace:
        matrix_trace.instance(instance)
    answer, trace = attack(
        instance,
        attack_rng,
        args.outer_budget,
        args.inner_budget,
        matrix_trace,
    )

    print(f"RSD({p.n},{p.k},{p.h}), beta={p.beta}")
    print(f"(eta1, eta, eta_pre)=({p.eta1},{p.eta1 + p.eta2},{p.eta_pre})")
    print(
        f"trials: outer={trace.outer_trials}, inner={trace.inner_trials}, "
        f"singular_B2={trace.singular_B2}"
    )
    print(
        f"rows: F1={trace.f1_rows}, after_row_cut={trace.row_cut_rows}, "
        f"after_preprocessing={trace.preprocessed_rows}, "
        f"final={trace.final_rows}"
    )
    if answer is None:
        print("FAILED within the trial budgets")
        return 1
    print(f"verified: {verify(instance, answer)}")
    print(f"exact planted solution: {answer == instance.planted}")
    if p.n <= 80:
        print("planted :", "".join(map(str, instance.planted)))
        print("recovered:", "".join(map(str, answer)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
