"""Parse and evaluate KEGG module DEFINITION expressions.

Grammar (KEGG Module documentation, "Module definition"; binding from tightest to loosest):

    unit         := KO | '(' alternatives ')' | '--'
    complex      := ['-'] unit (('+' | '-') unit)*        '+' = required subunit, '-' = optional
    sequence     := complex (' ' complex)*                 space = consecutive steps, all required
    alternatives := sequence (',' sequence)*               comma = alternatives, any one suffices

Examples verified against the KEGG definitions this module is used with:
    ((K00925 K00625),K01895)   AckA AND Pta, or Acs                         (M00357 step 1)
    K00200+K00201+K00202+K00203-K11261+(K00205,K11260,K00204)
                               FwdABCD, optional FwdE, and one of FwdF/G/H   (M00567 step 1)
    (((K03831,K03638) K03750),K15376)   (MogA or MoaB) then MoeA, or the fusion  (M00880 step 3)

`--` marks a step KEGG leaves without a KO; it can never be satisfied by a KO set and is
reported as such (none of the modules used by this pipeline contains one). A DEFINITION that
spans several lines is the space-joined sequence of its lines.

evaluate(present) reports each top-level step (the elements of the top-level sequence) as
satisfied or not, with the minimal set of KOs whose addition would satisfy it; a module is
complete when every step is satisfied (KEGG Mapper's "complete").
"""
import re

KO_RE = re.compile(r"K\d{5}")


class Node:
    kind = "node"


class KO(Node):
    kind = "ko"

    def __init__(self, ko):
        self.ko = ko

    def __repr__(self):
        return self.ko


class Undefined(Node):
    kind = "undefined"

    def __repr__(self):
        return "--"


class Complex(Node):
    """Subunits: list of (node, optional)."""
    kind = "complex"

    def __init__(self, parts):
        self.parts = parts

    def __repr__(self):
        out = []
        for i, (n, opt) in enumerate(self.parts):
            sep = "-" if opt else ("+" if i else "")
            out.append(f"{sep}{_wrap(n)}")
        return "".join(out)


class Seq(Node):
    kind = "seq"

    def __init__(self, items):
        self.items = items

    def __repr__(self):
        return " ".join(_wrap(i) for i in self.items)


class Alt(Node):
    kind = "alt"

    def __init__(self, items):
        self.items = items

    def __repr__(self):
        return ",".join(_wrap(i) for i in self.items)


def _wrap(n):
    """Parenthesise a child whose operators bind more loosely than its parent's, so that
    repr(parse(x)) re-parses to the same expression."""
    return f"({n!r})" if isinstance(n, (Seq, Alt)) else repr(n)


def _normalise(text: str) -> str:
    """Join lines, collapse runs of whitespace and drop spaces around ',', '+', '(' and ')'."""
    s = " ".join(line.strip() for line in str(text).splitlines() if line.strip())
    s = re.sub(r"\s+", " ", s).strip()
    s = re.sub(r"\s*,\s*", ",", s)
    s = re.sub(r"\s*\+\s*", "+", s)
    s = re.sub(r"\(\s+", "(", s)
    s = re.sub(r"\s+\)", ")", s)
    return s


class _Parser:
    def __init__(self, text):
        self.s = _normalise(text)
        self.i = 0
        if not self.s:
            raise ValueError("empty KEGG module definition")

    def peek(self, k=1):
        return self.s[self.i:self.i + k]

    def error(self, msg):
        raise ValueError(f"KEGG definition parse error at {self.i}: {msg} in {self.s!r}")

    def parse(self):
        node = self.alternatives()
        if self.i != len(self.s):
            self.error(f"unexpected {self.peek()!r}")
        return node

    def alternatives(self):
        items = [self.sequence()]
        while self.peek() == ",":
            self.i += 1
            items.append(self.sequence())
        return items[0] if len(items) == 1 else Alt(items)

    def sequence(self):
        items = [self.complex()]
        while self.peek() == " ":
            self.i += 1
            items.append(self.complex())
        return items[0] if len(items) == 1 else Seq(items)

    def complex(self):
        parts = []
        opt = False
        if self.peek(2) != "--" and self.peek() == "-":
            opt = True
            self.i += 1
        parts.append((self.unit(), opt))
        while True:
            c = self.peek()
            if c == "+" or (c == "-" and self.peek(2) != "--"):
                self.i += 1
                parts.append((self.unit(), c == "-"))
            else:
                break
        if len(parts) == 1 and not parts[0][1]:
            return parts[0][0]
        return Complex(parts)

    def unit(self):
        if self.peek(2) == "--":
            self.i += 2
            return Undefined()
        if self.peek() == "(":
            self.i += 1
            node = self.alternatives()
            if self.peek() != ")":
                self.error("missing ')'")
            self.i += 1
            return node
        m = KO_RE.match(self.s, self.i)
        if not m:
            self.error(f"expected KO, '(' or '--', got {self.peek()!r}")
        self.i = m.end()
        return KO(m.group(0))


def parse(text: str) -> Node:
    return _Parser(text).parse()


def kos(node: Node) -> set:
    """Every KO referenced by the expression (optional subunits included)."""
    if isinstance(node, KO):
        return {node.ko}
    if isinstance(node, Undefined):
        return set()
    if isinstance(node, Complex):
        return set().union(*(kos(n) for n, _ in node.parts))
    return set().union(*(kos(n) for n in node.items))


def satisfied(node: Node, present) -> bool:
    if isinstance(node, KO):
        return node.ko in present
    if isinstance(node, Undefined):
        return False
    if isinstance(node, Complex):
        return all(satisfied(n, present) for n, opt in node.parts if not opt)
    if isinstance(node, Seq):
        return all(satisfied(n, present) for n in node.items)
    return any(satisfied(n, present) for n in node.items)


def missing(node: Node, present) -> frozenset:
    """Minimal set of KOs whose addition satisfies the node (ties: first alternative).

    An Undefined step contributes the pseudo-KO '--' (it cannot be satisfied by KOs)."""
    if isinstance(node, KO):
        return frozenset() if node.ko in present else frozenset([node.ko])
    if isinstance(node, Undefined):
        return frozenset(["--"])
    if isinstance(node, (Complex, Seq)):
        subs = node.parts if isinstance(node, Complex) else [(n, False) for n in node.items]
        out = frozenset()
        for n, opt in subs:
            if not opt:
                out |= missing(n, present)
        return out
    best = None
    for n in node.items:
        m = missing(n, present)
        if best is None or len(m) < len(best):
            best = m
    return best


def steps(node: Node) -> list:
    """Top-level steps: the elements of a top-level sequence, else the node itself."""
    return list(node.items) if isinstance(node, Seq) else [node]


def evaluate(node: Node, present) -> dict:
    present = set(present)
    out = []
    for k, st in enumerate(steps(node), 1):
        out.append({"step": k, "expr": repr(st), "satisfied": satisfied(st, present),
                    "missing": sorted(missing(st, present))})
    n_ok = sum(s["satisfied"] for s in out)
    return {"steps": out, "n_steps": len(out), "n_satisfied": n_ok, "complete": n_ok == len(out)}


def complex_slots(node: Node):
    """(mandatory_slots, optional_slots) of a node that is a single complex, a flat sequence of
    KO steps, a single KO or a flat alternative of KOs; each slot is a sorted tuple of
    alternative KOs. Raises for any other shape (a relaxed rule needs a flat slot list).

    A flat `Seq` (consecutive steps, each one KO or one flat alternative) is treated as a list
    of mandatory slots, exactly like the mandatory subunits of a complex: both mean "all of
    these must be found", which is what a relaxation counts.
    """
    def slot(n):
        if isinstance(n, KO):
            return (n.ko,)
        if isinstance(n, Alt) and all(isinstance(x, KO) for x in n.items):
            return tuple(sorted(x.ko for x in n.items))
        raise ValueError(f"complex subunit {n!r} is not a KO or a flat KO alternative")
    if isinstance(node, Complex):
        mand = [slot(n) for n, opt in node.parts if not opt]
        opt = [slot(n) for n, o in node.parts if o]
        return mand, opt
    if isinstance(node, Seq):
        return [slot(n) for n in node.items], []
    return [slot(node)], []
