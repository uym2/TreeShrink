from __future__ import annotations

import math
from collections import defaultdict, deque
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, FrozenSet, List, Optional, Set, Tuple


@dataclass
class Node:
    label: Optional[str] = None
    length: Optional[float] = None
    children: List["Node"] = field(default_factory=list)


class NewickError(ValueError):
    pass


class NewickParser:
    def __init__(self, text: str):
        self.text = text
        self.i = 0

    def parse(self) -> Node:
        self._skip_ws_and_comments()
        if self.i >= len(self.text):
            raise NewickError("empty Newick string")
        tree = self._parse_subtree()
        self._skip_ws_and_comments()
        if self._peek() == ";":
            self.i += 1
        self._skip_ws_and_comments()
        if self.i != len(self.text):
            raise NewickError(f"unexpected text at position {self.i}: {self.text[self.i:self.i + 20]!r}")
        return tree

    def _parse_subtree(self) -> Node:
        self._skip_ws_and_comments()
        if self._peek() == "(":
            self.i += 1
            children = [self._parse_subtree()]
            while True:
                self._skip_ws_and_comments()
                char = self._peek()
                if char == ",":
                    self.i += 1
                    children.append(self._parse_subtree())
                elif char == ")":
                    self.i += 1
                    break
                else:
                    raise NewickError(f"expected ',' or ')' at position {self.i}")
            label = self._parse_label(optional=True)
            length = self._parse_length()
            return Node(label=label, length=length, children=children)

        label = self._parse_label(optional=False)
        length = self._parse_length()
        return Node(label=label, length=length)

    def _parse_label(self, optional: bool) -> Optional[str]:
        self._skip_ws_and_comments()
        char = self._peek()
        if char is None or char in ":,();":
            if optional:
                return None
            raise NewickError(f"expected label at position {self.i}")

        if char == "'":
            self.i += 1
            parts: List[str] = []
            while self.i < len(self.text):
                char = self.text[self.i]
                self.i += 1
                if char == "'":
                    if self._peek() == "'":
                        self.i += 1
                        parts.append("'")
                    else:
                        return "".join(parts)
                else:
                    parts.append(char)
            raise NewickError("unterminated quoted label")

        start = self.i
        while self.i < len(self.text) and self.text[self.i] not in ":,();[] \t\r\n":
            self.i += 1
        label = self.text[start:self.i].strip()
        if not label and not optional:
            raise NewickError(f"expected label at position {start}")
        return label or None

    def _parse_length(self) -> Optional[float]:
        self._skip_ws_and_comments()
        if self._peek() != ":":
            return None
        self.i += 1
        self._skip_ws_and_comments()
        start = self.i
        while self.i < len(self.text) and self.text[self.i] not in ",();[] \t\r\n":
            self.i += 1
        raw = self.text[start:self.i]
        if not raw:
            raise NewickError(f"missing branch length at position {start}")
        try:
            return float(raw)
        except ValueError as exc:
            raise NewickError(f"invalid branch length {raw!r} at position {start}") from exc

    def _skip_ws_and_comments(self) -> None:
        while self.i < len(self.text):
            if self.text[self.i].isspace():
                self.i += 1
            elif self.text[self.i] == "[":
                self.i += 1
                depth = 1
                while self.i < len(self.text) and depth:
                    if self.text[self.i] == "[":
                        depth += 1
                    elif self.text[self.i] == "]":
                        depth -= 1
                    self.i += 1
                if depth:
                    raise NewickError("unterminated comment")
            else:
                break

    def _peek(self) -> Optional[str]:
        return self.text[self.i] if self.i < len(self.text) else None


@dataclass(frozen=True)
class TreeSignature:
    taxa: FrozenSet[str]
    topology_splits: FrozenSet[FrozenSet[str]]
    edge_lengths: Dict[FrozenSet[str], float]


def tree_signature(text: str) -> TreeSignature:
    root = NewickParser(text).parse()
    adjacency: Dict[int, Dict[int, Optional[float]]] = defaultdict(dict)
    taxa_by_node: Dict[int, str] = {}
    nodes_by_id: Dict[int, Node] = {}

    def visit(node: Node, parent_id: Optional[int] = None) -> None:
        node_id = id(node)
        nodes_by_id[node_id] = node
        if not node.children:
            if not node.label:
                raise NewickError("leaf without taxon label")
            if node.label in taxa_by_node.values():
                raise NewickError(f"duplicate taxon label: {node.label}")
            taxa_by_node[node_id] = node.label
        if parent_id is not None:
            adjacency[parent_id][node_id] = node.length
            adjacency[node_id][parent_id] = node.length
        else:
            adjacency[node_id]
        for child in node.children:
            visit(child, node_id)

    visit(root)
    suppress_unlabeled_degree_two_nodes(adjacency, taxa_by_node, nodes_by_id)

    taxa = frozenset(taxa_by_node.values())
    if not taxa:
        raise NewickError("tree has no taxa")

    edge_lengths: Dict[FrozenSet[str], float] = {}
    topology_splits: Set[FrozenSet[str]] = set()

    for u, neighbors in list(adjacency.items()):
        for v, length in list(neighbors.items()):
            if u > v:
                continue
            side = taxa_on_side(adjacency, taxa_by_node, u, v)
            split = canonical_split(side, taxa)
            if not split:
                continue
            if 1 < len(split) < len(taxa) - 1:
                topology_splits.add(split)
            edge_lengths[split] = 0.0 if length is None else length

    return TreeSignature(taxa=taxa, topology_splits=frozenset(topology_splits), edge_lengths=edge_lengths)


def suppress_unlabeled_degree_two_nodes(
    adjacency: Dict[int, Dict[int, Optional[float]]],
    taxa_by_node: Dict[int, str],
    nodes_by_id: Dict[int, Node],
) -> None:
    changed = True
    while changed:
        changed = False
        for node_id in list(adjacency):
            if node_id in taxa_by_node:
                continue
            node = nodes_by_id[node_id]
            if node.label:
                continue
            neighbors = list(adjacency[node_id])
            if len(neighbors) != 2:
                continue
            left, right = neighbors
            combined = combine_lengths(adjacency[node_id][left], adjacency[node_id][right])
            del adjacency[left][node_id]
            del adjacency[right][node_id]
            del adjacency[node_id]
            adjacency[left][right] = combined
            adjacency[right][left] = combined
            changed = True
            break


def combine_lengths(left: Optional[float], right: Optional[float]) -> Optional[float]:
    if left is None and right is None:
        return None
    return (0.0 if left is None else left) + (0.0 if right is None else right)


def taxa_on_side(
    adjacency: Dict[int, Dict[int, Optional[float]]],
    taxa_by_node: Dict[int, str],
    start: int,
    blocked: int,
) -> FrozenSet[str]:
    seen = {blocked}
    queue = deque([start])
    labels: Set[str] = set()
    while queue:
        node = queue.popleft()
        if node in seen:
            continue
        seen.add(node)
        if node in taxa_by_node:
            labels.add(taxa_by_node[node])
        for neighbor in adjacency[node]:
            if neighbor not in seen:
                queue.append(neighbor)
    return frozenset(labels)


def canonical_split(side: FrozenSet[str], taxa: FrozenSet[str]) -> FrozenSet[str]:
    other = taxa - side
    if len(other) < len(side):
        return frozenset(other)
    if len(side) < len(other):
        return side
    return min(side, frozenset(other), key=lambda labels: tuple(sorted(labels)))


def compare_signatures(left: TreeSignature, right: TreeSignature, length_tol: float) -> Tuple[str, List[str]]:
    if left.taxa != right.taxa:
        return "DIFFERENT taxa", [
            f"taxa only in left: {sorted(left.taxa - right.taxa)[:8]}",
            f"taxa only in right: {sorted(right.taxa - left.taxa)[:8]}",
        ]

    if left.topology_splits != right.topology_splits:
        return "DIFFERENT topology", [
            f"splits only in left: {format_splits(left.topology_splits - right.topology_splits)}",
            f"splits only in right: {format_splits(right.topology_splits - left.topology_splits)}",
        ]

    if set(left.edge_lengths) != set(right.edge_lengths):
        return "SAME topology, different edge set", [
            f"edge splits only in left: {format_splits(set(left.edge_lengths) - set(right.edge_lengths))}",
            f"edge splits only in right: {format_splits(set(right.edge_lengths) - set(left.edge_lengths))}",
        ]

    length_diffs = []
    for split in sorted(left.edge_lengths, key=lambda s: (len(s), tuple(sorted(s)))):
        left_length = left.edge_lengths[split]
        right_length = right.edge_lengths[split]
        if not math.isclose(left_length, right_length, rel_tol=0.0, abs_tol=length_tol):
            length_diffs.append((split, left_length, right_length, abs(left_length - right_length)))

    if length_diffs:
        details = [
            f"{format_split(split)}: left={left_length:.12g} right={right_length:.12g} absdiff={diff:.3g}"
            for split, left_length, right_length, diff in length_diffs[:5]
        ]
        if len(length_diffs) > 5:
            details.append(f"... {len(length_diffs) - 5} more branch length difference(s)")
        return "SAME topology, different lengths", details

    return "SAME topology + lengths", []


def format_splits(splits, limit: int = 3) -> str:
    sorted_splits = sorted(splits, key=lambda split: (len(split), tuple(sorted(split))))
    rendered = [format_split(split) for split in sorted_splits[:limit]]
    if len(sorted_splits) > limit:
        rendered.append(f"... {len(sorted_splits) - limit} more")
    return "; ".join(rendered) if rendered else "none"


def format_split(split: FrozenSet[str]) -> str:
    labels = sorted(split)
    if len(labels) > 8:
        labels = labels[:8] + [f"... {len(split) - 8} more"]
    return "{" + ",".join(labels) + "}"


def read_trees(path: Path) -> List[str]:
    return [line.strip() for line in path.read_text().splitlines() if line.strip()]


def compare_tree_files(left_path: Path, right_path: Path, length_tol: float = 1e-9) -> List[str]:
    left_trees = read_trees(left_path)
    right_trees = read_trees(right_path)
    mismatches = []

    if len(left_trees) != len(right_trees):
        mismatches.append(f"tree count differs: left={len(left_trees)} right={len(right_trees)}")

    for index, (left_tree, right_tree) in enumerate(zip(left_trees, right_trees), start=1):
        try:
            left_sig = tree_signature(left_tree)
            right_sig = tree_signature(right_tree)
            status, details = compare_signatures(left_sig, right_sig, length_tol)
        except NewickError as exc:
            status, details = "PARSE ERROR", [str(exc)]

        if status != "SAME topology + lengths":
            detail = "; ".join(details)
            mismatches.append(f"tree {index}: {status}: {detail}")

    return mismatches
