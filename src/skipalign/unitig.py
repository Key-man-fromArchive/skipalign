"""Compacted De Bruijn graph construction and unitig extraction."""

from __future__ import annotations

from collections import defaultdict


def build_debruijn(kmers: list[str]) -> dict[str, list[str]]:
    """Build a directed De Bruijn graph from k-mers using (k-1) overlaps.

    Each k-mer is an edge: prefix → suffix where prefix and suffix are (k-1)-mers.
    """
    graph: dict[str, list[str]] = defaultdict(list)
    for kmer in kmers:
        prefix = kmer[:-1]
        suffix = kmer[1:]
        graph[prefix].append(suffix)
    return dict(graph)


def _compute_in_degree(graph: dict[str, list[str]]) -> dict[str, int]:
    """Compute in-degree for each node."""
    in_deg: dict[str, int] = defaultdict(int)
    for node in graph:
        if node not in in_deg:
            in_deg[node] = 0
        for neighbor in graph[node]:
            in_deg[neighbor] += 1
    return dict(in_deg)


def extract_unitigs(kmers: list[str]) -> list[str]:
    """Extract unitigs by compacting non-branching paths in the De Bruijn graph.

    A unitig is a maximal non-branching path: each internal node has
    in-degree == 1 and out-degree == 1.
    """
    if not kmers:
        return []

    graph = build_debruijn(kmers)
    in_deg = _compute_in_degree(graph)
    out_deg = {node: len(neighbors) for node, neighbors in graph.items()}

    # All nodes in the graph
    all_nodes = set(graph.keys())
    for neighbors in graph.values():
        all_nodes.update(neighbors)

    # A node is a "branch point" (unitig start) if it's NOT a simple pass-through
    def is_branch(node: str) -> bool:
        return in_deg.get(node, 0) != 1 or out_deg.get(node, 0) != 1

    visited_edges: set[tuple[str, str]] = set()
    unitigs: list[str] = []

    # Start unitig traversal from every branch point that has outgoing edges
    for start in all_nodes:
        if not is_branch(start):
            continue
        if start not in graph:
            continue
        for next_node in graph[start]:
            if (start, next_node) in visited_edges:
                continue
            # Walk the non-branching path
            path = [start, next_node]
            visited_edges.add((start, next_node))
            current = next_node
            while not is_branch(current) and current in graph:
                nxt = graph[current][0]
                visited_edges.add((current, nxt))
                path.append(nxt)
                current = nxt
            # Reconstruct sequence: first node + last char of each subsequent node
            k_minus_1 = len(path[0])
            seq = path[0] + "".join(node[-1] for node in path[1:])
            unitigs.append(seq)

    # Handle isolated cycles (nodes never visited because all have in/out = 1)
    for node in all_nodes:
        if node not in graph:
            continue
        for next_node in graph[node]:
            if (node, next_node) not in visited_edges:
                path = [node, next_node]
                visited_edges.add((node, next_node))
                current = next_node
                while current != node and current in graph:
                    nxt = graph[current][0]
                    if (current, nxt) in visited_edges:
                        break
                    visited_edges.add((current, nxt))
                    path.append(nxt)
                    current = nxt
                seq = path[0] + "".join(n[-1] for n in path[1:])
                unitigs.append(seq)

    return unitigs
