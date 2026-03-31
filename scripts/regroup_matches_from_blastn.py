#!/bin/python

from collections import defaultdict, deque

def build_graph(filename):
    graph = defaultdict(set)
    with open(filename) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) != 4:
                continue  # ignore invalid lines
            a = (parts[0], parts[1])
            b = (parts[2], parts[3])
            graph[a].add(b)
            graph[b].add(a)
    return graph

def find_connected_components(graph):
    visited = set()
    components = []

    for node in graph:
        if node not in visited:
            queue = deque([node])
            component = set()
            while queue:
                current = queue.popleft()
                if current in visited:
                    continue
                visited.add(current)
                component.add(current)
                queue.extend(graph[current] - visited)
            components.append(component)
    return components

def main(input_file, output_file):
    graph = build_graph(input_file)
    components = find_connected_components(graph)

    with open(output_file, "w") as out:
        for i, component in enumerate(components):
            for slice_name, species in sorted(component):
                out.write(f"{slice_name}\t{species}\tuce-{i}\n")

    print(f"{len(components)} unique UCEs detected. Results written in {output_file}")

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 3:
        print("Usage : python regroup_matches_from_blastn.py input.txt output.tsv")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])

