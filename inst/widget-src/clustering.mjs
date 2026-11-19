/**
 * @typedef {Object} ClusterNode
 * @property {number} id
 * @property {number} distance
 * @property {ClusterNode | null} left
 * @property {ClusterNode | null} right
 * @property {number | undefined} leafIndex
 * @property {number} leafStart
 * @property {number} leafEnd
 */

/**
 * @typedef {Object} Clustering
 * @property {ClusterNode | null} root
 * @property {number[]} order
 * @property {number} maxDistance
 */

/**
 * Complete-linkage hierarchical clustering with Euclidean leaf distances.
 * Missing and non-finite values contribute zero, matching the interactive
 * heatmap's previous clustering input treatment.
 *
 * @param {Array<Array<number | null | undefined>>} matrix
 * @returns {Clustering}
 */
export function completeLinkage(matrix) {
  const leafCount = matrix.length;
  if (!leafCount) return { root: null, order: [], maxDistance: 0 };

  /** @type {ClusterNode[]} */
  const nodes = Array.from({ length: leafCount }, (_, id) => ({
    id,
    distance: 0,
    left: null,
    right: null,
    leafIndex: id,
    leafStart: -1,
    leafEnd: -1,
  }));
  if (leafCount === 1) {
    nodes[0].leafStart = 0;
    nodes[0].leafEnd = 1;
    return { root: nodes[0], order: [0], maxDistance: 0 };
  }

  const distances = new Map();
  const key = (a, b) => a < b ? `${a},${b}` : `${b},${a}`;
  const setDistance = (a, b, value) => distances.set(key(a, b), value);
  const getDistance = (a, b) => distances.get(key(a, b));
  const heap = new MinHeap();
  const active = new Set(nodes.map((node) => node.id));

  for (let a = 0; a < leafCount; a += 1) {
    for (let b = a + 1; b < leafCount; b += 1) {
      let squaredDistance = 0;
      const aRow = matrix[a] ?? [];
      const bRow = matrix[b] ?? [];
      const columns = Math.max(aRow.length, bRow.length);
      for (let column = 0; column < columns; column += 1) {
        const rawA = aRow[column];
        const rawB = bRow[column];
        const valueA = typeof rawA === "number" && Number.isFinite(rawA) ? rawA : 0;
        const valueB = typeof rawB === "number" && Number.isFinite(rawB) ? rawB : 0;
        const difference = valueA - valueB;
        squaredDistance += difference * difference;
      }
      const distance = Math.sqrt(squaredDistance);
      setDistance(a, b, distance);
      heap.push({ a, b, distance });
    }
  }

  let nextId = leafCount;
  while (active.size > 1) {
    let closest = heap.pop();
    while (closest && (!active.has(closest.a) || !active.has(closest.b))) {
      closest = heap.pop();
    }
    if (!closest) throw new Error("Clustering stopped before reaching a root.");

    active.delete(closest.a);
    active.delete(closest.b);
    const merged = {
      id: nextId++,
      distance: closest.distance,
      left: nodes[closest.a],
      right: nodes[closest.b],
      leafIndex: undefined,
      leafStart: -1,
      leafEnd: -1,
    };
    nodes[merged.id] = merged;
    active.add(merged.id);

    for (const other of active) {
      if (other === merged.id) continue;
      const distance = Math.max(
        getDistance(closest.a, other),
        getDistance(closest.b, other),
      );
      setDistance(merged.id, other, distance);
      heap.push({ a: Math.min(merged.id, other), b: Math.max(merged.id, other), distance });
    }
  }

  const root = nodes[nextId - 1];
  const order = [];
  const assignLeafRanges = (node) => {
    node.leafStart = order.length;
    if (node.leafIndex !== undefined) {
      order.push(node.leafIndex);
    } else {
      assignLeafRanges(node.left);
      assignLeafRanges(node.right);
    }
    node.leafEnd = order.length;
  };
  assignLeafRanges(root);
  return { root, order, maxDistance: root.distance };
}

class MinHeap {
  constructor() {
    this.items = [];
  }

  push(candidate) {
    const items = this.items;
    items.push(candidate);
    let index = items.length - 1;
    while (index > 0) {
      const parent = Math.floor((index - 1) / 2);
      if (compareCandidates(items[parent], candidate) <= 0) break;
      items[index] = items[parent];
      index = parent;
    }
    items[index] = candidate;
  }

  pop() {
    const items = this.items;
    if (!items.length) return undefined;
    const first = items[0];
    const last = items.pop();
    if (items.length && last) {
      let index = 0;
      while (true) {
        const left = index * 2 + 1;
        const right = left + 1;
        if (left >= items.length) break;
        const child = right < items.length && compareCandidates(items[right], items[left]) < 0
          ? right
          : left;
        if (compareCandidates(last, items[child]) <= 0) break;
        items[index] = items[child];
        index = child;
      }
      items[index] = last;
    }
    return first;
  }
}

function compareCandidates(a, b) {
  return a.distance - b.distance || a.a - b.a || a.b - b.b;
}
