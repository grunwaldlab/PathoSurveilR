import test from "node:test";
import assert from "node:assert/strict";
import { completeLinkage } from "../../inst/widget-src/clustering.mjs";

function merges(root) {
  const result = [];
  const visit = (node) => {
    if (!node || node.leafIndex !== undefined) return;
    visit(node.left);
    visit(node.right);
    result.push(node);
  };
  visit(root);
  return result;
}

test("complete linkage merges the nearest Euclidean rows first", () => {
  const result = completeLinkage([
    [0, 0, 0, 0],
    [0, 1, 0, 0],
    [10, 10, 10, 10],
    [10, 11, 10, 10],
  ]);
  const treeMerges = merges(result.root);

  assert.equal(treeMerges.length, 3);
  assert.deepEqual(
    [treeMerges[0].left.leafIndex, treeMerges[0].right.leafIndex],
    [0, 1],
  );
  assert.equal(treeMerges[0].distance, 1);
  assert.ok(treeMerges[2].distance > treeMerges[0].distance);
  assert.deepEqual([...result.order].sort((a, b) => a - b), [0, 1, 2, 3]);
  assert.equal(result.root.leafStart, 0);
  assert.equal(result.root.leafEnd, 4);
});

test("cluster-to-cluster distances use the complete-linkage maximum", () => {
  const result = completeLinkage([[0], [2], [3]]);

  assert.equal(result.root.left.leafIndex, 0);
  assert.equal(result.root.right.left.leafIndex, 1);
  assert.equal(result.root.right.right.leafIndex, 2);
  assert.equal(result.root.distance, 3);
});

test("missing and non-finite values contribute zero", () => {
  const withMissing = completeLinkage([
    [null, 1, Number.NaN],
    [0, 1, 0],
    [5, 5, 5],
  ]);
  const withZeros = completeLinkage([
    [0, 1, 0],
    [0, 1, 0],
    [5, 5, 5],
  ]);

  assert.equal(withMissing.root.distance, withZeros.root.distance);
  assert.equal(merges(withMissing.root)[0].distance, 0);
});

test("equal distances have deterministic ordering", () => {
  const matrix = [
    [0, 0, 0],
    [0, 0, 0],
    [0, 0, 0],
  ];
  const first = completeLinkage(matrix);
  const second = completeLinkage(matrix);

  assert.deepEqual(first.order, second.order);
  assert.equal(first.root.distance, 0);
  assert.deepEqual([...first.order].sort((a, b) => a - b), [0, 1, 2]);
});

test("a one-row matrix returns a single leaf without a merge", () => {
  const result = completeLinkage([[1]]);

  assert.equal(result.root.leafIndex, 0);
  assert.deepEqual(result.order, [0]);
  assert.equal(result.maxDistance, 0);
});

test("an empty matrix returns an empty clustering", () => {
  assert.deepEqual(completeLinkage([]), {
    root: null,
    order: [],
    maxDistance: 0,
  });
});
