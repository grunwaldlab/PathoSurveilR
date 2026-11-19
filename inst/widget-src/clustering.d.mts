export type ClusterNode = {
  id: number;
  distance: number;
  left: ClusterNode | null;
  right: ClusterNode | null;
  leafIndex: number | undefined;
  leafStart: number;
  leafEnd: number;
};

export type Clustering = {
  root: ClusterNode | null;
  order: number[];
  maxDistance: number;
};

export function completeLinkage(
  matrix: Array<Array<number | null | undefined>>,
): Clustering;
