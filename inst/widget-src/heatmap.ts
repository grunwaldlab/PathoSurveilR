import {
  interpolateViridis,
  pointer,
  scaleSequential,
  select,
  zoom,
  zoomIdentity,
  type ZoomTransform,
} from "d3";
import {
  completeLinkage,
  type Clustering,
  type ClusterNode,
} from "./clustering.mjs";

type HeatmapPayload = {
  labels?: string[];
  display_labels?: string[];
  full_labels?: string[];
  matrix?: Array<Array<number | null>>;
  sample_ids?: string[];
  metric?: "ani" | "pocp" | "snp_distance";
  metric_label?: string;
};

type Cell = {
  row: number;
  column: number;
  value: number | null;
};

type HeatmapInstance = {
  dispose: () => void;
  resize: () => void;
  getZoomFocus: () => ZoomFocus;
  getClustering: () => Clustering | undefined;
};

type ZoomFocus = {
  scale: number;
  centerX: number;
  centerY: number;
};

let chartSequence = 0;

const mountHeatmap = (
  container: HTMLElement,
  payload: unknown,
  {
    signal,
    zoomFocus,
    clustering: cachedClustering,
  }: { signal: AbortSignal; zoomFocus?: ZoomFocus; clustering?: Clustering },
): HeatmapInstance => {
  if (signal.aborted) {
    return {
      dispose: () => {},
      resize: () => {},
      getZoomFocus: () => ({ scale: 1, centerX: 0.5, centerY: 0.5 }),
      getClustering: () => undefined,
    };
  }

  const output = (payload ?? {}) as HeatmapPayload;
  const isDistance = output.metric === "snp_distance";
  const metricName = isDistance ? "SNP distance" : output.metric === "pocp" ? "POCP" : "ANI";
  const labels = Array.isArray(output.labels) ? output.labels : [];
  const matrix = Array.isArray(output.matrix) ? output.matrix : [];
  const validMatrix = matrix.length === labels.length
    && matrix.every((row) => Array.isArray(row) && row.length >= labels.length);

  if (!labels.length || !validMatrix) {
    const dispose = showEmptyState(container, signal, output.metric_label);
    return {
      dispose,
      resize: () => {},
      getZoomFocus: () => ({ scale: 1, centerX: 0.5, centerY: 0.5 }),
      getClustering: () => undefined,
    };
  }

  const displayLabels = labels.map((label, index) => {
    const displayLabel = output.display_labels?.[index];
    return typeof displayLabel === "string" ? displayLabel : label;
  });
  const fullLabels = labels.map((label, index) => {
    const fullLabel = output.full_labels?.[index];
    return typeof fullLabel === "string" ? fullLabel : label;
  });
  const clustering = cachedClustering ?? completeLinkage(matrix);
  const order = clustering.order;
  const plotLabels = order.map((index) => labels[index]);
  const plotDisplayLabels = order.map((index) => displayLabels[index]);
  const plotFullLabels = order.map((index) => fullLabels[index]);
  const plotMatrix = order.map((row) => order.map((column) => {
    const value = matrix[row]?.[column];
    return typeof value === "number" && Number.isFinite(value) ? value : null;
  }));
  const cells: Cell[] = plotMatrix.flatMap((row, rowIndex) =>
    row.map((value, columnIndex) => ({ row: rowIndex, column: columnIndex, value })),
  );
  const values = cells.flatMap((cell) => cell.value === null ? [] : [cell.value]);

  if (!values.length) {
    const dispose = showEmptyState(container, signal, output.metric_label);
    return {
      dispose,
      resize: () => {},
      getZoomFocus: () => ({ scale: 1, centerX: 0.5, centerY: 0.5 }),
      getClustering: () => clustering,
    };
  }

  const sampleIds = new Set(Array.isArray(output.sample_ids) ? output.sample_ids : []);
  const minValue = Math.min(...values);
  const maxValue = Math.max(...values);
  const colorScale = scaleSequential(interpolateViridis).domain(
    minValue === maxValue
      ? isDistance ? [maxValue + 0.5, minValue - 0.5] : [minValue - 0.5, maxValue + 0.5]
      : isDistance ? [maxValue, minValue] : [minValue, maxValue],
  );
  const missingColor = "#b8b8b8";
  const toolbarHost = container.parentElement ?? container;
  const previousToolbarHostPosition = toolbarHost.style.position;
  const positionToolbarHost = getComputedStyle(toolbarHost).position === "static";
  if (positionToolbarHost) toolbarHost.style.position = "relative";
  const previousContainerDisplay = container.style.display;
  const previousContainerAlignItems = container.style.alignItems;
  const previousContainerJustifyContent = container.style.justifyContent;
  container.style.display = "flex";
  container.style.alignItems = "center";
  container.style.justifyContent = "center";
  const root = document.createElement("div");
  root.className = "ani-heatmap";
  Object.assign(root.style, {
    position: "relative",
    display: "flex",
    flexDirection: "column",
    alignItems: "center",
    width: "max-content",
    maxWidth: "100%",
    color: "#202124",
    font: "13px/1.35 system-ui, sans-serif",
  });
  const toolbar = document.createElement("div");
  toolbar.className = "ani-heatmap-controls";
  Object.assign(toolbar.style, {
    position: "absolute",
    zIndex: "2",
    top: "0",
    display: "flex",
    flexDirection: "column",
    alignItems: "flex-start",
    gap: "5px",
  });
  const resetButton = document.createElement("button");
  resetButton.type = "button";
  resetButton.className = "ani-heatmap-reset";
  resetButton.setAttribute("aria-label", "Reset zoom and pan");
  resetButton.title = "Reset zoom and pan";
  Object.assign(resetButton.style, {
    display: "grid",
    placeItems: "center",
    width: "34px",
    height: "34px",
    padding: "5px",
    border: "1px solid transparent",
    borderRadius: "4px",
    background: "transparent",
    color: "#9ca3af",
    cursor: "pointer",
  });
  const resetIcon = document.createElementNS("http://www.w3.org/2000/svg", "svg");
  resetIcon.setAttribute("viewBox", "0 0 24 24");
  resetIcon.setAttribute("aria-hidden", "true");
  Object.assign(resetIcon.style, { width: "22px", height: "22px" });
  const resetPath = document.createElementNS("http://www.w3.org/2000/svg", "path");
  resetPath.setAttribute("d", "M20 11a8 8 0 1 1-2.34-5.66L20 8m0-4v4h-4");
  resetPath.setAttribute("fill", "none");
  resetPath.setAttribute("stroke", "currentColor");
  resetPath.setAttribute("stroke-width", "1.8");
  resetPath.setAttribute("stroke-linecap", "round");
  resetPath.setAttribute("stroke-linejoin", "round");
  resetIcon.append(resetPath);
  resetButton.append(resetIcon);
  const zoomToggle = document.createElement("button");
  zoomToggle.type = "button";
  zoomToggle.className = "ani-heatmap-zoom-toggle";
  zoomToggle.setAttribute("aria-label", "Enable zoom and pan");
  zoomToggle.setAttribute("aria-pressed", "false");
  zoomToggle.title = "Enable zoom and pan";
  Object.assign(zoomToggle.style, {
    display: "grid",
    placeItems: "center",
    width: "34px",
    height: "34px",
    padding: "5px",
    border: "1px solid transparent",
    borderRadius: "4px",
    background: "transparent",
    color: "#9ca3af",
    opacity: "0.55",
    cursor: "pointer",
  });
  const zoomIcon = document.createElementNS("http://www.w3.org/2000/svg", "svg");
  zoomIcon.setAttribute("viewBox", "0 0 24 24");
  zoomIcon.setAttribute("aria-hidden", "true");
  Object.assign(zoomIcon.style, { width: "22px", height: "22px" });
  const iconCircle = document.createElementNS("http://www.w3.org/2000/svg", "circle");
  iconCircle.setAttribute("cx", "9.5");
  iconCircle.setAttribute("cy", "9.5");
  iconCircle.setAttribute("r", "6.5");
  const iconPath = document.createElementNS("http://www.w3.org/2000/svg", "path");
  iconPath.setAttribute("d", "m14.5 14.5 6 6M9.5 6.5v6M6.5 9.5h6");
  for (const shape of [iconCircle, iconPath]) {
    shape.setAttribute("fill", "none");
    shape.setAttribute("stroke", "currentColor");
    shape.setAttribute("stroke-width", "1.8");
    shape.setAttribute("stroke-linecap", "round");
    shape.setAttribute("stroke-linejoin", "round");
    zoomIcon.append(shape);
  }
  zoomToggle.append(zoomIcon);
  const zoomSlider = document.createElement("input");
  zoomSlider.type = "range";
  zoomSlider.className = "ani-heatmap-zoom-slider";
  zoomSlider.min = "0";
  zoomSlider.max = "1000";
  zoomSlider.step = "1";
  zoomSlider.value = "0";
  zoomSlider.setAttribute("aria-label", "Heatmap zoom level");
  zoomSlider.title = "Adjust heatmap zoom";
  Object.assign(zoomSlider.style, {
    width: "22px",
    height: "112px",
    margin: "0 6px",
    writingMode: "vertical-lr",
    direction: "rtl",
    accentColor: "#4f46e5",
    cursor: "pointer",
  });
  toolbar.append(resetButton, zoomToggle, zoomSlider);

  const legend = document.createElement("div");
  legend.className = "ani-heatmap-legend";
  legend.setAttribute(
    "aria-label",
    isDistance ? "SNP distance color legend" : `${metricName} similarity color legend`,
  );
  Object.assign(legend.style, {
    display: "flex",
    alignItems: "center",
    justifyContent: "center",
    flexWrap: "wrap",
    gap: "14px",
    minHeight: "42px",
    width: "max-content",
    maxWidth: "calc(100% - 16px)",
    boxSizing: "border-box",
    padding: "4px 8px",
  });
  const scaleLegend = document.createElement("div");
  Object.assign(scaleLegend.style, {
    display: "flex",
    flexDirection: "column",
    alignItems: "center",
    gap: "3px",
  });
  const legendDetails = document.createElement("div");
  Object.assign(legendDetails.style, {
    display: "flex",
    alignItems: "center",
    justifyContent: "center",
    flexWrap: "nowrap",
    gap: "14px",
  });
  const scaleTitle = document.createElement("span");
  scaleTitle.textContent = output.metric_label ?? (isDistance ? "SNP differences (count)" : "Similarity (%)");
  const scaleRange = document.createElement("div");
  Object.assign(scaleRange.style, {
    display: "grid",
    gridTemplateColumns: "auto minmax(48px, 220px) auto",
    alignItems: "center",
    gap: "7px",
  });
  const scaleMinimum = document.createElement("span");
  scaleMinimum.textContent = formatValue(minValue, isDistance);
  const gradient = document.createElement("span");
  gradient.setAttribute("aria-hidden", "true");
  Object.assign(gradient.style, {
    height: "12px",
    border: "1px solid #777",
    background: isDistance
      ? `linear-gradient(to right, ${interpolateViridis(1)}, ${interpolateViridis(0.75)}, ${interpolateViridis(0.5)}, ${interpolateViridis(0.25)}, ${interpolateViridis(0)})`
      : `linear-gradient(to right, ${interpolateViridis(0)}, ${interpolateViridis(0.25)}, ${interpolateViridis(0.5)}, ${interpolateViridis(0.75)}, ${interpolateViridis(1)})`,
  });
  const scaleMaximum = document.createElement("span");
  scaleMaximum.textContent = formatValue(maxValue, isDistance);
  scaleRange.append(scaleMinimum, gradient, scaleMaximum);
  legendDetails.append(scaleRange);
  const missingLegend = document.createElement("div");
  Object.assign(missingLegend.style, { display: "flex", alignItems: "center", gap: "6px" });
  const missingSwatch = document.createElement("span");
  Object.assign(missingSwatch.style, {
    display: "inline-block",
    width: "13px",
    height: "13px",
    border: "1px solid #777",
    background: missingColor,
  });
  const missingLabel = document.createElement("span");
  missingLabel.textContent = `Missing ${metricName}`;
  missingLegend.append(missingSwatch, missingLabel);
  legendDetails.append(missingLegend);
  scaleLegend.append(scaleTitle, legendDetails);
  legend.append(scaleLegend);
  root.append(legend);
  container.replaceChildren(root);
  toolbarHost.append(toolbar);

  const maxWidgetWidth = cssDimension(container, "--widget-max-width", 1600);
  const maxWidgetHeight = cssDimension(container, "--widget-max-height", 1440);
  const legendHeight = 60;
  const measurementSvg = document.createElementNS("http://www.w3.org/2000/svg", "svg");
  Object.assign(measurementSvg.style, {
    position: "absolute",
    left: "-10000px",
    top: "0",
    width: "1px",
    height: "1px",
    overflow: "visible",
    visibility: "hidden",
    pointerEvents: "none",
  });
  const measurementText = document.createElementNS("http://www.w3.org/2000/svg", "text");
  measurementSvg.append(measurementText);
  root.append(measurementSvg);
  const measureLabels = (fontSize: number) => Math.max(...plotDisplayLabels.map((label, index) => {
    measurementText.setAttribute("font-size", String(fontSize));
    measurementText.setAttribute("font-weight", sampleIds.has(plotLabels[index]) ? "700" : "400");
    measurementText.textContent = label;
    return measurementText.getComputedTextLength();
  }));
  const desiredPlotSide = Math.max(
    1,
    Math.min(maxWidgetWidth - 16, maxWidgetHeight - legendHeight),
  );
  const calculateLayout = (showTrees: boolean) => {
    let plotSide = desiredPlotSide;
    let rowLabelSpace = 0;
    let columnLabelSpace = 0;
    let rowTreeSpace = 0;
    let columnTreeSpace = 0;
    for (let attempt = 0; attempt < 5; attempt += 1) {
      const fontSize = Math.min(14, plotSide / labels.length);
      const widestLabel = measureLabels(fontSize);
      rowLabelSpace = Math.min(maxWidgetWidth * 0.42, widestLabel + 8);
      columnLabelSpace = Math.min(maxWidgetHeight * 0.24, widestLabel + 8);
      rowTreeSpace = showTrees ? rowLabelSpace : 0;
      columnTreeSpace = showTrees ? columnLabelSpace : 0;
      plotSide = Math.max(1, Math.min(
        desiredPlotSide,
        maxWidgetWidth - rowLabelSpace - rowTreeSpace - 16,
        maxWidgetHeight - columnLabelSpace - columnTreeSpace - legendHeight,
      ));
    }
    return { plotSide, rowLabelSpace, columnLabelSpace, rowTreeSpace, columnTreeSpace };
  };
  const noTreeLayout = calculateLayout(false);
  const candidateTreeLayout = labels.length > 1 ? calculateLayout(true) : noTreeLayout;
  const minimumTreePlotSide = 300;
  const showDendrograms = labels.length > 1 && candidateTreeLayout.plotSide >= minimumTreePlotSide;
  const layout = showDendrograms ? candidateTreeLayout : noTreeLayout;
  const { plotSide, rowLabelSpace, columnLabelSpace, rowTreeSpace, columnTreeSpace } = layout;
  measurementSvg.remove();
  const plotLeft = rowTreeSpace + rowLabelSpace;
  const plotTop = 0;
  const svgWidth = plotLeft + plotSide + 16;
  const treeTop = plotTop + plotSide + columnLabelSpace;
  const svgHeight = treeTop + columnTreeSpace;
  root.style.width = `${svgWidth}px`;
  toolbar.style.left = "0";
  toolbar.style.top = "0";

  const svg = document.createElementNS("http://www.w3.org/2000/svg", "svg");
  svg.setAttribute("role", "img");
  svg.setAttribute(
    "aria-label",
    isDistance ? "Pairwise SNP distance heatmap" : `Pairwise ${metricName} similarity heatmap`,
  );
  svg.setAttribute("viewBox", `0 0 ${svgWidth} ${svgHeight}`);
  svg.setAttribute("width", String(svgWidth));
  svg.setAttribute("height", String(svgHeight));
  svg.setAttribute("data-plot-side", String(plotSide));
  svg.setAttribute("data-plot-left", String(plotLeft));
  svg.setAttribute("data-plot-top", String(plotTop));
  svg.setAttribute("data-row-label-space", String(rowLabelSpace));
  svg.setAttribute("data-column-label-space", String(columnLabelSpace));
  svg.setAttribute("data-row-tree-space", String(rowTreeSpace));
  svg.setAttribute("data-column-tree-space", String(columnTreeSpace));
  svg.setAttribute("data-dendrograms", String(showDendrograms));
  Object.assign(svg.style, {
    display: "block",
    flex: "0 0 auto",
    width: `${svgWidth}px`,
    height: "auto",
    maxWidth: "100%",
    maxHeight: `calc(var(--widget-max-height, ${maxWidgetHeight}px) - ${legendHeight}px)`,
    overflow: "hidden",
    touchAction: "auto",
  });
  root.insertBefore(svg, legend);
  const tooltip = document.createElement("div");
  tooltip.className = "ani-heatmap-tooltip";
  tooltip.setAttribute("role", "tooltip");
  tooltip.setAttribute("aria-hidden", "true");
  Object.assign(tooltip.style, {
    position: "absolute",
    zIndex: "3",
    display: "none",
    maxWidth: "min(360px, calc(100% - 16px))",
    padding: "6px 8px",
    border: "1px solid #737373",
    borderRadius: "4px",
    background: "rgba(255,255,255,0.97)",
    color: "#171717",
    boxShadow: "0 2px 6px rgba(0,0,0,0.2)",
    font: "12px/1.35 system-ui, sans-serif",
    whiteSpace: "pre-line",
    overflowWrap: "anywhere",
    pointerEvents: "none",
  });
  root.append(tooltip);

  const svgSelection = select(svg);
  const clipId = `ani-heatmap-clip-${++chartSequence}`;
  const defs = svgSelection.append("defs");
  const clipRect = defs
    .append("clipPath")
    .attr("id", clipId)
    .append("rect");
  const rowTreeClipId = `ani-heatmap-row-tree-clip-${chartSequence}`;
  const columnTreeClipId = `ani-heatmap-column-tree-clip-${chartSequence}`;
  if (showDendrograms) {
    defs.append("clipPath")
      .attr("id", rowTreeClipId)
      .append("rect")
      .attr("x", 0)
      .attr("y", plotTop)
      .attr("width", rowTreeSpace)
      .attr("height", plotSide);
    defs.append("clipPath")
      .attr("id", columnTreeClipId)
      .append("rect")
      .attr("x", plotLeft)
      .attr("y", treeTop)
      .attr("width", plotSide)
      .attr("height", columnTreeSpace);
  }
  const matrixViewport = svgSelection.append("g").attr("clip-path", `url(#${clipId})`);
  const grid = matrixViewport.append("g").attr("class", "ani-heatmap-grid");
  const valueLabels = matrixViewport.append("g").attr("class", "ani-heatmap-values");
  const selectionLines = matrixViewport.append("g").attr("class", "ani-heatmap-selection-lines");
  const interaction = svgSelection.append("rect")
    .attr("class", "ani-heatmap-interaction")
    .attr("fill", "transparent")
    .attr("pointer-events", "all");
  const rowLabels = svgSelection.append("g").attr("class", "ani-heatmap-row-labels");
  const columnLabels = svgSelection.append("g").attr("class", "ani-heatmap-column-labels");
  const rowTreePath = showDendrograms
    ? svgSelection.append("g")
      .attr("class", "ani-heatmap-row-dendrogram")
      .attr("clip-path", `url(#${rowTreeClipId})`)
      .append("path")
      .attr("fill", "none")
      .attr("stroke", "#555")
      .attr("stroke-width", 1.25)
      .attr("vector-effect", "non-scaling-stroke")
      .attr("pointer-events", "none")
    : null;
  const columnTreePath = showDendrograms
    ? svgSelection.append("g")
      .attr("class", "ani-heatmap-column-dendrogram")
      .attr("clip-path", `url(#${columnTreeClipId})`)
      .append("path")
      .attr("fill", "none")
      .attr("stroke", "#555")
      .attr("stroke-width", 1.25)
      .attr("vector-effect", "non-scaling-stroke")
      .attr("pointer-events", "none")
    : null;
  const valueMeasurement = svgSelection.append("text")
    .attr("x", -10000)
    .attr("y", -10000)
    .attr("font-weight", "600")
    .attr("visibility", "hidden")
    .attr("pointer-events", "none");
  const selectionRect = svgSelection.append("rect")
    .attr("class", "ani-heatmap-selection")
    .attr("fill", "rgba(255,255,255,0.18)")
    .attr("stroke", "#111")
    .attr("stroke-dasharray", "4 3")
    .attr("stroke-width", 1.5)
    .attr("pointer-events", "none")
    .attr("display", "none");
  clipRect.attr("x", plotLeft).attr("y", plotTop).attr("width", plotSide).attr("height", plotSide);
  interaction
    .attr("x", plotLeft)
    .attr("y", plotTop)
    .attr("width", plotSide)
    .attr("height", plotSide);

  const plotCells = grid.selectAll<SVGRectElement, Cell>("rect.ani-heatmap-cell")
    .data(cells)
    .join("rect")
    .attr("class", "ani-heatmap-cell")
    .attr("data-row", (cell) => cell.row)
    .attr("data-column", (cell) => cell.column)
    .attr("x", (cell) => plotLeft + cell.column * (plotSide / labels.length))
    .attr("y", (cell) => plotTop + cell.row * (plotSide / labels.length))
    .attr("width", plotSide / labels.length)
    .attr("height", plotSide / labels.length)
    .attr("fill", (cell) => cell.value === null ? missingColor : colorScale(cell.value))
    .attr("shape-rendering", "crispEdges");
  plotCells.append("title");
  plotCells.select("title").text((cell) =>
    `${plotFullLabels[cell.row]} versus ${plotFullLabels[cell.column]}\n${output.metric_label ?? (isDistance ? "SNP differences" : "Similarity")}: ${cell.value === null ? "Missing" : formatValue(cell.value, isDistance)}`,
  );

  const cellSize = plotSide / labels.length;
  let transform: ZoomTransform = zoomFocus
    ? zoomIdentity.translate(
      plotLeft + plotSide / 2 - (plotLeft + zoomFocus.centerX * plotSide) * zoomFocus.scale,
      plotTop + plotSide / 2 - (plotTop + zoomFocus.centerY * plotSide) * zoomFocus.scale,
    ).scale(zoomFocus.scale)
    : zoomIdentity;
  let disposed = false;
  let zoomEnabled = true;
  let dragStart: [number, number] | undefined;
  let suppressClick = false;
  let selectedRow: number | null = null;
  let selectedColumn: number | null = null;
  const maxZoomScale = Math.max(2, labels.length);
  const hideCellTooltip = () => {
    tooltip.style.display = "none";
    tooltip.setAttribute("aria-hidden", "true");
  };
  const cellAt = (event: MouseEvent | PointerEvent): Cell | undefined => {
    const [screenX, screenY] = pointer(event, svg);
    if (
      screenX < plotLeft
      || screenX > plotLeft + plotSide
      || screenY < plotTop
      || screenY > plotTop + plotSide
    ) {
      return;
    }
    const row = Math.floor((transform.invertY(screenY) - plotTop) / cellSize);
    const column = Math.floor((transform.invertX(screenX) - plotLeft) / cellSize);
    return row >= 0 && row < labels.length && column >= 0 && column < labels.length
      ? cells[row * labels.length + column]
      : undefined;
  };
  const showCellTooltip = (event: PointerEvent) => {
    if (event.pointerType !== "mouse" || dragStart) {
      hideCellTooltip();
      return;
    }
    const cell = cellAt(event);
    if (!cell) {
      hideCellTooltip();
      return;
    }
    tooltip.textContent = `Row: ${plotFullLabels[cell.row]}\nColumn: ${plotFullLabels[cell.column]}\n${output.metric_label ?? (isDistance ? "SNP differences" : "Similarity")}: ${cell.value === null ? "Missing" : formatValue(cell.value, isDistance)}`;
    tooltip.style.display = "block";
    tooltip.setAttribute("aria-hidden", "false");
    const rootBounds = root.getBoundingClientRect();
    const tooltipBounds = tooltip.getBoundingClientRect();
    const left = Math.max(4, Math.min(event.clientX - rootBounds.left + 12, rootBounds.width - tooltipBounds.width - 4));
    const top = Math.max(4, Math.min(event.clientY - rootBounds.top + 12, rootBounds.height - tooltipBounds.height - 4));
    tooltip.style.left = `${left}px`;
    tooltip.style.top = `${top}px`;
  };

  const zoomBehavior = zoom<SVGRectElement, unknown>()
    .extent([[plotLeft, plotTop], [plotLeft + plotSide, plotTop + plotSide]])
    .translateExtent([[plotLeft, plotTop], [plotLeft + plotSide, plotTop + plotSide]])
    .scaleExtent([1, maxZoomScale])
    .filter((event) => {
      if (event.type === "wheel") return zoomEnabled;
      if (event.type.startsWith("touch")) return zoomEnabled || transform.k > 1.000001;
      if (event.type === "dblclick") return false;
      return event.type === "mousedown" && event.button === 0 && transform.k > 1.000001;
    })
    .on("zoom", (event) => {
      transform = event.transform;
      zoomSlider.value = String(Math.round(
        (Math.log(Math.max(1, transform.k)) / Math.log(maxZoomScale)) * 1000,
      ));
      hideCellTooltip();
      renderZoom();
    });
  interaction.call(zoomBehavior);
  const dendrogramNodes = showDendrograms ? internalNodes(clustering.root) : [];

  const renderZoom = () => {
    if (disposed) return;
    svg.setAttribute("data-selected-row", selectedRow === null ? "" : String(selectedRow));
    svg.setAttribute("data-selected-column", selectedColumn === null ? "" : String(selectedColumn));
    grid.attr("transform", transform.toString());
    valueLabels.attr("transform", transform.toString());
    selectionLines.attr("transform", transform.toString());
    if (rowTreePath && columnTreePath && clustering.root) {
      rowTreePath.attr("d", rowDendrogramPath(
        dendrogramNodes,
        labels.length,
        clustering.maxDistance,
        transform,
        plotTop,
        cellSize,
        rowTreeSpace,
      ));
      columnTreePath.attr("d", columnDendrogramPath(
        dendrogramNodes,
        labels.length,
        clustering.maxDistance,
        transform,
        plotLeft,
        plotTop + plotSide,
        columnLabelSpace,
        columnTreeSpace,
        cellSize,
      ));
    }
    const visibleCellSize = cellSize * transform.k;
    const fontSize = Math.min(14, visibleCellSize);
    const displayScale = svg.getBoundingClientRect().width / svgWidth;
    const displayedCellSize = visibleCellSize * displayScale;
    const hasSelection = selectedRow !== null || selectedColumn !== null;
    plotCells.attr("fill", (cell) => cell.value === null ? missingColor : colorScale(cell.value));
    interaction.attr("cursor", transform.k > 1.000001 ? "grab" : zoomEnabled ? "crosshair" : "default");

    const visibleRows = plotLabels.flatMap((id, index) => {
      const center = transform.applyY(plotTop + (index + 0.5) * cellSize);
      return center >= plotTop && center <= plotTop + plotSide
        ? [{ id, label: plotDisplayLabels[index], fullLabel: plotFullLabels[index], index, center }]
        : [];
    });
    rowLabels.selectAll<SVGTextElement, (typeof visibleRows)[number]>("text")
      .data(visibleCellSize >= 4 ? visibleRows : [], (datum) => String(datum.index))
      .join(
        (enter) => {
          const text = enter.append("text").attr("class", "ani-heatmap-label");
          text.append("title");
          return text;
        },
        (update) => update,
        (exit) => exit.remove(),
      )
      .attr("data-label-id", (datum) => datum.id)
      .attr("data-row-index", (datum) => datum.index)
      .attr("data-axis", "row")
      .attr("class", (datum) => `ani-heatmap-label${sampleIds.has(datum.id) ? " sample" : " reference"}`)
      .attr("x", plotLeft - 4)
      .attr("y", (datum) => datum.center + fontSize * 0.35)
      .attr("text-anchor", "end")
      .attr("font-size", fontSize)
      .attr("font-weight", (datum) => sampleIds.has(datum.id) ? "700" : "400")
      .attr("fill", "#202124")
      .attr("role", "button")
      .attr("tabindex", 0)
      .attr("aria-label", (datum) => `Select row ${datum.fullLabel}`)
      .attr("cursor", "pointer")
      .on("click.ani-selection", (event, datum) => {
        event.stopPropagation();
        selectRow(datum.index);
      })
      .on("keydown.ani-selection", (event, datum) => {
        if (event.key !== "Enter" && event.key !== " ") return;
        event.preventDefault();
        selectRow(datum.index);
      })
      .text((datum) => truncateLabel(datum.label, rowLabelSpace - 8, fontSize));
    rowLabels.selectAll<SVGTextElement, (typeof visibleRows)[number]>("text")
      .append("title")
      .text((datum) => datum.fullLabel);

    const visibleColumns = plotLabels.flatMap((id, index) => {
      const center = transform.applyX(plotLeft + (index + 0.5) * cellSize);
      return center >= plotLeft && center <= plotLeft + plotSide
        ? [{ id, label: plotDisplayLabels[index], fullLabel: plotFullLabels[index], index, center }]
        : [];
    });
    columnLabels.selectAll<SVGTextElement, (typeof visibleColumns)[number]>("text")
      .data(visibleCellSize >= 4 ? visibleColumns : [], (datum) => String(datum.index))
      .join(
        (enter) => {
          const text = enter.append("text").attr("class", "ani-heatmap-label");
          text.append("title");
          return text;
        },
        (update) => update,
        (exit) => exit.remove(),
      )
      .attr("data-label-id", (datum) => datum.id)
      .attr("data-column-index", (datum) => datum.index)
      .attr("data-axis", "column")
      .attr("class", (datum) => `ani-heatmap-label${sampleIds.has(datum.id) ? " sample" : " reference"}`)
      .attr("transform", (datum) => `translate(${datum.center},${plotTop + plotSide + 4}) rotate(-90)`)
      .attr("x", 0)
      .attr("dy", fontSize * 0.35)
      .attr("text-anchor", "end")
      .attr("font-size", fontSize)
      .attr("font-weight", (datum) => sampleIds.has(datum.id) ? "700" : "400")
      .attr("fill", "#202124")
      .attr("role", "button")
      .attr("tabindex", 0)
      .attr("aria-label", (datum) => `Select column ${datum.fullLabel}`)
      .attr("cursor", "pointer")
      .on("click.ani-selection", (event, datum) => {
        event.stopPropagation();
        selectColumn(datum.index);
      })
      .on("keydown.ani-selection", (event, datum) => {
        if (event.key !== "Enter" && event.key !== " ") return;
        event.preventDefault();
        selectColumn(datum.index);
      })
      .text((datum) => truncateLabel(datum.label, columnLabelSpace - 8, fontSize));
    columnLabels.selectAll<SVGTextElement, (typeof visibleColumns)[number]>("text")
      .append("title")
      .text((datum) => datum.fullLabel);

    const selectionLineData: Array<{
      key: string;
      x1: number;
      y1: number;
      x2: number;
      y2: number;
    }> = [];
    if (selectedRow !== null) {
      const y = plotTop + selectedRow * cellSize;
      selectionLineData.push(
        { key: "row-top", x1: plotLeft, y1: y, x2: plotLeft + plotSide, y2: y },
        { key: "row-bottom", x1: plotLeft, y1: y + cellSize, x2: plotLeft + plotSide, y2: y + cellSize },
      );
    }
    if (selectedColumn !== null) {
      const x = plotLeft + selectedColumn * cellSize;
      selectionLineData.push(
        { key: "column-left", x1: x, y1: plotTop, x2: x, y2: plotTop + plotSide },
        { key: "column-right", x1: x + cellSize, y1: plotTop, x2: x + cellSize, y2: plotTop + plotSide },
      );
    }
    selectionLines.selectAll<SVGLineElement, (typeof selectionLineData)[number]>("line")
      .data(selectionLineData, (datum) => datum.key)
      .join("line")
      .attr("class", (datum) => `ani-heatmap-selection-line ${datum.key}`)
      .attr("x1", (datum) => datum.x1)
      .attr("y1", (datum) => datum.y1)
      .attr("x2", (datum) => datum.x2)
      .attr("y2", (datum) => datum.y2)
      .attr("stroke", "#111")
      .attr("stroke-width", 1.5)
      .attr("stroke-dasharray", "2 2")
      .attr("vector-effect", "non-scaling-stroke")
      .attr("pointer-events", "none");

    valueMeasurement
      .attr("font-size", 1)
      .attr("font-weight", "600")
      .text("0000");
    const valueFontWidth = valueMeasurement.node()?.getComputedTextLength() ?? Infinity;
    const availableValueWidth = displayedCellSize - 4;
    const transformScale = displayScale * transform.k;
    const screenFontSize = availableValueWidth / valueFontWidth;
    const valueFontSize = screenFontSize / transformScale;
    svg.setAttribute("data-value-font-size", String(screenFontSize));
    const valueLabelsToShow = hasSelection && screenFontSize >= 4
      ? cells.flatMap((cell) => {
        if (
          cell.value === null
          || (cell.row !== selectedRow && cell.column !== selectedColumn)
        ) {
          return [];
        }
        return [{ ...cell, label: formatCellValue(cell.value, isDistance) }];
      })
      : [];
    valueLabels.selectAll<SVGTextElement, (typeof valueLabelsToShow)[number]>("text")
      .data(valueLabelsToShow, (datum) => `${datum.row}:${datum.column}`)
      .join("text")
      .attr("class", "ani-heatmap-value")
      .attr("data-row", (datum) => datum.row)
      .attr("data-column", (datum) => datum.column)
      .attr("x", (datum) => plotLeft + (datum.column + 0.5) * cellSize)
      .attr("y", (datum) => plotTop + (datum.row + 0.5) * cellSize)
      .attr("text-anchor", "middle")
      .attr("dominant-baseline", "central")
      .attr("font-size", valueFontSize)
      .attr("font-weight", "600")
      .attr("fill", "#111")
      .attr("pointer-events", "none")
      .text((datum) => datum.label);
  };

  const selectRow = (row: number) => {
    hideCellTooltip();
    const selected = selectedRow === row && selectedColumn === row;
    selectedRow = selected ? null : row;
    selectedColumn = selected ? null : row;
    renderZoom();
  };
  const selectColumn = (column: number) => {
    hideCellTooltip();
    const selected = selectedRow === column && selectedColumn === column;
    selectedRow = selected ? null : column;
    selectedColumn = selected ? null : column;
    renderZoom();
  };
  const selectCell = (row: number, column: number) => {
    hideCellTooltip();
    const selected = selectedRow === row && selectedColumn === column;
    selectedRow = selected ? null : row;
    selectedColumn = selected ? null : column;
    renderZoom();
  };

  const selectionPoint = (event: PointerEvent): [number, number] => {
    const [x, y] = pointer(event, svg);
    return [
      Math.max(plotLeft, Math.min(plotLeft + plotSide, x)),
      Math.max(plotTop, Math.min(plotTop + plotSide, y)),
    ];
  };
  const onPointerDown = (event: PointerEvent) => {
    if (event.pointerType !== "mouse" || event.button !== 0 || transform.k > 1.000001) return;
    hideCellTooltip();
    suppressClick = false;
    event.preventDefault();
    dragStart = selectionPoint(event);
    selectionRect.attr("display", null).attr("x", dragStart[0]).attr("y", dragStart[1]).attr("width", 0).attr("height", 0);
    interaction.node()?.setPointerCapture(event.pointerId);
  };
  const onPointerMove = (event: PointerEvent) => {
    if (!dragStart) return;
    const [x, y] = selectionPoint(event);
    const directionX = x >= dragStart[0] ? 1 : -1;
    const directionY = y >= dragStart[1] ? 1 : -1;
    const maxX = directionX > 0 ? plotLeft + plotSide - dragStart[0] : dragStart[0] - plotLeft;
    const maxY = directionY > 0 ? plotTop + plotSide - dragStart[1] : dragStart[1] - plotTop;
    const side = Math.min(Math.abs(x - dragStart[0]), Math.abs(y - dragStart[1]), maxX, maxY);
    const left = directionX > 0 ? dragStart[0] : dragStart[0] - side;
    const top = directionY > 0 ? dragStart[1] : dragStart[1] - side;
    selectionRect.attr("x", left).attr("y", top).attr("width", side).attr("height", side);
  };
  const onPointerUp = (event: PointerEvent) => {
    if (!dragStart) return;
    onPointerMove(event);
    const selection = selectionRect.node();
    const x = Number(selection?.getAttribute("x"));
    const y = Number(selection?.getAttribute("y"));
    const side = Number(selection?.getAttribute("width"));
    dragStart = undefined;
    selectionRect.attr("display", "none");
    if (side < 5) return;
    suppressClick = true;
    const scale = Math.min(labels.length, plotSide / side);
    interaction.call(
      zoomBehavior.transform,
      zoomIdentity.translate(plotLeft - x * scale, plotTop - y * scale).scale(scale),
    );
  };
  interaction
    .on("pointerdown.ani-selection", onPointerDown)
    .on("pointermove.ani-selection", onPointerMove)
    .on("pointermove.ani-tooltip", showCellTooltip)
    .on("pointerup.ani-selection", onPointerUp)
    .on("pointerleave.ani-tooltip", hideCellTooltip)
    .on("pointercancel.ani-selection", () => {
      dragStart = undefined;
      suppressClick = false;
      selectionRect.attr("display", "none");
    });
  interaction.on("click.ani-selection", (event: MouseEvent) => {
    if (suppressClick) {
      suppressClick = false;
      event.preventDefault();
      return;
    }
    const cell = cellAt(event);
    if (cell) selectCell(cell.row, cell.column);
  });
  const updateZoomToggle = () => {
    const label = zoomEnabled ? "Disable zoom and pan" : "Enable zoom and pan";
    zoomToggle.setAttribute("aria-label", label);
    zoomToggle.setAttribute("aria-pressed", String(zoomEnabled));
    zoomToggle.title = label;
    zoomToggle.style.background = "transparent";
    zoomToggle.style.borderColor = "transparent";
    zoomToggle.style.color = zoomEnabled ? "#1f2937" : "#9ca3af";
    zoomToggle.style.opacity = zoomEnabled ? "1" : "0.55";
    for (const shape of [iconCircle, iconPath]) {
      shape.setAttribute("stroke-width", zoomEnabled ? "2.6" : "1.8");
    }
    svg.style.touchAction = zoomEnabled ? "none" : "auto";
    interaction.attr("cursor", transform.k > 1.000001 ? "grab" : zoomEnabled ? "crosshair" : "default");
  };
  updateZoomToggle();
  const toggleZoom = () => {
    zoomEnabled = !zoomEnabled;
    updateZoomToggle();
  };
  const setSliderZoom = () => {
    const scale = Math.exp(Math.log(maxZoomScale) * (Number(zoomSlider.value) / 1000));
    const centerX = plotLeft + plotSide / 2;
    const centerY = plotTop + plotSide / 2;
    interaction.call(
      zoomBehavior.transform,
      zoomIdentity.translate(centerX * (1 - scale), centerY * (1 - scale)).scale(scale),
    );
  };
  const resetZoom = () => {
    zoomEnabled = false;
    updateZoomToggle();
    interaction.call(zoomBehavior.transform, zoomIdentity);
  };
  zoomToggle.addEventListener("click", toggleZoom);
  zoomSlider.addEventListener("input", setSliderZoom);
  resetButton.addEventListener("click", resetZoom);
  interaction.call(zoomBehavior.transform, transform);

  const dispose = () => {
    if (disposed) return;
    disposed = true;
    interaction.on(".zoom", null).on(".ani-selection", null).on(".ani-tooltip", null);
    zoomToggle.removeEventListener("click", toggleZoom);
    zoomSlider.removeEventListener("input", setSliderZoom);
    resetButton.removeEventListener("click", resetZoom);
    signal.removeEventListener("abort", dispose);
    if (root.parentElement === container) root.remove();
    if (toolbar.parentElement === toolbarHost) toolbar.remove();
    if (positionToolbarHost && toolbarHost.style.position === "relative") {
      toolbarHost.style.position = previousToolbarHostPosition;
    }
    container.style.display = previousContainerDisplay;
    container.style.alignItems = previousContainerAlignItems;
    container.style.justifyContent = previousContainerJustifyContent;
  };
  signal.addEventListener("abort", dispose, { once: true });
  return {
    dispose,
    resize: renderZoom,
    getZoomFocus: () => ({
      scale: transform.k,
      centerX: (transform.invertX(plotLeft + plotSide / 2) - plotLeft) / plotSide,
      centerY: (transform.invertY(plotTop + plotSide / 2) - plotTop) / plotSide,
    }),
    getClustering: () => clustering,
  };
};

function internalNodes(root: ClusterNode | null): ClusterNode[] {
  if (!root || root.leafIndex !== undefined) return [];
  const nodes: ClusterNode[] = [];
  const visit = (node: ClusterNode | null) => {
    if (!node || node.leafIndex !== undefined) return;
    visit(node.left);
    visit(node.right);
    nodes.push(node);
  };
  visit(root);
  return nodes;
}

function clusterDepth(node: ClusterNode, leafCount: number, maxDistance: number): number {
  if (node.leafIndex !== undefined) return 0;
  if (maxDistance > 0) return Math.max(0, Math.min(1, node.distance / maxDistance));
  return Math.max(0, Math.min(1, (node.id - leafCount + 1) / Math.max(1, leafCount - 1)));
}

function rowDendrogramPath(
  nodes: ClusterNode[],
  leafCount: number,
  maxDistance: number,
  transform: ZoomTransform,
  plotTop: number,
  cellSize: number,
  treeWidth: number,
): string {
  const nodeX = (node: ClusterNode) =>
    treeWidth * (1 - clusterDepth(node, leafCount, maxDistance));
  const nodeY = (node: ClusterNode) =>
    transform.applyY(plotTop + ((node.leafStart + node.leafEnd) / 2) * cellSize);
  const parts: string[] = [];
  for (const node of nodes) {
    if (!node.left || !node.right) continue;
    const parentX = nodeX(node);
    const leftY = nodeY(node.left);
    const rightY = nodeY(node.right);
    parts.push(
      `M${parentX},${leftY}V${rightY}`,
      `M${parentX},${leftY}H${nodeX(node.left)}`,
      `M${parentX},${rightY}H${nodeX(node.right)}`,
    );
  }
  return parts.join("");
}

function columnDendrogramPath(
  nodes: ClusterNode[],
  leafCount: number,
  maxDistance: number,
  transform: ZoomTransform,
  plotLeft: number,
  plotBottom: number,
  labelSpace: number,
  treeHeight: number,
  cellSize: number,
): string {
  const treeTop = plotBottom + labelSpace;
  const nodeY = (node: ClusterNode) =>
    treeTop + clusterDepth(node, leafCount, maxDistance) * treeHeight;
  const nodeX = (node: ClusterNode) =>
    transform.applyX(plotLeft + ((node.leafStart + node.leafEnd) / 2) * cellSize);
  const parts: string[] = [];
  for (const node of nodes) {
    if (!node.left || !node.right) continue;
    const parentY = nodeY(node);
    const leftX = nodeX(node.left);
    const rightX = nodeX(node.right);
    parts.push(
      `M${leftX},${nodeY(node.left)}V${parentY}`,
      `M${rightX},${nodeY(node.right)}V${parentY}`,
      `M${leftX},${parentY}H${rightX}`,
    );
  }
  return parts.join("");
}

function truncateLabel(label: string, maxWidth: number, fontSize: number): string {
  const maxCharacters = Math.floor(maxWidth / Math.max(fontSize * 0.54, 1));
  if (maxCharacters <= 0) return "";
  if (label.length <= maxCharacters) return label;
  return maxCharacters <= 1 ? "…" : `${label.slice(0, maxCharacters - 1)}…`;
}

function formatValue(value: number, isDistance = false): string {
  return isDistance ? String(Math.round(value)) : `${value.toFixed(2)}%`;
}

function formatCellValue(value: number, isDistance = false): string {
  return isDistance ? String(Math.round(value)) : Number(value.toPrecision(3)).toString();
}

function cssDimension(element: HTMLElement, property: string, fallback: number): number {
  const value = Number.parseFloat(getComputedStyle(element).getPropertyValue(property));
  return Number.isFinite(value) && value > 0 ? value : fallback;
}

function showEmptyState(container: HTMLElement, signal: AbortSignal, metricLabel?: string): () => void {
  const message = document.createElement("p");
  message.textContent = `${metricLabel ?? "Similarity"} data not available.`;
  container.replaceChildren(message);
  const cleanup = () => {
    signal.removeEventListener("abort", cleanup);
    if (message.parentElement === container) message.remove();
  };
  signal.addEventListener("abort", cleanup, { once: true });
  return cleanup;
}

declare global {
  interface Window {
    HTMLWidgets?: {
      widget: (definition: {
          name: string;
          type: "output";
          factory: (element: HTMLElement) => {
            renderValue: (payload: unknown) => void;
          resize: (width?: number, height?: number) => void;
        };
      }) => void;
    };
  }
}

if (window.HTMLWidgets) {
  window.HTMLWidgets.widget({
    name: "psHeatmap",
    type: "output",
    factory: (element) => {
      let instance: HeatmapInstance | undefined;
      let controller: AbortController | undefined;
      let payload: unknown;
      let clustering: Clustering | undefined;
      const render = (width?: number, height?: number, preserveZoom = false) => {
        const zoomFocus = preserveZoom ? instance?.getZoomFocus() : undefined;
        if (!preserveZoom) clustering = undefined;
        if (Number.isFinite(width) && Number(width) > 0) {
          element.style.setProperty("--widget-max-width", `${Number(width)}px`);
        }
        if (Number.isFinite(height) && Number(height) > 0) {
          element.style.setProperty("--widget-max-height", `${Number(height)}px`);
        }
        controller?.abort();
        controller = new AbortController();
        instance = mountHeatmap(element, payload, {
          signal: controller.signal,
          zoomFocus,
          clustering,
        });
        clustering = instance.getClustering();
      };
      return {
        renderValue: (value) => {
          payload = value;
          render(element.clientWidth, element.clientHeight, false);
        },
        resize: (width, height) => render(width, height, true),
      };
    },
  });
}
