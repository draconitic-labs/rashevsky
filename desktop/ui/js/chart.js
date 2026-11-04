// Small plotting library over canvas: lines with bands, scatter and heatmaps, with crisp axes.

const INK = '#a3a9b0';
const FAINT = '#2a3036';
const LINE = '#444b53';

function ticks(lo, hi, count = 6) {
  if (!(hi > lo)) return [lo];
  const span = hi - lo;
  const step = 10 ** Math.floor(Math.log10(span / count));
  const scaled = span / count / step;
  const nice = scaled < 1.5 ? 1 : scaled < 3 ? 2 : scaled < 7 ? 5 : 10;
  const inc = nice * step;
  const out = [];
  for (let v = Math.ceil(lo / inc) * inc; v <= hi + inc * 1e-6; v += inc) out.push(Math.abs(v) < inc * 1e-9 ? 0 : v);
  return out;
}

function label(v) {
  const a = Math.abs(v);
  if (a === 0) return '0';
  if (a >= 10000 || a < 0.001) return v.toExponential(1).replace('e+', 'e').replace('e-', 'e−').replace('-', '−');
  const text = a >= 100 ? v.toFixed(0) : a >= 10 ? v.toFixed(1) : a >= 1 ? v.toFixed(2) : v.toPrecision(2);
  const trimmed = text.includes('.') ? text.replace(/\.?0+$/, '') : text;
  return trimmed.replace('-', '−');
}

function prepare(canvas) {
  const box = canvas.getBoundingClientRect();
  const dpr = Math.min(window.devicePixelRatio || 1, 2);
  canvas.width = Math.max(1, Math.round(box.width * dpr));
  canvas.height = Math.max(1, Math.round(box.height * dpr));
  const g = canvas.getContext('2d');
  g.setTransform(dpr, 0, 0, dpr, 0, 0);
  g.clearRect(0, 0, box.width, box.height);
  return { g, w: box.width, h: box.height };
}

/** Draws axes and returns scale functions. */
function axes(g, w, h, xr, yr, spec) {
  const m = { l: 66, r: 18, t: 14, b: 46, ...(spec.margin ?? {}) };
  const sx = (v) => m.l + ((v - xr[0]) / (xr[1] - xr[0] || 1)) * (w - m.l - m.r);
  const sy = (v) => h - m.b - ((v - yr[0]) / (yr[1] - yr[0] || 1)) * (h - m.t - m.b);
  g.font = '12px "Public Sans", sans-serif';
  g.lineWidth = 1;
  g.textBaseline = 'middle';
  g.textAlign = 'right';
  for (const t of ticks(yr[0], yr[1], Math.max(3, Math.floor((h - m.t - m.b) / 56)))) {
    const y = Math.round(sy(t)) + 0.5;
    g.strokeStyle = FAINT;
    g.beginPath();
    g.moveTo(m.l, y);
    g.lineTo(w - m.r, y);
    g.stroke();
    g.fillStyle = INK;
    g.fillText(label(t), m.l - 9, y);
  }
  g.textAlign = 'center';
  g.textBaseline = 'top';
  for (const t of ticks(xr[0], xr[1], Math.max(3, Math.floor((w - m.l - m.r) / 92)))) {
    const x = Math.round(sx(t)) + 0.5;
    g.strokeStyle = FAINT;
    g.beginPath();
    g.moveTo(x, m.t);
    g.lineTo(x, h - m.b);
    g.stroke();
    g.fillStyle = INK;
    g.fillText(label(t), x, h - m.b + 8);
  }
  g.strokeStyle = LINE;
  g.strokeRect(m.l + 0.5, m.t + 0.5, w - m.l - m.r - 1, h - m.t - m.b - 1);
  g.fillStyle = '#c8ccd1';
  g.font = '500 13px "Public Sans", sans-serif';
  if (spec.xlabel) g.fillText(spec.xlabel, (m.l + w - m.r) / 2, h - 20);
  if (spec.ylabel) {
    g.save();
    g.translate(15, (m.t + h - m.b) / 2);
    g.rotate(-Math.PI / 2);
    g.textBaseline = 'middle';
    g.fillText(spec.ylabel, 0, 0);
    g.restore();
  }
  return { sx, sy, m };
}

function extent(values, pad = 0.04) {
  let lo = Infinity;
  let hi = -Infinity;
  for (const v of values) {
    if (Number.isFinite(v)) {
      if (v < lo) lo = v;
      if (v > hi) hi = v;
    }
  }
  if (!Number.isFinite(lo)) return [0, 1];
  if (lo === 0 && hi === 0) return [0, 1];
  if (lo === hi) {
    const d = Math.abs(lo) * 0.1;
    return [lo - d, hi + d];
  }
  const d = (hi - lo) * pad;
  return [lo - d, hi + d];
}

/** `series`: {x, y, colour, name}; `bands`: {x, low, high, colour}. */
export function lineChart(canvas, spec) {
  const { g, w, h } = prepare(canvas);
  const xs = spec.series.flatMap((s) => s.x).concat((spec.bands ?? []).flatMap((b) => b.x));
  const ys = spec.series.flatMap((s) => s.y).concat((spec.bands ?? []).flatMap((b) => [...b.low, ...b.high]));
  const xr = spec.xrange ?? extent(xs, 0);
  const yr = spec.yrange ?? extent(ys);
  const { sx, sy, m } = axes(g, w, h, xr, yr, spec);
  g.save();
  g.beginPath();
  g.rect(m.l, m.t, w - m.l - m.r, h - m.t - m.b);
  g.clip();
  for (const b of spec.bands ?? []) {
    g.beginPath();
    b.x.forEach((v, i) => (i ? g.lineTo(sx(v), sy(b.high[i])) : g.moveTo(sx(v), sy(b.high[i]))));
    for (let i = b.x.length - 1; i >= 0; i--) g.lineTo(sx(b.x[i]), sy(b.low[i]));
    g.closePath();
    g.globalAlpha = 0.22;
    g.fillStyle = b.colour;
    g.fill();
    g.globalAlpha = 1;
  }
  g.lineJoin = 'round';
  for (const s of spec.series) {
    g.beginPath();
    let started = false;
    s.x.forEach((v, i) => {
      if (!Number.isFinite(s.y[i])) return;
      if (started) g.lineTo(sx(v), sy(s.y[i]));
      else {
        g.moveTo(sx(v), sy(s.y[i]));
        started = true;
      }
    });
    g.lineWidth = s.width ?? 2;
    g.strokeStyle = s.colour;
    g.stroke();
  }
  g.restore();
  // Legend inside the plot, top right, only when it names more than one thing.
  const named = spec.series.filter((s) => s.name);
  if (named.length > 1 && named.length <= 12) {
    g.font = '12px "Public Sans", sans-serif';
    g.textAlign = 'left';
    g.textBaseline = 'middle';
    const width = Math.max(...named.map((s) => g.measureText(s.name).width)) + 26;
    let y = m.t + 14;
    g.fillStyle = 'rgba(20,23,26,0.82)';
    g.fillRect(w - m.r - width - 10, m.t + 4, width + 6, named.length * 19 + 6);
    for (const s of named) {
      g.fillStyle = s.colour;
      g.fillRect(w - m.r - width - 4, y - 1.5, 14, 3);
      g.fillStyle = '#c8ccd1';
      g.fillText(s.name, w - m.r - width + 14, y);
      y += 19;
    }
  }
  return { sx, sy, xr, yr, m };
}

/** `points`: {x, y, colour, r}; returns a function to find the point near a screen position. */
export function scatter(canvas, spec) {
  const { g, w, h } = prepare(canvas);
  const xr = spec.xrange ?? extent(spec.points.map((p) => p.x));
  const yr = spec.yrange ?? extent(spec.points.map((p) => p.y));
  const { sx, sy, m } = axes(g, w, h, xr, yr, spec);
  g.save();
  g.beginPath();
  g.rect(m.l, m.t, w - m.l - m.r, h - m.t - m.b);
  g.clip();
  for (const p of spec.points) {
    g.beginPath();
    g.arc(sx(p.x), sy(p.y), p.r ?? 3, 0, Math.PI * 2);
    g.globalAlpha = p.alpha ?? 0.72;
    g.fillStyle = p.colour;
    g.fill();
  }
  g.globalAlpha = 1;
  if (spec.mark) {
    g.beginPath();
    g.arc(sx(spec.mark.x), sy(spec.mark.y), 7, 0, Math.PI * 2);
    g.lineWidth = 2;
    g.strokeStyle = '#ece8df';
    g.stroke();
  }
  g.restore();
  return (px, py) => {
    let best = null;
    let bestD = 12;
    spec.points.forEach((p, i) => {
      const d = Math.hypot(sx(p.x) - px, sy(p.y) - py);
      if (d < bestD) {
        bestD = d;
        best = i;
      }
    });
    return best;
  };
}

/** Diverging heatmap of a matrix, with row and column names, blue below zero and yellow above. */
export function heatmap(canvas, spec) {
  const { g, w, h } = prepare(canvas);
  const { values, rows, cols } = spec;
  const left = spec.left ?? 76;
  const top = spec.top ?? 74;
  const cw = (w - left - 12) / cols.length;
  const ch = (h - top - 12) / rows.length;
  const limit = spec.limit ?? 1;
  const colour = (v) => {
    if (v === null || !Number.isFinite(v)) return '#191c20';
    const t = Math.max(-1, Math.min(1, v / limit));
    const a = Math.abs(t) ** 0.6;
    const base = t >= 0 ? [244, 214, 78] : [92, 132, 255];
    const dark = [27, 31, 35];
    return `rgb(${[0, 1, 2].map((k) => Math.round(dark[k] + a * (base[k] - dark[k]))).join(',')})`;
  };
  values.forEach((row, i) => {
    row.forEach((v, j) => {
      g.fillStyle = colour(v);
      g.fillRect(left + j * cw, top + i * ch, Math.max(1, cw - 1), Math.max(1, ch - 1));
    });
  });
  g.fillStyle = INK;
  g.font = '11.5px "Public Sans", sans-serif';
  // The colour scale, in the free corner above the row names.
  const barWidth = left - 22;
  for (let x = 0; x < barWidth; x++) {
    g.fillStyle = colour(limit * ((2 * x) / (barWidth - 1) - 1));
    g.fillRect(10 + x, 12, 1, 9);
  }
  g.fillStyle = INK;
  g.textBaseline = 'top';
  g.textAlign = 'left';
  g.fillText(`\u2212${limit}`, 10, 25);
  g.textAlign = 'right';
  g.fillText(`+${limit}`, 10 + barWidth, 25);
  g.textBaseline = 'middle';
  rows.forEach((name, i) => g.fillText(name, left - 7, top + (i + 0.5) * ch));
  cols.forEach((name, j) => {
    g.save();
    g.translate(left + (j + 0.5) * cw, top - 7);
    g.rotate(-Math.PI / 3.2);
    g.textAlign = 'left';
    g.fillText(name, 0, 0);
    g.restore();
  });
  return (px, py) => {
    const j = Math.floor((px - left) / cw);
    const i = Math.floor((py - top) / ch);
    if (i < 0 || j < 0 || i >= rows.length || j >= cols.length) return null;
    return { row: rows[i], col: cols[j], value: values[i][j] };
  };
}

/** Histogram bars. */
export function bars(canvas, counts, colour) {
  const { g, w, h } = prepare(canvas);
  const top = Math.max(1, ...counts);
  const bw = w / counts.length;
  g.fillStyle = colour;
  counts.forEach((c, i) => {
    const bh = (c / top) * (h - 2);
    g.fillRect(i * bw, h - bh, Math.max(1, bw - 1), bh);
  });
}
