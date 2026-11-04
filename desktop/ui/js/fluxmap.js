// The flux map: a canvas that draws the network and the flow through it. Line width and colour show
// the size of a flux, and particles travel in its direction at a speed that grows with it. A static
// layer holds lines and nodes, and only the particles draw every frame.

import { ramp, REDUCED_MOTION } from './core.js';

const SPRITE_STEPS = 12;
const BACKGROUND = '#14171a';
const IDLE_LINE = '#343a41';
const CUT = '#ef5a6f';

function makeSprites() {
  const sprites = [];
  for (let k = 0; k < SPRITE_STEPS; k++) {
    const c = document.createElement('canvas');
    c.width = c.height = 28;
    const g = c.getContext('2d');
    const colour = ramp(k / (SPRITE_STEPS - 1));
    const gradient = g.createRadialGradient(14, 14, 0, 14, 14, 14);
    gradient.addColorStop(0, '#fffdf2');
    gradient.addColorStop(0.16, colour);
    gradient.addColorStop(0.42, colour.replace('rgb', 'rgba').replace(')', ',0.32)'));
    gradient.addColorStop(1, colour.replace('rgb', 'rgba').replace(')', ',0)'));
    g.fillStyle = gradient;
    g.fillRect(0, 0, 28, 28);
    sprites.push(c);
  }
  return sprites;
}

function quad(e, t) {
  const u = 1 - t;
  return [u * u * e.x0 + 2 * u * t * e.cx + t * t * e.x1, u * u * e.y0 + 2 * u * t * e.cy + t * t * e.y1];
}

function curveLength(e) {
  let total = 0;
  let [px, py] = quad(e, 0);
  for (let k = 1; k <= 8; k++) {
    const [x, y] = quad(e, k / 8);
    total += Math.hypot(x - px, y - py);
    px = x;
    py = y;
  }
  return total;
}

export class FluxMap {
  constructor(canvas, handlers = {}) {
    this.canvas = canvas;
    this.ctx = canvas.getContext('2d', { alpha: false });
    this.handlers = handlers;
    this.layer = document.createElement('canvas');
    this.sprites = makeSprites();
    this.view = { x: 0, y: 0, k: 1 };
    this.size = { w: 1, h: 1, dpr: 1 };
    this.reactions = [];
    this.nodes = [];
    this.nodeById = new Map();
    this.particles = [];
    this.shown = [];
    this.target = [];
    this.disabled = new Set();
    this.hover = null;
    this.selected = null;
    this.layerDirty = true;
    this.transition = 1;
    this.reveal = REDUCED_MOTION ? 1 : 0;
    this.revealStart = null;
    this.flashes = [];
    this.frames = 0;
    this.frameClock = performance.now();
    this.fps = 0;
    this.running = false;
    this.biomass = null;
    this.growth = null;
    this.bindEvents();
    new ResizeObserver(() => this.resize()).observe(canvas);
  }

  /** Builds geometry from a layout, in the reaction order of the model. */
  setLayout(layout, reactionIds, keepView = false) {
    this.layout = layout;
    this.nodes = Object.entries(layout.nodes).map(([id, n]) => ({ id, ...n, flux: 0 }));
    this.nodeById = new Map(this.nodes.map((n) => [n.id, n]));
    const [W, H] = layout.world;
    this.world = { w: W, h: H };
    this.reactions = [];
    this.indexOf = new Map();
    const position = (id) => this.nodeById.get(id);
    reactionIds.forEach((id, modelIndex) => {
      const entry = layout.reactions[id];
      if (!entry) return;
      const subs = entry.s.map(position).filter(Boolean);
      const prods = entry.p.map(position).filter(Boolean);
      if (entry.biomass) {
        this.biomass = { modelIndex, x: entry.m[0], y: entry.m[1], id };
        return;
      }
      const all = [...subs, ...prods];
      if (all.length === 0) return;
      let mx = all.reduce((a, n) => a + n.x, 0) / all.length;
      let my = all.reduce((a, n) => a + n.y, 0) / all.length;
      if (entry.m) [mx, my] = entry.m;
      const incoming = [];
      const outgoing = [];
      const stub = (n) => {
        // Exchange reactions leave the map towards its nearest edge.
        const dLeft = n.x;
        const dRight = W - n.x;
        const dTop = n.y;
        const dBottom = H - n.y;
        const least = Math.min(dLeft, dRight, dTop, dBottom);
        const reach = Math.max(40, least);
        if (least === dLeft) return { x: n.x - reach, y: n.y };
        if (least === dRight) return { x: n.x + reach, y: n.y };
        if (least === dTop) return { x: n.x, y: n.y - reach };
        return { x: n.x, y: n.y + reach };
      };
      if (entry.boundary && subs.length + prods.length === 1) {
        const n = all[0];
        const end = stub(n);
        const edge = { x0: n.x, y0: n.y, cx: (n.x + end.x) / 2, cy: (n.y + end.y) / 2, x1: end.x, y1: end.y };
        edge.length = curveLength(edge);
        // A positive flux leaves the cell when the metabolite is a substrate of the exchange.
        (subs.length ? outgoing : incoming).push(subs.length ? edge : { x0: end.x, y0: end.y, cx: edge.cx, cy: edge.cy, x1: n.x, y1: n.y, length: edge.length });
        mx = (n.x + end.x) / 2;
        my = (n.y + end.y) / 2;
      } else {
        for (const n of subs) {
          const bend = subs.length === 1 && prods.length === 1 ? 0 : 0.42;
          const cx = mx + (n.x - mx) * bend;
          const cy = my + (n.y - my) * bend;
          const edge = { x0: n.x, y0: n.y, cx: subs.length === 1 && prods.length === 1 ? (n.x + mx) / 2 : cx, cy: subs.length === 1 && prods.length === 1 ? (n.y + my) / 2 : cy, x1: mx, y1: my };
          edge.length = curveLength(edge);
          incoming.push(edge);
        }
        for (const n of prods) {
          const bend = subs.length === 1 && prods.length === 1 ? 0 : 0.42;
          const cx = mx + (n.x - mx) * bend;
          const cy = my + (n.y - my) * bend;
          const edge = { x0: mx, y0: my, cx: subs.length === 1 && prods.length === 1 ? (n.x + mx) / 2 : cx, cy: subs.length === 1 && prods.length === 1 ? (n.y + my) / 2 : cy, x1: n.x, y1: n.y };
          edge.length = curveLength(edge);
          outgoing.push(edge);
        }
      }
      const samples = [];
      for (const e of [...incoming, ...outgoing]) {
        for (let k = 0; k <= 8; k++) samples.push(quad(e, k / 8));
      }
      this.indexOf.set(modelIndex, this.reactions.length);
      this.reactions.push({ id, modelIndex, incoming, outgoing, mx, my, samples, subs, prods, boundary: !!entry.boundary });
    });
    this.shown = new Array(reactionIds.length).fill(0);
    this.target = new Array(reactionIds.length).fill(0);
    this.maxFlux = 1;
    this.ranked = null;
    this.particles = [];
    this.layerDirty = true;
    this.computeReveal();
    if (!keepView) this.fit();
  }

  /** Distance from the carbon source, in reactions, for the opening reveal. */
  computeReveal() {
    const start = this.nodeById.get('glc__D_e') ?? this.nodes[0];
    const dist = new Map();
    if (!start) return;
    dist.set(start.id, 0);
    let frontier = [start.id];
    let depth = 0;
    const neighbours = new Map();
    for (const r of this.reactions) {
      const members = [...r.subs, ...r.prods].map((n) => n.id);
      for (const a of members) {
        if (!neighbours.has(a)) neighbours.set(a, new Set());
        for (const b of members) if (a !== b) neighbours.get(a).add(b);
      }
    }
    while (frontier.length) {
      depth += 1;
      const next = [];
      for (const a of frontier) {
        for (const b of neighbours.get(a) ?? []) {
          if (!dist.has(b)) {
            dist.set(b, depth);
            next.push(b);
          }
        }
      }
      frontier = next;
    }
    const most = Math.max(1, depth);
    for (const n of this.nodes) n.order = (dist.get(n.id) ?? most) / most;
    for (const r of this.reactions) {
      const members = [...r.subs, ...r.prods];
      r.order = members.length ? Math.min(...members.map((n) => n.order ?? 1)) : 1;
    }
  }

  /** Sets new fluxes; lines and particles move to them smoothly. */
  setFluxes(fluxes, disabled, growth) {
    this.target = fluxes.slice();
    this.disabled = new Set(disabled);
    this.growth = growth;
    let top = 1e-9;
    const carrying = [];
    for (const r of this.reactions) {
      const size = Math.abs(fluxes[r.modelIndex] ?? 0);
      top = Math.max(top, size);
      if (size > 1e-7) carrying.push(size);
    }
    this.maxFlux = top;
    // Fluxes on a large map cluster, so its colour shows rank among the fluxes on view.
    this.ranked = this.large ? carrying.sort((a, b) => a - b) : null;
    this.transition = 0;
    this.seedParticles();
    this.layerDirty = true;
  }

  get large() {
    return this.reactions.length > 150;
  }

  magnitude(v) {
    if (this.maxFlux <= 0) return 0;
    if (this.ranked?.length > 1) {
      const size = Math.abs(v);
      let low = 0;
      let high = this.ranked.length;
      while (low < high) {
        const middle = (low + high) >> 1;
        if (this.ranked[middle] < size) low = middle + 1;
        else high = middle;
      }
      return low / (this.ranked.length - 1);
    }
    return Math.min(1, Math.log1p(Math.abs(v) * 8) / Math.log1p(this.maxFlux * 8));
  }

  seedParticles() {
    const old = this.particles;
    this.particles = [];
    let kept = 0;
    for (const r of this.reactions) {
      const v = this.target[r.modelIndex] ?? 0;
      if (Math.abs(v) < 1e-7 || this.disabled.has(r.id)) continue;
      const m = this.magnitude(v);
      const count = Math.max(1, Math.round(1 + m * 7));
      const pathLength = (r.incoming[0]?.length ?? 0) + (r.outgoing[0]?.length ?? 0);
      const dense = Math.max(1, Math.round(pathLength / 90));
      // Large maps thin out, so that flow stays readable and the frame rate stays up.
      const thin = this.large ? 0.3 : 1;
      for (let k = 0; k < Math.max(1, Math.round(count * Math.min(dense, 4) * thin)); k++) {
        const previous = old[kept];
        kept += 1;
        this.particles.push({
          r,
          t: previous?.r === r ? previous.t : Math.random(),
          inEdge: r.incoming.length ? r.incoming[k % r.incoming.length] : null,
          outEdge: r.outgoing.length ? r.outgoing[k % r.outgoing.length] : null,
        });
      }
    }
  }

  fit() {
    // The space left free by the panels that float over the map.
    const inset = this.handlers.insets?.() ?? { left: 30, top: 30, right: 30, bottom: 30 };
    const w = Math.max(100, this.size.w - inset.left - inset.right);
    const hh = Math.max(100, this.size.h - inset.top - inset.bottom);
    const k = Math.min(w / this.world.w, hh / this.world.h);
    this.view.k = k;
    this.view.x = inset.left + (w - this.world.w * k) / 2;
    this.view.y = inset.top + (hh - this.world.h * k) / 2;
    this.layerDirty = true;
  }

  /** Moves the view so a world point is centred. */
  focus(x, y, k = null) {
    if (k) this.view.k = k;
    this.view.x = this.size.w / 2 - x * this.view.k;
    this.view.y = this.size.h / 2 - y * this.view.k;
    this.layerDirty = true;
  }

  resize() {
    const box = this.canvas.getBoundingClientRect();
    const dpr = Math.min(window.devicePixelRatio || 1, 2);
    const first = this.size.w <= 1;
    this.size = { w: Math.max(1, box.width), h: Math.max(1, box.height), dpr };
    this.canvas.width = Math.round(this.size.w * dpr);
    this.canvas.height = Math.round(this.size.h * dpr);
    this.layer.width = this.canvas.width;
    this.layer.height = this.canvas.height;
    if (this.world && first) this.fit();
    this.layerDirty = true;
  }

  start() {
    if (this.running) return;
    this.running = true;
    this.last = performance.now();
    const loop = (now) => {
      if (!this.running) return;
      // Scheduled first, so one failed frame cannot stop the map for good.
      requestAnimationFrame(loop);
      this.frame(now);
    };
    requestAnimationFrame(loop);
  }

  stop() {
    this.running = false;
  }

  worldToScreen(x, y) {
    return [x * this.view.k + this.view.x, y * this.view.k + this.view.y];
  }

  screenToWorld(x, y) {
    return [(x - this.view.x) / this.view.k, (y - this.view.y) / this.view.k];
  }

  frame(now) {
    const dt = Math.min(0.05, (now - this.last) / 1000);
    this.last = now;
    if (!this.layout) return;
    if (this.revealStart === null) this.revealStart = now;
    if (this.reveal < 1) {
      this.reveal = Math.min(1, (now - this.revealStart) / 1300);
      this.layerDirty = true;
    }
    if (this.transition < 1) {
      this.transition = Math.min(1, this.transition + dt / 0.45);
      const ease = 1 - (1 - this.transition) ** 3;
      for (const r of this.reactions) {
        const i = r.modelIndex;
        this.shown[i] += (this.target[i] - this.shown[i]) * Math.min(1, ease * 0.35 + 0.08);
        if (this.transition === 1) this.shown[i] = this.target[i];
      }
      this.layerDirty = true;
    }
    if (this.layerDirty) this.drawLayer();
    const g = this.ctx;
    g.setTransform(1, 0, 0, 1, 0, 0);
    g.drawImage(this.layer, 0, 0);
    if (!REDUCED_MOTION) this.drawParticles(g, dt);
    this.drawFlashes(g, now);
    this.frames += 1;
    if (now - this.frameClock > 1000) {
      this.fps = (this.frames * 1000) / (now - this.frameClock);
      this.frames = 0;
      this.frameClock = now;
      this.handlers.onFps?.(this.fps);
    }
  }

  reachedByReveal(order) {
    return this.reveal >= 1 ? 1 : Math.max(0, Math.min(1, (this.reveal * 1.5 - (order ?? 0)) * 3));
  }

  drawPath(g, r, width, colour, alpha = 1, dash = null) {
    g.lineWidth = width;
    g.strokeStyle = colour;
    g.globalAlpha = alpha;
    g.setLineDash(dash ?? []);
    g.beginPath();
    for (const e of [...r.incoming, ...r.outgoing]) {
      g.moveTo(e.x0, e.y0);
      g.quadraticCurveTo(e.cx, e.cy, e.x1, e.y1);
    }
    g.stroke();
  }

  drawLayer() {
    const g = this.layer.getContext('2d', { alpha: false });
    const { dpr } = this.size;
    g.setTransform(1, 0, 0, 1, 0, 0);
    g.fillStyle = BACKGROUND;
    g.fillRect(0, 0, this.layer.width, this.layer.height);
    g.setTransform(dpr * this.view.k, 0, 0, dpr * this.view.k, dpr * this.view.x, dpr * this.view.y);
    const k = this.view.k;
    g.lineCap = 'round';

    // Region names sit under everything, large and quiet.
    g.textBaseline = 'middle';
    g.font = `500 ${Math.max(15, 15 / Math.sqrt(k))}px "Public Sans", sans-serif`;
    g.fillStyle = '#4a515a';
    for (const region of this.layout.regions ?? []) {
      g.textAlign = region.anchor === 'end' ? 'right' : region.anchor === 'start' ? 'left' : 'center';
      g.fillText(region.label, region.x, region.y);
    }

    // Throughput per node, for the rings.
    for (const n of this.nodes) n.flux = 0;
    for (const r of this.reactions) {
      const v = Math.abs(this.shown[r.modelIndex] ?? 0);
      for (const n of [...r.subs, ...r.prods]) n.flux = Math.max(n.flux, v);
    }

    const widthOf = (m) => (1.2 + m * 7.5) / Math.max(0.65, Math.sqrt(k));
    // Idle reactions first, then the ones carrying flux, largest last, so that big flows read on top.
    const ordered = [...this.reactions].sort((a, b) => Math.abs(this.shown[a.modelIndex] ?? 0) - Math.abs(this.shown[b.modelIndex] ?? 0));
    for (const r of ordered) {
      const shown = this.reachedByReveal(r.order);
      if (shown <= 0) continue;
      const v = this.shown[r.modelIndex] ?? 0;
      if (this.disabled.has(r.id)) {
        this.drawPath(g, r, 2 / Math.max(0.65, Math.sqrt(k)), CUT, 0.72 * shown, [8, 7]);
      } else if (Math.abs(v) < 1e-7) {
        this.drawPath(g, r, 1.3 / Math.max(0.65, Math.sqrt(k)), IDLE_LINE, shown);
      } else {
        const m = this.magnitude(v);
        this.drawPath(g, r, widthOf(m) * 2.4, ramp(m), 0.1 * shown);
        this.drawPath(g, r, widthOf(m), ramp(m), 0.92 * shown);
      }
    }
    g.setLineDash([]);
    g.globalAlpha = 1;

    // Arrowhead at each product on carrying reactions, showing the net direction.
    for (const r of this.reactions) {
      const v = this.shown[r.modelIndex] ?? 0;
      if (Math.abs(v) < 1e-7 || this.disabled.has(r.id) || this.reachedByReveal(r.order) < 1) continue;
      const edges = v > 0 ? r.outgoing : r.incoming;
      const m = this.magnitude(v);
      g.fillStyle = ramp(m);
      for (const e of edges) {
        const from = v > 0 ? quad(e, 0.8) : quad(e, 0.2);
        const to = v > 0 ? quad(e, 0.93) : quad(e, 0.07);
        const angle = Math.atan2(to[1] - from[1], to[0] - from[0]);
        const size = (5 + m * 7) / Math.max(0.65, Math.sqrt(k));
        g.beginPath();
        g.moveTo(to[0] + Math.cos(angle) * size * 0.6, to[1] + Math.sin(angle) * size * 0.6);
        g.lineTo(to[0] + Math.cos(angle + 2.5) * size, to[1] + Math.sin(angle + 2.5) * size);
        g.lineTo(to[0] + Math.cos(angle - 2.5) * size, to[1] + Math.sin(angle - 2.5) * size);
        g.closePath();
        g.fill();
      }
    }

    // Selection and hover.
    for (const [target, width, alpha] of [[this.hover, 3, 0.75], [this.selected, 3.5, 1]]) {
      if (target?.kind === 'reaction') {
        const r = this.reactions[this.indexOf.get(this.handlers.indexOfId?.(target.id))];
        if (r) this.drawPath(g, r, (width + 2) / Math.max(0.65, Math.sqrt(k)), '#ece8df', alpha);
      }
    }

    // Cut marks on deleted reactions.
    g.setLineDash([]);
    g.strokeStyle = CUT;
    g.lineWidth = 3 / Math.max(0.65, Math.sqrt(k));
    for (const r of this.reactions) {
      if (this.disabled.has(r.id)) {
        const s = 8 / Math.max(0.65, Math.sqrt(k));
        g.beginPath();
        g.moveTo(r.mx - s, r.my - s);
        g.lineTo(r.mx + s, r.my + s);
        g.moveTo(r.mx + s, r.my - s);
        g.lineTo(r.mx - s, r.my + s);
        g.stroke();
      }
    }

    // Nodes and labels. Labels keep a constant size on screen until zoomed well in.
    const fontSize = Math.min(15, 12.5 / Math.max(0.55, k));
    const namesFit = fontSize * k >= 8;
    g.font = `550 ${fontSize}px "Public Sans", sans-serif`;
    g.textAlign = 'center';
    g.textBaseline = 'alphabetic';
    for (const n of this.nodes) {
      const shown = this.reachedByReveal(n.order);
      if (shown <= 0) continue;
      g.globalAlpha = shown;
      const m = n.flux > 1e-7 ? this.magnitude(n.flux) : 0;
      const radius = (n.external ? 5.5 : 7.5) / Math.max(0.65, Math.sqrt(k));
      g.beginPath();
      g.arc(n.x, n.y, radius, 0, Math.PI * 2);
      g.fillStyle = n.external ? BACKGROUND : '#20252a';
      g.fill();
      g.lineWidth = 2 / Math.max(0.65, Math.sqrt(k));
      g.strokeStyle = n.flux > 1e-7 ? ramp(m) : '#566070';
      g.stroke();
      const selected = this.selected?.kind === 'metabolite' && this.selected.id === n.id;
      const hovered = this.hover?.kind === 'metabolite' && this.hover.id === n.id;
      if (selected || hovered) {
        g.beginPath();
        g.arc(n.x, n.y, radius * 1.9, 0, Math.PI * 2);
        g.strokeStyle = '#ece8df';
        g.stroke();
      }
      // Names too small to read are left out until the map is zoomed in.
      if (!namesFit && !selected && !hovered) continue;
      // A dark halo keeps labels legible where lines pass.
      const gap = radius + 5 / Math.max(0.65, Math.sqrt(k));
      let lx = n.x;
      let ly = n.y - gap;
      g.textAlign = 'center';
      g.textBaseline = 'alphabetic';
      if (n.side === 'l') {
        lx = n.x - gap;
        ly = n.y;
        g.textAlign = 'right';
        g.textBaseline = 'middle';
      } else if (n.side === 'r') {
        lx = n.x + gap;
        ly = n.y;
        g.textAlign = 'left';
        g.textBaseline = 'middle';
      } else if (n.side === 'b') {
        ly = n.y + gap;
        g.textBaseline = 'top';
      }
      g.lineWidth = fontSize * 0.42;
      g.strokeStyle = BACKGROUND;
      g.lineJoin = 'round';
      g.strokeText(n.label, lx, ly);
      g.fillStyle = n.flux > 1e-7 ? '#ece8df' : '#8b929b';
      g.fillText(n.label, lx, ly);
    }
    g.textAlign = 'center';
    g.textBaseline = 'alphabetic';
    g.globalAlpha = 1;

    // Reaction names and values appear once there is room, and always for selection and hover.
    const labelSize = Math.min(13, 11 / Math.max(0.55, k));
    g.font = `500 ${labelSize}px "Public Sans", sans-serif`;
    for (const r of this.reactions) {
      const isOn = this.selected?.id === r.id || this.hover?.id === r.id;
      if (!(k > 1.05 || isOn)) continue;
      const v = this.shown[r.modelIndex] ?? 0;
      const text = `${r.id}  ${Math.abs(v) < 1e-7 ? '0' : v.toFixed(Math.abs(v) >= 100 ? 0 : 2)}`;
      g.lineWidth = labelSize * 0.42;
      g.strokeStyle = BACKGROUND;
      g.strokeText(text, r.mx, r.my + 16 / Math.max(0.65, k));
      g.fillStyle = isOn ? '#ece8df' : '#a3a9b0';
      g.fillText(text, r.mx, r.my + 16 / Math.max(0.65, k));
    }

    // Biomass, the sink the whole map serves: its ring grows with the growth rate.
    if (this.biomass) {
      const b = this.biomass;
      const value = this.growth ?? 0;
      const reference = this.handlers.referenceGrowth?.() || value || 1;
      const share = Math.max(0, Math.min(1.2, value / reference));
      const shown = this.reachedByReveal(1);
      g.globalAlpha = shown;
      g.beginPath();
      g.arc(b.x, b.y, 44, 0, Math.PI * 2);
      g.fillStyle = '#1b1f23';
      g.fill();
      g.lineWidth = 2;
      g.strokeStyle = '#343a41';
      g.stroke();
      g.beginPath();
      g.arc(b.x, b.y, 44, -Math.PI / 2, -Math.PI / 2 + Math.PI * 2 * Math.min(1, share));
      g.lineWidth = 6;
      g.strokeStyle = value > 1e-9 ? ramp(0.95) : CUT;
      g.stroke();
      g.fillStyle = '#ece8df';
      g.font = '600 20px "Public Sans", sans-serif';
      g.textBaseline = 'middle';
      g.fillText(value > 1e-9 ? value.toFixed(3) : '0', b.x, b.y - 3);
      g.font = '500 12px "Public Sans", sans-serif';
      g.fillStyle = '#a3a9b0';
      g.fillText('biomass', b.x, b.y + 18);
      g.textBaseline = 'alphabetic';
      g.globalAlpha = 1;
    }
    this.layerDirty = false;
  }

  drawParticles(g, dt) {
    const { dpr } = this.size;
    const k = this.view.k;
    g.setTransform(dpr, 0, 0, dpr, 0, 0);
    g.globalCompositeOperation = 'lighter';
    // Particles add up where they overlap, so a dense map uses smaller, fainter ones.
    const [baseSize, growSize, alpha] = this.large ? [4, 7, 0.5] : [7, 13, 0.92];
    for (const p of this.particles) {
      const v = this.target[p.r.modelIndex] ?? 0;
      if (this.reachedByReveal(p.r.order) < 1) continue;
      const m = this.magnitude(v);
      const speed = (70 + 420 * m);
      const total = (p.inEdge?.length ?? 0) + (p.outEdge?.length ?? 0);
      if (total <= 0) continue;
      p.t += (dt * speed) / total;
      if (p.t > 1) p.t -= 1;
      // Position along the chain substrate -> midpoint -> product, reversed for a negative flux.
      const along = (v >= 0 ? p.t : 1 - p.t) * total;
      let point;
      const first = p.inEdge?.length ?? 0;
      if (along < first) point = quad(p.inEdge, along / first);
      else point = quad(p.outEdge ?? p.inEdge, (along - first) / Math.max(1e-9, total - first));
      const [sx, sy] = this.worldToScreen(point[0], point[1]);
      if (sx < -20 || sy < -20 || sx > this.size.w + 20 || sy > this.size.h + 20) continue;
      const sprite = this.sprites[Math.round(m * (SPRITE_STEPS - 1))];
      const size = (baseSize + m * growSize) * Math.min(1.5, Math.max(0.75, Math.sqrt(k)));
      g.globalAlpha = alpha;
      g.drawImage(sprite, sx - size / 2, sy - size / 2, size, size);
    }
    g.globalCompositeOperation = 'source-over';
    g.globalAlpha = 1;
  }

  flash(id) {
    const r = this.reactions[this.indexOf.get(this.handlers.indexOfId?.(id))];
    if (r && !REDUCED_MOTION) this.flashes.push({ x: r.mx, y: r.my, start: performance.now() });
  }

  drawFlashes(g, now) {
    if (!this.flashes.length) return;
    const { dpr } = this.size;
    g.setTransform(dpr, 0, 0, dpr, 0, 0);
    this.flashes = this.flashes.filter((f) => now - f.start < 600);
    for (const f of this.flashes) {
      const t = (now - f.start) / 600;
      const [sx, sy] = this.worldToScreen(f.x, f.y);
      g.beginPath();
      g.arc(sx, sy, 6 + t * 46, 0, Math.PI * 2);
      g.lineWidth = 3 * (1 - t);
      g.strokeStyle = CUT;
      g.globalAlpha = 1 - t;
      g.stroke();
    }
    g.globalAlpha = 1;
  }

  pick(sx, sy) {
    const [x, y] = this.screenToWorld(sx, sy);
    const k = this.view.k;
    for (const n of this.nodes) {
      if (Math.hypot(n.x - x, n.y - y) < 13 / Math.max(0.5, k)) return { kind: 'metabolite', id: n.id };
    }
    if (this.biomass && Math.hypot(this.biomass.x - x, this.biomass.y - y) < 46) return { kind: 'reaction', id: this.biomass.id };
    let best = null;
    let bestD = 11 / Math.max(0.5, k);
    for (const r of this.reactions) {
      for (const [px, py] of r.samples) {
        const d = Math.hypot(px - x, py - y);
        if (d < bestD) {
          bestD = d;
          best = r;
        }
      }
    }
    return best ? { kind: 'reaction', id: best.id } : null;
  }

  bindEvents() {
    const c = this.canvas;
    let drag = null;
    c.addEventListener('pointerdown', (e) => {
      c.setPointerCapture(e.pointerId);
      drag = { x: e.clientX, y: e.clientY, vx: this.view.x, vy: this.view.y, moved: false };
    });
    c.addEventListener('pointermove', (e) => {
      const box = c.getBoundingClientRect();
      if (drag) {
        const dx = e.clientX - drag.x;
        const dy = e.clientY - drag.y;
        if (Math.abs(dx) + Math.abs(dy) > 3) drag.moved = true;
        if (drag.moved) {
          this.view.x = drag.vx + dx;
          this.view.y = drag.vy + dy;
          this.layerDirty = true;
          c.style.cursor = 'grabbing';
        }
        return;
      }
      const hit = this.pick(e.clientX - box.left, e.clientY - box.top);
      if (hit?.id !== this.hover?.id) {
        this.hover = hit;
        this.layerDirty = true;
        c.style.cursor = hit ? 'pointer' : 'grab';
        this.handlers.onHover?.(hit, e.clientX, e.clientY);
      } else if (hit) {
        this.handlers.onHover?.(hit, e.clientX, e.clientY);
      }
    });
    c.addEventListener('pointerleave', () => {
      this.hover = null;
      this.layerDirty = true;
      this.handlers.onHover?.(null);
    });
    c.addEventListener('pointerup', (e) => {
      const moved = drag?.moved;
      drag = null;
      c.style.cursor = this.hover ? 'pointer' : 'grab';
      if (moved) return;
      const box = c.getBoundingClientRect();
      const hit = this.pick(e.clientX - box.left, e.clientY - box.top);
      this.handlers.onSelect?.(hit, e);
    });
    c.addEventListener('dblclick', (e) => {
      const box = c.getBoundingClientRect();
      const hit = this.pick(e.clientX - box.left, e.clientY - box.top);
      if (hit?.kind === 'reaction') this.handlers.onCut?.(hit.id);
    });
    c.addEventListener('wheel', (e) => {
      e.preventDefault();
      const box = c.getBoundingClientRect();
      const mx = e.clientX - box.left;
      const my = e.clientY - box.top;
      const factor = Math.exp(-e.deltaY * 0.0014);
      const k = Math.min(6, Math.max(0.2, this.view.k * factor));
      const [wx, wy] = this.screenToWorld(mx, my);
      this.view.k = k;
      this.view.x = mx - wx * k;
      this.view.y = my - wy * k;
      this.layerDirty = true;
    }, { passive: false });
  }

  select(target) {
    this.selected = target;
    this.layerDirty = true;
  }

  zoomBy(factor) {
    const cx = this.size.w / 2;
    const cy = this.size.h / 2;
    const [wx, wy] = this.screenToWorld(cx, cy);
    this.view.k = Math.min(6, Math.max(0.2, this.view.k * factor));
    this.view.x = cx - wx * this.view.k;
    this.view.y = cy - wy * this.view.k;
    this.layerDirty = true;
  }

  /** A slight move, used when fonts arrive after the first draw. */
  redraw() {
    this.layerDirty = true;
  }
}
