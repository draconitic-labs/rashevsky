// The open model's fluxes on the metabolism map of the Biochemistry site.

import { emit, fmt, h, invoke, on, perturbation, ramp, setModel, solve, state, toast } from '../core.js';

let refs = {};
let data = null;
const view = { x: 0, y: 0, k: 0.1 };
let drawing = false;

function fit() {
  const box = refs.canvas.getBoundingClientRect();
  const world = data.scene.world;
  view.k = Math.min(box.width / world.w, box.height / world.h) * 0.98;
  view.x = (box.width - world.w * view.k) / 2;
  view.y = (box.height - world.h * view.k) / 2;
}

function draw() {
  if (!data || drawing) return;
  drawing = true;
  requestAnimationFrame(() => {
    drawing = false;
    const c = refs.canvas;
    const box = c.getBoundingClientRect();
    const dpr = Math.min(window.devicePixelRatio || 1, 2);
    c.width = Math.max(1, box.width * dpr);
    c.height = Math.max(1, box.height * dpr);
    const g = c.getContext('2d', { alpha: false });
    g.setTransform(1, 0, 0, 1, 0, 0);
    g.fillStyle = '#14171a';
    g.fillRect(0, 0, c.width, c.height);
    g.setTransform(dpr * view.k, 0, 0, dpr * view.k, dpr * view.x, dpr * view.y);
    const { scene, overlay } = data;
    const pxl = 1 / view.k;
    const events = overlay.events;
    const top = Math.max(1e-9, ...Object.values(events).map((e) => e.flux));
    const measure = (flux) => Math.min(1, Math.log1p(flux * 8) / Math.log1p(top * 8));
    const fluxOf = (list) => {
      let total = 0;
      let kind = null;
      for (const e of list) {
        const o = events[e];
        if (o && o.flux > 1e-9) {
          total += o.flux;
          if (kind === null || o.match !== 'ec') kind = o.match;
        }
      }
      return [total, kind];
    };
    // Domain names behind everything.
    g.textAlign = 'center';
    g.textBaseline = 'middle';
    g.font = `600 ${Math.max(46, 15 * pxl)}px "Public Sans", sans-serif`;
    for (const d of scene.domains) {
      g.fillStyle = d.color;
      g.globalAlpha = 0.34;
      g.fillText(d.label, d.x, d.y);
    }
    g.globalAlpha = 1;
    // Idle first, then flux.
    g.lineCap = 'round';
    g.strokeStyle = '#2b3036';
    g.lineWidth = 1.2 * pxl;
    g.beginPath();
    for (const r of scene.reactions) {
      if (fluxOf(r.e)[0] > 0) continue;
      for (const k of [...r.i, ...r.o]) {
        const n = scene.nodes[k];
        if (n) {
          g.moveTo(r.x, r.y);
          g.lineTo(n.x, n.y);
        }
      }
    }
    g.stroke();
    const carrying = scene.reactions.map((r) => ({ r, f: fluxOf(r.e) })).filter((x) => x.f[0] > 0).sort((a, b) => a.f[0] - b.f[0]);
    for (const { r, f } of carrying) {
      const m = measure(f[0]);
      g.strokeStyle = ramp(m);
      g.lineWidth = (1.6 + 4.5 * m) * pxl;
      g.setLineDash(f[1] === 'ec' ? [7 * pxl, 6 * pxl] : []);
      g.globalAlpha = f[1] === 'ec' ? 0.75 : 1;
      g.beginPath();
      for (const k of [...r.i, ...r.o]) {
        const n = scene.nodes[k];
        if (n) {
          g.moveTo(r.x, r.y);
          g.lineTo(n.x, n.y);
        }
      }
      g.stroke();
    }
    g.setLineDash([]);
    g.globalAlpha = 1;
    g.fillStyle = '#737b85';
    for (const n of scene.nodes) {
      g.beginPath();
      g.arc(n.x, n.y, 1.8 * pxl, 0, Math.PI * 2);
      g.fill();
    }
    if (view.k > 0.32) {
      g.font = `500 ${11 * pxl}px "Public Sans", sans-serif`;
      g.fillStyle = '#c8ccd1';
      g.textAlign = 'left';
      for (const n of scene.nodes) {
        if (n.l) g.fillText(n.l, n.x + 4 * pxl, n.y);
      }
    }
  });
}

async function load() {
  if (!state.model) return;
  try {
    refs.status.textContent = 'Projecting fluxes';
    data = await invoke('atlas', { perturbation: perturbation() });
    const c = data.overlay.coverage;
    refs.status.textContent = '';
    refs.tiles.replaceChildren(
      tile(state.model.name || state.model.id, 'Model', 'wide'),
      tile(c.events_with_flux.toLocaleString('en'), 'Events carrying flux'),
      tile(c.events_exact.toLocaleString('en'), 'Events matched by identifier'),
      tile(c.events_ec_only.toLocaleString('en'), 'Events matched by EC number'),
      tile(c.events_in_scene.toLocaleString('en'), 'Events on the map'),
    );
    refs.empty.hidden = c.events_exact + c.events_ec_only > 0;
    fit();
    draw();
  } catch (error) {
    refs.status.textContent = '';
    toast(String(error));
  }
}

function tile(value, label, extra = '') {
  return h('div', { class: `tile ${extra}` }, h('strong', {}, value), h('span', {}, label));
}

async function openHuman() {
  refs.status.textContent = 'Getting Human-GEM';
  try {
    const path = await invoke('download_model', { name: 'Human-GEM' });
    refs.status.textContent = 'Opening Human-GEM';
    const model = await invoke('open_file', { path });
    setModel(model);
    await solve();
    await load();
  } catch (error) {
    refs.status.textContent = '';
    toast(String(error));
  }
}

export default {
  id: 'atlas',
  label: 'Atlas',
  icon: ['M12 21a9 9 0 1 0 0-18a9 9 0 1 0 0 18', 'M3 12h18', 'M12 3c3 3 3 15 0 18', 'M12 3c-3 3-3 15 0 18'],
  mount(el) {
    refs.canvas = h('canvas', { class: 'atlas-canvas', 'aria-label': 'Fluxes on the human metabolism map' });
    refs.tiles = h('div', { class: 'tiles' });
    refs.status = h('span', { class: 'hint' });
    refs.empty = h('div', { class: 'floating atlas-empty' }, h('h2', {}, 'This model has no human identifiers'), h('p', {}, 'The map draws human Reactome events. Human-GEM reaches them.'), h('button', { class: 'button button-primary', type: 'button', onclick: openHuman }, 'Open Human-GEM'));
    refs.empty.hidden = true;
    el.append(
      refs.canvas,
      h('header', { class: 'view-head floating-head' },
        h('h1', {}, 'Atlas'),
        refs.tiles,
        h('div', { class: 'controls' }, h('button', { class: 'button button-quiet', type: 'button', onclick: openHuman }, 'Open Human-GEM'), refs.status)),
      h('div', { class: 'legend atlas-legend' }, h('span', { class: 'legend-solid' }, 'Identifier match'), h('span', { class: 'legend-dash' }, 'EC number match')),
      refs.empty,
    );
    let drag = null;
    refs.canvas.addEventListener('pointerdown', (e) => {
      refs.canvas.setPointerCapture(e.pointerId);
      drag = { x: e.clientX, y: e.clientY, vx: view.x, vy: view.y };
    });
    refs.canvas.addEventListener('pointermove', (e) => {
      if (drag) {
        view.x = drag.vx + e.clientX - drag.x;
        view.y = drag.vy + e.clientY - drag.y;
        draw();
      }
    });
    refs.canvas.addEventListener('pointerup', () => { drag = null; });
    refs.canvas.addEventListener('wheel', (e) => {
      e.preventDefault();
      const box = refs.canvas.getBoundingClientRect();
      const mx = e.clientX - box.left;
      const my = e.clientY - box.top;
      const factor = Math.exp(-e.deltaY * 0.0015);
      const wx = (mx - view.x) / view.k;
      const wy = (my - view.y) / view.k;
      view.k = Math.min(3, Math.max(0.03, view.k * factor));
      view.x = mx - wx * view.k;
      view.y = my - wy * view.k;
      draw();
    }, { passive: false });
    new ResizeObserver(draw).observe(refs.canvas);
    on('model', () => { data = null; });
  },
  show() {
    load();
  },
};
void emit;
void fmt;
