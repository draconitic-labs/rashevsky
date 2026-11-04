// Shared helpers: element building, calls into the engine, and application state.

export function h(tag, props = {}, ...children) {
  const el = tag === 'svg' || props.svg ? document.createElementNS('http://www.w3.org/2000/svg', tag) : document.createElement(tag);
  for (const [key, value] of Object.entries(props)) {
    if (value === undefined || value === null || value === false || key === 'svg') continue;
    if (key === 'class') el.setAttribute('class', value);
    else if (key === 'dataset') Object.assign(el.dataset, value);
    else if (key === 'style') Object.assign(el.style, value);
    else if (key.startsWith('on')) el.addEventListener(key.slice(2).toLowerCase(), value);
    else if (key === 'value' || key === 'checked' || key === 'disabled' || key === 'selected') el[key] = value;
    else el.setAttribute(key, value === true ? '' : value);
  }
  for (const child of children.flat(Infinity)) {
    if (child === null || child === undefined || child === false) continue;
    el.append(child.nodeType ? child : document.createTextNode(String(child)));
  }
  return el;
}

export function icon(path) {
  const svg = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
  svg.setAttribute('viewBox', '0 0 24 24');
  svg.setAttribute('aria-hidden', 'true');
  svg.setAttribute('class', 'glyph');
  for (const d of [].concat(path)) {
    const p = document.createElementNS('http://www.w3.org/2000/svg', 'path');
    p.setAttribute('d', d);
    svg.append(p);
  }
  return svg;
}

const tauri = window.__TAURI__;

export function invoke(command, args = {}) {
  return tauri.core.invoke(command, args);
}

export async function pickFile(title, extensions = ['xml', 'sbml', 'gz']) {
  return tauri.dialog.open({ title, multiple: false, filters: [{ name: 'SBML model', extensions }] });
}

export async function pickSave(title, defaultPath, extensions) {
  return tauri.dialog.save({ title, defaultPath, filters: [{ name: 'File', extensions }] });
}

// Application state, with listeners by topic.
const listeners = new Map();

export const state = {
  model: null,
  byId: new Map(),
  layout: null,
  method: 'pfba',
  genes: [],          // deleted genes
  reactions: [],      // deleted reactions
  bounds: new Map(),  // medium and manual bounds: id -> [lb, ub]
  solution: null,
  reference: null,    // growth with nothing deleted
  selected: null,
  kinetic: null,
  startup: {},
};

export function on(topic, fn) {
  if (!listeners.has(topic)) listeners.set(topic, new Set());
  listeners.get(topic).add(fn);
}

export function emit(topic, value) {
  for (const fn of listeners.get(topic) ?? []) fn(value);
}

export function perturbation(includeKnockouts = true) {
  return {
    genes: includeKnockouts ? state.genes : [],
    reactions: includeKnockouts ? state.reactions : [],
    bounds: [...state.bounds.entries()].map(([id, [lb, ub]]) => [id, lb, ub]),
  };
}

let solveTicket = 0;
export async function solve() {
  if (!state.model) return;
  const ticket = ++solveTicket;
  const request = { method: state.method, ...perturbation() };
  try {
    const reply = await invoke('solve', { request });
    if (ticket !== solveTicket) return;
    state.solution = reply;
    if (state.genes.length === 0 && state.reactions.length === 0 && reply.growth !== null) {
      state.reference = reply.growth;
    } else if (state.reference === null) {
      const base = await invoke('solve', { request: { method: 'fba', ...perturbation(false) } });
      state.reference = base.growth;
    }
    emit('solution', reply);
  } catch (error) {
    toast(String(error));
  }
}

let referenceTicket = 0;
export async function refreshReference() {
  const ticket = ++referenceTicket;
  const base = await invoke('solve', { request: { method: 'fba', ...perturbation(false) } });
  if (ticket === referenceTicket) state.reference = base.growth;
}

export function setModel(view, curated) {
  state.model = view;
  state.byId = new Map(view.reactions.map((r, i) => [r.id, i]));
  state.genes = [];
  state.reactions = [];
  state.bounds = new Map();
  state.selected = null;
  state.solution = null;
  state.reference = null;
  state.curated = curated;
  emit('model', view);
}

let toastTimer = null;
export function toast(message) {
  const node = document.getElementById('toast');
  if (!node) return;
  node.textContent = message;
  node.hidden = false;
  clearTimeout(toastTimer);
  toastTimer = setTimeout(() => { node.hidden = true; }, 6000);
}

export const fmt = {
  flux(v) {
    if (v === null || v === undefined || Number.isNaN(v)) return 'n/a';
    if (!Number.isFinite(v)) return v > 0 ? '∞' : '−∞';
    const a = Math.abs(v);
    // Below the solver's tolerance a value is zero; printing it would only show rounding noise.
    if (a < 1e-9) return '0';
    const text = a >= 1000 ? a.toFixed(0) : a >= 100 ? a.toFixed(1) : a >= 1 ? a.toFixed(2) : a >= 0.001 ? a.toFixed(3) : a.toExponential(1);
    return (v < 0 ? '−' : '') + text;
  },
  sci(v) {
    if (v === null || v === undefined || Number.isNaN(v)) return 'n/a';
    if (v === 0) return '0';
    return v.toExponential(1).replace('e-', 'e−').replace('e+', 'e');
  },
  percent(v) {
    return `${v >= 0 ? '+' : '−'}${Math.abs(v).toFixed(1)} %`;
  },
};

// The redox fluorescence scale: NADH blue for little flux through to flavin yellow for much.
const RAMP = [[92, 132, 255], [122, 190, 214], [196, 214, 128], [244, 214, 78]];

export function ramp(t) {
  const x = Math.min(1, Math.max(0, t)) * (RAMP.length - 1);
  const i = Math.min(RAMP.length - 2, Math.floor(x));
  const w = x - i;
  const c = [0, 1, 2].map((k) => Math.round(RAMP[i][k] + w * (RAMP[i + 1][k] - RAMP[i][k])));
  return `rgb(${c[0]},${c[1]},${c[2]})`;
}

export const REDUCED_MOTION = window.matchMedia('(prefers-reduced-motion: reduce)').matches;

export function debounce(fn, ms) {
  let timer = null;
  return (...args) => {
    clearTimeout(timer);
    timer = setTimeout(() => fn(...args), ms);
  };
}

/** Polls a background analysis until it finishes. */
export async function follow(id, onProgress) {
  for (;;) {
    const status = await invoke('poll_job', { id });
    onProgress?.(status);
    if (status.finished) {
      if (status.error) throw new Error(status.error);
      return status.result;
    }
    await new Promise((resolve) => setTimeout(resolve, 120));
  }
}
