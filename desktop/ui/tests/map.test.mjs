// The map workspace and its renderer, with the real modules on a stand-in page.
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { test } from 'node:test';

class Element {
  constructor() {
    this.children = [];
    this.style = {};
    this.dataset = {};
    this.hidden = false;
    this.value = '';
    this.textContent = '';
    this.nodeType = 1;
    this.width = 1;
    this.height = 1;
    this.classList = { toggle() {}, add() {}, remove() {}, contains: () => false };
  }

  append(...items) { this.children.push(...items); }

  replaceChildren(...items) { this.children = items; }

  setAttribute() {}

  removeAttribute() {}

  addEventListener() {}

  querySelector() { return new Element(); }

  querySelectorAll() { return []; }

  getBoundingClientRect() { return { left: 0, top: 0, width: 1130, height: 900 }; }

  setPointerCapture() {}

  focus() {}

  matches() { return false; }

  getContext() {
    return new Proxy({}, {
      get: (target, key) => (key === 'createRadialGradient' || key === 'createLinearGradient' ? () => ({ addColorStop() {} }) : key in target ? target[key] : () => {}),
      set: (target, key, value) => { target[key] = value; return true; },
    });
  }
}

const bundled = JSON.parse(readFileSync(new URL('../../../models/ich360/layout.json', import.meta.url)));
const reaction = (id) => ({ id, name: id, lb: 0, ub: 10, equation: 'a → b', rule: '', genes: [], boundary: false });
const view = (id, ids, curated = false) => ({
  id, name: id, curated_layout: curated, reactions: ids.map(reaction), metabolites: [], genes: [], objective: [ids[0]],
  compartments: 1, subsystems: [], medium: [], balance: { balanced: 0, imbalanced: 0, formula_missing: 0 },
});
const models = {
  iCH360: view('iCH360', Object.keys(bundled.reactions), true),
  a: view('a', ['PGI', 'PFK', 'FBA']),
  b: view('b', ['ZZZ1', 'PFK', 'PGI', 'ZZZ2']),
};
let open = null;
function automaticLayout(model) {
  const nodes = {};
  const reactions = {};
  model.reactions.forEach((r, i) => {
    nodes[`m${i}`] = { x: 100 * i, y: 100, label: `m${i}`, external: false };
    nodes[`n${i}`] = { x: 100 * i, y: 300, label: `n${i}`, external: false };
    reactions[r.id] = { s: [`m${i}`], p: [`n${i}`], m: [100 * i, 200], boundary: false };
  });
  return { world: [1200, 800], nodes, regions: [], reactions, source: 'force-directed' };
}
globalThis.window = {
  __TAURI__: {
    core: {
      invoke: async (command, args) => {
        if (command === 'layout') return args.curated ? bundled : automaticLayout(open);
        if (command === 'solve') {
          return { status: 'Optimal', growth: 0.5, fluxes: open.reactions.map((_, i) => i + 1), disabled: [], milliseconds: 1, total_flux: 1 };
        }
        if (command === 'reaction_range') return [0, 1];
        return null;
      },
    },
  },
  matchMedia: () => ({ matches: true }),
  addEventListener() {},
  innerWidth: 1560,
  devicePixelRatio: 1,
};
globalThis.document = {
  createTextNode: (text) => ({ nodeType: 3, text }),
  createElement: () => new Element(),
  createElementNS: () => new Element(),
  getElementById: () => new Element(),
};
globalThis.ResizeObserver = class { observe() {} };
globalThis.requestAnimationFrame = () => 0;

const core = await import('../js/core.js');
const { FluxMap } = await import('../js/fluxmap.js');
let map = null;
const setLayout = FluxMap.prototype.setLayout;
FluxMap.prototype.setLayout = function capture(...args) {
  map = this;
  return setLayout.apply(this, args);
};
const workspace = (await import('../js/views/map.js')).default;
workspace.mount(new Element());

async function openModel(name, path) {
  open = models[name];
  core.setModel(open, path);
  await core.solve();
  await new Promise((resolve) => { setTimeout(resolve, 20); });
}

test('each automatic model gets its own map, with its own fluxes on it', async () => {
  await openModel('iCH360');
  await openModel('a');
  await openModel('b');
  const pgi = map.reactions.find((r) => r.id === 'PGI');
  assert.equal(pgi.modelIndex, 2);
  assert.equal(map.target[pgi.modelIndex], 3);
  assert.deepEqual(map.reactions.map((r) => r.id).sort(), ['PFK', 'PGI', 'ZZZ1', 'ZZZ2']);
});

test('the bundled growth ring goes with the bundled map', async () => {
  await openModel('iCH360');
  assert.equal(map.biomass.id, 'Biomass');
  await openModel('a');
  assert.equal(map.biomass, null);
  map.view = { x: 0, y: 0, k: 0.3 };
  assert.equal(map.pick(2559 * 0.3, 1111 * 0.3), null);
});

test('the growth ring covers no reaction at any zoom', async () => {
  await openModel('iCH360');
  for (const k of [0.1, 0.2, 0.4, 1, 2]) {
    map.view = { x: 0, y: 0, k };
    const radius = 46 * map.ringScale();
    for (const r of map.reactions) {
      for (const [x, y] of r.samples) {
        assert.ok(Math.hypot(x - map.biomass.x, y - map.biomass.y) > radius, `${r.id} at zoom ${k}`);
      }
    }
  }
});

test('a model read from a file remembers that file and no other', async () => {
  await openModel('a', '/models/a.xml');
  assert.equal(core.state.modelPath, '/models/a.xml');
  await openModel('iCH360');
  assert.equal(core.state.modelPath, null);
});

test('rounding noise around zero shows no sign', () => {
  assert.equal(core.fixed(-9.99e-16, 4), '0.0000');
  assert.equal(core.fixed(-0.25, 2), '-0.25');
  assert.equal(core.fixed(0.8656, 4), '0.8656');
});
