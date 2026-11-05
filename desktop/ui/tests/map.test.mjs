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
    this.attributes = {};
    this.listeners = {};
  }

  append(...items) { this.children.push(...items); }

  replaceChildren(...items) { this.children = items; }

  setAttribute(key, value) { this.attributes[key] = value; }

  removeAttribute() {}

  addEventListener(name, fn) { this.listeners[name] = fn; }

  querySelector() { return new Element(); }

  querySelectorAll() { return []; }

  getBoundingClientRect() { return { left: 0, top: 0, width: 1130, height: 900 }; }

  setPointerCapture() {}

  focus() {}

  matches() { return false; }

  getContext() {
    return new Proxy({}, {
      get: (target, key) => {
        if (key === 'createRadialGradient' || key === 'createLinearGradient') return () => ({ addColorStop() {} });
        if (key === 'fillText') return (text) => drawn.push(text);
        return key in target ? target[key] : () => {};
      },
      set: (target, key, value) => { target[key] = value; return true; },
    });
  }
}

const drawn = [];
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
  unsolvable: { ...view('unsolvable', ['V1', 'V2']), objective: [] },
  parts: { ...view('parts', ['PGI', 'PFK']), subsystems: ['Glycolysis'] },
};
let open = null;
// The layout of this model waits for this promise, as it does on a large model.
let gatedModel = null;
let layoutGate = null;
let layoutFailure = null;
let saveAsked = false;
const calls = [];
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
        calls.push({ command, args });
        if (command === 'layout') {
          const model = open;
          if (model === gatedModel) await layoutGate;
          if (layoutFailure) throw layoutFailure;
          return args.curated ? bundled : automaticLayout(model);
        }
        if (command === 'solve') {
          if (open === models.unsolvable) throw new Error('the model has no objective');
          const shift = 100 * args.request.reactions.length;
          return { status: open.status ?? 'Optimal', growth: 0.5, fluxes: open.reactions.map((_, i) => shift + i + 1), disabled: [], milliseconds: 1, total_flux: 1 };
        }
        if (command === 'bundled_terms') return [];
        if (command === 'open_file') return models.a;
        if (command === 'download_model') return `/data/${args.name}.xml`;
        if (command === 'reaction_range') return [0, 1];
        return null;
      },
    },
    dialog: {
      open: async () => '/models/picked.xml',
      save: async () => {
        saveAsked = true;
        return '/out/map.json';
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
const page = new Element();
workspace.mount(page);
const modelsPage = new Element();
(await import('../js/views/models.js')).default.mount(modelsPage);
const atlasPage = new Element();
(await import('../js/views/atlas.js')).default.mount(atlasPage);

function words(el) {
  return (el.children ?? []).map((c) => (c.nodeType === 3 ? c.text : words(c))).join(' ');
}

function find(root, wanted) {
  if (wanted(root)) return root;
  for (const child of root.children ?? []) {
    if (child.nodeType !== 1) continue;
    const hit = find(child, wanted);
    if (hit) return hit;
  }
  return null;
}

const withClass = (name) => find(page, (el) => el.attributes.class === name);
const button = (root, text) => find(root, (el) => el.listeners.click && words(el).includes(text));
const settle = () => new Promise((resolve) => { setTimeout(resolve, 20); });

async function openModel(name, path) {
  open = models[name];
  core.setModel(open, path);
  await core.solve();
  await settle();
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
  assert.equal(core.fmt.percent(-6.3e-13), '0.0 %');
  assert.equal(core.fmt.percent(-1.24), '−1.2 %');
  assert.equal(core.fmt.percent(2.5), '+2.5 %');
});

test('a model that cannot be solved takes the last map off view and draws its own', async () => {
  await openModel('iCH360');
  map.select({ kind: 'reaction', id: 'PGI' });
  await openModel('unsolvable');
  assert.equal(map.biomass, null);
  assert.equal(map.selected, null);
  assert.deepEqual(map.reactions.map((r) => r.id).sort(), ['V1', 'V2']);
  assert.equal(core.state.solution, null);
  assert.equal(withClass('growth-value').textContent, 'None');
  assert.equal(withClass('growth-delta').textContent, 'Not solved');
  // The animation runs on what it was given.
  map.last = 0;
  map.frame(16);
  assert.deepEqual(map.shown, [0, 0]);
});

test('a failed solve after deletions keeps them on view', async () => {
  await openModel('a');
  core.state.reactions = ['PFK'];
  core.state.genes = ['g1'];
  const solvable = open;
  open = models.unsolvable;
  await core.solve();
  await settle();
  assert.equal(core.state.solution, null);
  assert.ok(map.disabled.has('PFK'));
  assert.ok(!words(withClass('tray')).includes('reactions off'));
  open = solvable;
  core.state.reactions = [];
  core.state.genes = [];
});

test('after a failed first solve, the first reply ranks the map by its fluxes', async () => {
  await openModel('unsolvable');
  // The same model, solvable from now on.
  open = { ...models.a, reactions: models.unsolvable.reactions };
  calls.length = 0;
  await core.solve();
  await settle();
  assert.ok(calls.some((c) => c.command === 'layout' && c.args.fluxes !== null));
  assert.equal(map.target.length, 2);
});

test('one failed layout does not stop the next reply reaching the map', async () => {
  await openModel('iCH360');
  layoutFailure = new Error('the layout failed');
  open = models.a;
  core.setModel(open);
  await core.solve();
  await settle();
  layoutFailure = null;
  core.state.reactions = ['FBA'];
  await core.solve();
  await settle();
  assert.equal(withClass('growth-value').textContent, '0.5000');
  core.state.reactions = [];
});

test('another model opening hides nothing behind a failed reply', async () => {
  await openModel('iCH360');
  let release;
  layoutGate = new Promise((resolve) => { release = resolve; });
  gatedModel = models.unsolvable;
  open = models.unsolvable;
  core.setModel(open);
  await core.solve();
  await settle();
  gatedModel = null;
  open = models.b;
  core.setModel(open);
  release();
  await settle();
  // Its own reply has not come in yet, so nothing is claimed about it.
  assert.equal(withClass('growth-delta').textContent, '');
  await core.solve();
  await settle();
  assert.equal(withClass('growth-value').textContent, '0.5000');
});

test('an export before any result says so, before asking where to save', async () => {
  await openModel('unsolvable');
  calls.length = 0;
  saveAsked = false;
  await button(page, 'Fluxes as CSV').listeners.click();
  await button(page, 'Escher map with data').listeners.click();
  assert.equal(saveAsked, false);
  assert.ok(!calls.some((c) => c.command.startsWith('export')));
});

test('a reply that is not optimal says what happened', async () => {
  const expected = { Infeasible: 'No steady state', Unbounded: 'No limit on the objective', Unknown: 'Solver stopped: Unknown' };
  for (const [status, text] of Object.entries(expected)) {
    models.a.status = status;
    await openModel('a');
    assert.equal(withClass('growth-delta').textContent, text);
  }
  delete models.a.status;
});

test('a slow layout does not put an older reply on the map', async () => {
  let release;
  layoutGate = new Promise((resolve) => { release = resolve; });
  gatedModel = models.a;
  open = models.a;
  core.setModel(open);
  const first = core.solve();
  await settle();
  core.state.reactions = ['FBA'];
  const second = core.solve();
  await Promise.all([first, second]);
  await settle();
  gatedModel = null;
  release();
  await settle();
  assert.equal(core.state.solution.fluxes[0], 101);
  assert.equal(map.target[0], 101);
  core.state.reactions = [];
});

test('another model opened during a layout keeps the map to itself', async () => {
  let release;
  layoutGate = new Promise((resolve) => { release = resolve; });
  gatedModel = models.a;
  open = models.a;
  core.setModel(open);
  await core.solve();
  await settle();
  // The second model's layout comes back first, and the first model's after it.
  open = models.b;
  core.setModel(open);
  await core.solve();
  await settle();
  gatedModel = null;
  release();
  await settle();
  assert.deepEqual(map.reactions.map((r) => r.id).sort(), ['PFK', 'PGI', 'ZZZ1', 'ZZZ2']);
  assert.deepEqual(map.target, [1, 2, 3, 4]);
});

test('a reply to the last model is dropped when another has opened', async () => {
  await openModel('a');
  const late = core.solve();
  open = models.b;
  core.setModel(open);
  await late;
  assert.equal(core.state.solution, null);
  await core.solve();
  await settle();
  assert.equal(core.state.solution.fluxes.length, 4);
});

test('an export asks for no part of the last model', async () => {
  await openModel('parts');
  find(page, (el) => el.listeners.change).value = 'Glycolysis';
  await openModel('iCH360');
  calls.length = 0;
  await button(page, 'Escher map with data').listeners.click();
  assert.equal(calls.find((c) => c.command === 'export_escher').args.subsystem, null);
});

test('a frame draws both signs of flux, the ring, and no sign on noise', async () => {
  await openModel('iCH360');
  drawn.length = 0;
  map.view = { x: 0, y: 0, k: 2 };
  const fluxes = map.target.map((_, i) => (i % 2 ? -3 : 3));
  fluxes[map.reactions[0].modelIndex] = -0.00116;
  map.setFluxes(fluxes, [], 0.5);
  map.shown = fluxes.slice();
  map.transition = 1;
  map.layerDirty = true;
  map.last = 0;
  map.frame(16);
  assert.ok(drawn.includes(`${map.reactions[0].id}  0.00`), drawn.slice(0, 3).join());
  assert.ok(!drawn.some((text) => String(text).includes('-0.00')));
  assert.ok(drawn.includes('biomass'));
});

test('the Models workspace passes on the file a model came from', async () => {
  await openModel('iCH360');
  open = models.a;
  await button(modelsPage, 'Open SBML file').listeners.click();
  await settle();
  assert.equal(core.state.modelPath, '/models/picked.xml');
  await button(modelsPage, 'E. coli core').listeners.click();
  await settle();
  assert.equal(core.state.modelPath, '/data/e_coli_core.xml');
});

test('the Atlas workspace passes on the file a model came from', async () => {
  await openModel('iCH360');
  open = models.a;
  await button(atlasPage, 'Open Hu').listeners.click();
  await settle();
  assert.equal(core.state.modelPath, '/data/Human-GEM.xml');
});
