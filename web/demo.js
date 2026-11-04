// The flux map of iCH360, solved in the page.
import init, { Core } from './rashevsky_wasm.js';
import { FluxMap } from './fluxmap.js';
import { fixed } from './core.js';

const $ = (id) => document.getElementById(id);
const canvas = $('map');
const tip = $('tip');
const tray = $('tray');
let core;
let reactions = [];
let byId = new Map();
let cut = [];
let reference = null;
const references = new Map();
let map;
let reply = null;

function solve(glucose, oxygen, list) {
  try {
    return JSON.parse(core.solve(glucose, oxygen, list));
  } catch (error) {
    return { growth: null, fluxes: new Array(reactions.length).fill(0), error: String(error.message ?? error) };
  }
}

function compute() {
  // Any value on show is about to change.
  tip.hidden = true;
  const glucose = Number($('glc').value);
  const oxygenValue = Number($('o2').value);
  const oxygen = oxygenValue >= 60 ? 1000 : oxygenValue;
  $('glcOut').textContent = glucose.toFixed(1);
  $('o2Out').textContent = oxygenValue >= 60 ? 'Unlimited' : oxygenValue.toFixed(0);
  reply = solve(glucose, oxygen, cut.join(','));
  // Growth with nothing deleted, kept for each pair of limits.
  const limits = `${glucose} ${oxygen}`;
  if (cut.length === 0) {
    if (!reply.error) references.set(limits, reply.growth);
  } else if (!references.has(limits)) {
    const whole = solve(glucose, oxygen, '');
    if (!whole.error) references.set(limits, whole.growth);
  }
  reference = references.get(limits) ?? null;
  map.setFluxes(reply.fluxes, cut, reply.growth);
  const g = $('growth');
  g.textContent = reply.growth === null ? 'None' : fixed(reply.growth, 4);
  g.classList.toggle('is-dead', !(reply.growth > 1e-9));
  $('delta').textContent = reply.error
    ? `Solver stopped: ${reply.error}`
    : reply.growth === null
    ? 'No steady state'
    : cut.length && reference > 1e-9
      ? `${reply.growth >= reference ? '+' : '−'}${Math.abs((100 * (reply.growth - reference)) / reference).toFixed(1)} % from unperturbed`
      : '';
  renderTray();
}

function toggle(id) {
  const at = cut.indexOf(id);
  if (at >= 0) cut.splice(at, 1);
  else {
    cut.push(id);
    map.flash(id);
  }
  compute();
}

function renderTray() {
  tray.replaceChildren();
  tray.hidden = cut.length === 0;
  for (const id of cut) {
    const b = document.createElement('button');
    b.className = 'chip';
    b.type = 'button';
    b.textContent = `${id} ×`;
    b.title = `Restore ${byId.get(id)?.name ?? id}`;
    b.addEventListener('click', () => toggle(id));
    tray.append(b);
  }
  if (cut.length) {
    const all = document.createElement('button');
    all.className = 'restore';
    all.type = 'button';
    all.textContent = 'Restore all';
    all.addEventListener('click', () => { cut = []; compute(); });
    tray.append(all);
  }
}

async function main() {
  await init();
  core = new Core();
  reactions = JSON.parse(core.reactions());
  byId = new Map(reactions.map((r, i) => [r.id, { ...r, index: i }]));
  map = new FluxMap(canvas, {
    indexOfId: (id) => byId.get(id)?.index,
    referenceGrowth: () => reference,
    insets: () => (window.innerWidth > 700 ? { left: 240, top: 20, right: 20, bottom: 20 } : { left: 10, top: 110, right: 10, bottom: 130 }),
    onSelect: (hit) => { if (hit?.kind === 'reaction' && hit.id !== map.biomass?.id) toggle(hit.id); },
    onHover: (hit, x, y) => {
      if (!hit) {
        tip.hidden = true;
        return;
      }
      tip.hidden = false;
      tip.style.left = `${x + 14}px`;
      tip.style.top = `${y + 14}px`;
      const reaction = byId.get(hit.id);
      const name = document.createElement('strong');
      name.textContent = reaction?.name ?? hit.id;
      const value = document.createElement('span');
      value.textContent = reaction ? `${hit.id}, flux ${fixed(reply?.fluxes[reaction.index] ?? 0, 2)}` : hit.id;
      tip.replaceChildren(name, value);
    },
  });
  const layout = JSON.parse(core.layout());
  map.resize();
  map.setLayout(layout, reactions.map((r) => r.id));
  map.start();
  compute();
  $('glc').addEventListener('input', compute);
  $('o2').addEventListener('input', compute);
  document.fonts?.ready.then(() => map.redraw());
}

main();
