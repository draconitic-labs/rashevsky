// The flux map workspace: the network, what is deleted from it, and what it does in response.

import { debounce, emit, fmt, h, icon, invoke, on, perturbation, pickSave, refreshReference, solve, state, toast } from '../core.js';
import { FluxMap } from '../fluxmap.js';

const METHODS = [
  ['pfba', 'pFBA', 'Optimal growth with the least total flux'],
  ['fba', 'FBA', 'Any flux distribution with optimal growth'],
  ['moma', 'MOMA', 'Closest fluxes to the unperturbed cell, by Euclidean distance'],
  ['lmoma', 'Linear MOMA', 'Closest fluxes to the unperturbed cell, by total change'],
  ['room', 'ROOM', 'Fewest significant flux changes from the unperturbed cell'],
  ['loopless', 'Loopless', 'Optimal growth with no flux around a closed loop'],
];

let map;
let refs = {};
let visible = false;

function reactionOf(id) {
  return state.model?.reactions[state.byId.get(id)];
}

function nameOf(id) {
  const r = reactionOf(id);
  return r?.name || id;
}

function exchangeLabel(id) {
  const r = reactionOf(id);
  return (r?.name || id).replace(/\s*exchange\s*$/i, '').replace(/^exchange of\s*/i, '').replace(/^D-|^L-/, '');
}

function toggleReaction(id) {
  const at = state.reactions.indexOf(id);
  if (at >= 0) state.reactions.splice(at, 1);
  else {
    state.reactions.push(id);
    map.flash(id);
  }
  changed();
}

function toggleGene(gene) {
  const at = state.genes.indexOf(gene);
  if (at >= 0) state.genes.splice(at, 1);
  else state.genes.push(gene);
  changed();
}

function changed() {
  renderTray();
  renderInspector();
  solve();
}

function growthReadout(reply) {
  const value = reply?.growth;
  const box = refs.growth;
  box.value.textContent = value === null || value === undefined ? 'None' : value.toFixed(4);
  box.value.classList.toggle('is-dead', !(value > 1e-9));
  let delta = '';
  if (reply?.status !== 'Optimal') delta = 'No steady state';
  else if ((state.genes.length || state.reactions.length) && state.reference > 1e-9) {
    delta = `${fmt.percent((100 * (value - state.reference)) / state.reference)} from unperturbed`;
  } else delta = reply ? `${(reply.total_flux ?? 0).toFixed(0)} total flux` : '';
  box.delta.textContent = delta;
  box.time.textContent = reply ? `${reply.milliseconds.toFixed(reply.milliseconds < 10 ? 1 : 0)} ms` : '';
}

function renderTray() {
  const tray = refs.tray;
  tray.replaceChildren();
  const items = [
    ...state.genes.map((g) => ({ label: g, title: `Gene ${g}`, undo: () => toggleGene(g) })),
    ...state.reactions.map((r) => ({ label: r, title: nameOf(r), undo: () => toggleReaction(r) })),
  ];
  tray.hidden = items.length === 0;
  if (!items.length) return;
  tray.append(h('span', { class: 'tray-label' }, 'Deleted'));
  for (const item of items) {
    tray.append(h('button', { class: 'chip chip-cut', type: 'button', title: `Restore ${item.title}`, onclick: item.undo }, item.label, h('span', { class: 'chip-x', 'aria-hidden': 'true' }, '×')));
  }
  const silent = (state.solution?.disabled.length ?? 0) - state.reactions.length;
  if (state.genes.length) tray.append(h('span', { class: 'tray-note' }, `${state.solution?.disabled.length ?? 0} reactions off`));
  void silent;
  tray.append(h('button', { class: 'link', type: 'button', onclick: () => { state.genes = []; state.reactions = []; changed(); } }, 'Restore all'));
}

async function renderInspector() {
  const pane = refs.inspector;
  const target = state.selected;
  pane.classList.toggle('is-open', !!target);
  pane.setAttribute('aria-hidden', target ? 'false' : 'true');
  pane.replaceChildren();
  if (!target) return;
  const close = h('button', { class: 'icon-button', type: 'button', title: 'Close (Esc)', 'aria-label': 'Close', onclick: () => select(null) }, icon('M6 6l12 12M18 6L6 18'));
  if (target.kind === 'metabolite') return renderMetabolite(pane, target.id, close);
  const j = state.byId.get(target.id);
  const r = state.model.reactions[j];
  const flux = state.solution?.fluxes[j] ?? 0;
  const off = state.reactions.includes(r.id);
  const silenced = !off && state.solution?.disabled.includes(r.id);
  const manual = state.bounds.get(r.id);
  const lb = manual ? manual[0] : r.lb;
  const ub = manual ? manual[1] : r.ub;
  const rangeOut = h('span', { class: 'range-text' }, 'Working out range');
  const rangeBar = h('div', { class: 'range-bar' }, h('div', { class: 'range-fill' }), h('div', { class: 'range-now' }));
  pane.append(
    h('header', { class: 'inspector-head' },
      h('div', {}, h('code', { class: 'ident' }, r.id), h('h2', {}, r.name || r.id)),
      close),
    h('p', { class: 'equation' }, r.equation),
    h('dl', { class: 'facts' },
      h('div', {}, h('dt', {}, 'Flux'), h('dd', { class: 'big-number' }, fmt.flux(flux))),
      r.subsystem && h('div', {}, h('dt', {}, 'Subsystem'), h('dd', {}, r.subsystem))),
    h('div', { class: 'block' }, h('h3', {}, 'Range at optimal growth'), rangeOut, rangeBar),
    h('div', { class: 'block' },
      h('h3', {}, 'Bounds'),
      h('div', { class: 'bounds-row' },
        h('label', {}, 'Lower', numberBox(lb, (v) => setBounds(r.id, v, null))),
        h('label', {}, 'Upper', numberBox(ub, (v) => setBounds(r.id, null, v))),
        manual && h('button', { class: 'link', type: 'button', onclick: () => { state.bounds.delete(r.id); refreshReference().then(changed); } }, 'Reset'))),
    r.genes.length > 0 && h('div', { class: 'block' },
      h('h3', {}, 'Genes'),
      h('p', { class: 'rule' }, r.rule),
      h('div', { class: 'chips' }, r.genes.map((g) => h('button', {
        class: `chip ${state.genes.includes(g) ? 'chip-cut' : ''}`,
        type: 'button',
        'aria-pressed': String(state.genes.includes(g)),
        title: state.genes.includes(g) ? `Restore ${g}` : `Delete ${g}`,
        onclick: () => toggleGene(g),
      }, g)))),
    h('div', { class: 'inspector-actions' },
      h('button', { class: off ? 'button' : 'button button-cut', type: 'button', onclick: () => toggleReaction(r.id) }, off ? 'Restore reaction' : 'Delete reaction'),
      silenced && h('span', { class: 'hint' }, 'Off through a gene deletion')),
  );
  try {
    const [low, high] = await invoke('reaction_range', { id: r.id, fraction: 1, perturbation: perturbation() });
    rangeOut.textContent = `${fmt.flux(low)} to ${fmt.flux(high)}`;
    const span = Math.max(Math.abs(lb), Math.abs(ub), Math.abs(low), Math.abs(high), 1e-9);
    const pos = (v) => `${50 + (50 * v) / span}%`;
    rangeBar.querySelector('.range-fill').style.left = pos(low);
    rangeBar.querySelector('.range-fill').style.width = `${(50 * (high - low)) / span}%`;
    rangeBar.querySelector('.range-now').style.left = pos(flux);
  } catch (error) {
    rangeOut.textContent = String(error);
  }
}

function renderMetabolite(pane, id, close) {
  const met = state.model.metabolites.find((m) => m.id === id);
  const rows = [];
  state.model.reactions.forEach((r, j) => {
    const flux = state.solution?.fluxes[j] ?? 0;
    const sides = r.equation.split(/\s[→⇌←]\s/);
    const inLeft = sides[0]?.split(' + ').some((t) => t.split(' ').pop() === id);
    const inRight = sides[1]?.split(' + ').some((t) => t.split(' ').pop() === id);
    if (inLeft || inRight) rows.push({ r, flux, makes: inRight ? flux > 0 : flux < 0 });
  });
  rows.sort((a, b) => Math.abs(b.flux) - Math.abs(a.flux));
  const active = rows.filter((x) => Math.abs(x.flux) > 1e-7);
  const line = ({ r, flux, makes }) => h('button', { class: 'met-row', type: 'button', onclick: () => select({ kind: 'reaction', id: r.id }) },
    h('code', { class: 'ident' }, r.id), h('span', { class: 'met-flux' }, `${makes ? 'makes' : 'uses'} ${fmt.flux(Math.abs(flux))}`));
  pane.append(
    h('header', { class: 'inspector-head' }, h('div', {}, h('code', { class: 'ident' }, id), h('h2', {}, met?.name || id)), close),
    h('dl', { class: 'facts' },
      met?.formula && h('div', {}, h('dt', {}, 'Formula'), h('dd', {}, met.formula)),
      met?.charge !== null && met?.charge !== undefined && h('div', {}, h('dt', {}, 'Charge'), h('dd', {}, String(met.charge))),
      h('div', {}, h('dt', {}, 'Reactions'), h('dd', {}, String(rows.length)))),
    h('div', { class: 'block' }, h('h3', {}, 'Carrying flux now'), active.length ? h('div', { class: 'met-list' }, active.map(line)) : h('p', { class: 'hint' }, 'None')),
    rows.length > active.length && h('div', { class: 'block' }, h('h3', {}, 'Idle'), h('div', { class: 'met-list' }, rows.filter((x) => Math.abs(x.flux) <= 1e-7).slice(0, 30).map(line))),
  );
}

function numberBox(value, apply) {
  return h('input', {
    type: 'number',
    step: 'any',
    value: Number.isFinite(value) ? value : '',
    placeholder: value > 0 ? '∞' : '−∞',
    onchange: (e) => {
      const text = e.target.value.trim();
      apply(text === '' ? (value > 0 ? Infinity : -Infinity) : Number(text));
    },
  });
}

async function setBounds(id, lb, ub) {
  const r = reactionOf(id);
  const current = state.bounds.get(id) ?? [r.lb, r.ub];
  state.bounds.set(id, [lb ?? current[0], ub ?? current[1]]);
  await refreshReference();
  changed();
}

export function select(target) {
  state.selected = target;
  map.select(target);
  renderInspector();
  if (target?.kind === 'reaction') {
    const r = map.reactions[map.indexOf.get(state.byId.get(target.id))];
    if (r && (r.mx < 0)) map.focus(r.mx, r.my);
  }
}

function renderMedium() {
  const list = refs.medium;
  list.replaceChildren();
  const model = state.model;
  if (!model) return;
  // Carbon sources first, then oxygen and nitrogen. Open uptakes of anything else, and water, protons,
  // ions and carbon dioxide, would only crowd the panel.
  const rank = (m) => {
    const label = exchangeLabel(m.id).toLowerCase();
    const id = m.id.toLowerCase();
    if (/co2|carbon dioxide/.test(`${id} ${label}`)) return null;
    if (/(^|_)o2(_|$)|^oxygen$/.test(id) || /^(o2|oxygen)$/.test(label)) return 1;
    if (/nh4|ammoni/.test(`${id} ${label}`)) return 2;
    if (/^(d-)?glucose$/.test(label) || /^ex_glc/.test(id)) return 0;
    if (m.limit >= 200) return null;
    if (/h2o|water|\bh\+|proton|phosphate|sulfate|sulphate|^pi$|^fe|iron|mg2|ca2|^k\b|^na|^cl|cobalt|copper|zinc|mn2|\bmo\b|ni2|tungstate|selen/.test(`${id} ${label}`)) return null;
    return 0;
  };
  const shown = model.medium
    .map((m) => ({ m, r: rank(m) }))
    .filter((x) => x.r !== null)
    .sort((a, b) => a.r - b.r)
    .slice(0, 4)
    .map((x) => x.m);
  const add = (m) => {
    const r = reactionOf(m.id);
    if (!r) return;
    const finite = m.limit < 200;
    const top = finite ? Math.max(30, m.limit * 2) : 60;
    const current = state.bounds.get(m.id);
    const value = current ? -current[0] : m.limit;
    const out = h('output', {}, value >= top ? 'Unlimited' : value.toFixed(1));
    const slider = h('input', {
      type: 'range', min: 0, max: top, step: 0.5, value: Math.min(value, top),
      'aria-label': `${exchangeLabel(m.id)} uptake`,
      oninput: (e) => {
        const v = Number(e.target.value);
        const unlimited = !finite && v >= top;
        out.textContent = unlimited ? 'Unlimited' : v.toFixed(1);
        state.bounds.set(m.id, [unlimited ? -m.limit : -v, r.ub]);
        mediumChanged();
      },
    });
    list.append(h('label', { class: 'slider' }, h('span', { class: 'slider-name' }, exchangeLabel(m.id), out), slider));
  };
  shown.forEach(add);
  refs.mediumPanel.hidden = shown.length === 0;
}

const mediumChanged = debounce(async () => {
  await refreshReference();
  solve();
}, 60);

function renderFind(query) {
  const box = refs.results;
  box.replaceChildren();
  const q = query.trim().toLowerCase();
  if (!q || !state.model) {
    box.hidden = true;
    return;
  }
  const hits = [];
  for (const r of state.model.reactions) {
    if (r.id.toLowerCase().includes(q) || r.name.toLowerCase().includes(q)) hits.push({ kind: 'reaction', id: r.id, text: r.name });
    if (hits.length > 12) break;
  }
  for (const m of state.model.metabolites) {
    if (hits.length > 18) break;
    if (m.id.toLowerCase().includes(q) || m.name.toLowerCase().includes(q)) hits.push({ kind: 'metabolite', id: m.id, text: m.name });
  }
  for (const g of state.model.genes) {
    if (hits.length > 22) break;
    if (g.id.toLowerCase().includes(q) || (g.name || '').toLowerCase().includes(q)) hits.push({ kind: 'gene', id: g.id, text: g.name || 'Gene' });
  }
  box.hidden = hits.length === 0;
  for (const hit of hits) {
    box.append(h('button', {
      class: 'find-hit',
      type: 'button',
      onclick: () => {
        refs.find.value = '';
        box.hidden = true;
        if (hit.kind === 'gene') toggleGene(hit.id);
        else {
          select(hit);
          const n = hit.kind === 'metabolite' ? map.nodeById.get(hit.id) : map.reactions[map.indexOf.get(state.byId.get(hit.id))];
          if (n) map.focus(n.x ?? n.mx, n.y ?? n.my);
        }
      },
    }, h('code', { class: 'ident' }, hit.id), h('span', {}, hit.kind === 'gene' ? (state.genes.includes(hit.id) ? `Restore gene ${hit.id}` : `Delete gene ${hit.text === 'Gene' ? '' : hit.text}`) : hit.text)));
  }
}

async function loadLayout() {
  const model = state.model;
  if (!model) return;
  const curated = !!state.curated;
  refs.viewSelect.hidden = curated;
  if (curated) {
    const layout = await invoke('layout', { curated: true });
    map.setLayout(layout, model.reactions.map((r) => r.id));
  } else {
    refs.viewSelect.replaceChildren(
      h('option', { value: '' }, 'Reactions with most flux'),
      ...model.subsystems.map((s) => h('option', { value: s }, s)),
    );
    const layout = await invoke('layout', { curated: false, fluxes: state.solution?.fluxes ?? null, count: 110 });
    map.setLayout(layout, model.reactions.map((r) => r.id));
  }
  if (state.solution) {
    map.setFluxes(state.solution.fluxes, state.solution.disabled, state.solution.growth);
  }
}

export default {
  id: 'map',
  label: 'Map',
  icon: ['M5 6a2 2 0 1 0 0 .01', 'M19 6a2 2 0 1 0 0 .01', 'M12 18a2 2 0 1 0 0 .01', 'M7 6h10', 'M6 8l5 8', 'M18 8l-5 8'],
  mount(el) {
    const canvas = h('canvas', { class: 'flux-canvas', tabindex: 0, 'aria-label': 'Flux map. Use Find to reach a reaction by name.' });
    const tooltip = h('div', { class: 'tooltip', hidden: true });
    tooltip.hidden = true;
    refs.inspector = h('aside', { class: 'inspector', 'aria-label': 'Selection', 'aria-hidden': 'true' });
    refs.tray = h('div', { class: 'tray', hidden: true });
    refs.tray.hidden = true;
    refs.results = h('div', { class: 'find-results', hidden: true });
    refs.results.hidden = true;
    refs.find = h('input', { class: 'find', type: 'search', placeholder: 'Find a reaction, metabolite or gene', 'aria-label': 'Find', oninput: (e) => renderFind(e.target.value) });
    refs.viewSelect = h('select', {
      class: 'select', 'aria-label': 'Part of the model to draw', hidden: true,
      onchange: async (e) => {
        const subsystem = e.target.value || null;
        const layout = await invoke('layout', { curated: false, subsystem, fluxes: state.solution?.fluxes ?? null, count: 110 });
        map.setLayout(layout, state.model.reactions.map((r) => r.id));
        if (state.solution) map.setFluxes(state.solution.fluxes, state.solution.disabled, state.solution.growth);
      },
    });
    refs.viewSelect.hidden = true;
    const methods = h('div', { class: 'segmented', role: 'radiogroup', 'aria-label': 'Flux method' },
      METHODS.map(([id, label, title]) => h('button', {
        class: 'segment', type: 'button', role: 'radio', title, dataset: { method: id },
        'aria-checked': String(state.method === id),
        onclick: (e) => {
          state.method = id;
          e.currentTarget.parentElement.querySelectorAll('.segment').forEach((b) => b.setAttribute('aria-checked', String(b.dataset.method === id)));
          solve();
        },
      }, label)));
    refs.growth = {
      value: h('span', { class: 'growth-value' }, '0'),
      delta: h('span', { class: 'growth-delta' }),
      time: h('span', { class: 'growth-time' }),
    };
    refs.title = h('h1', { class: 'model-title' }, '');
    refs.growthLabel = h('span', { class: 'growth-label' }, 'Growth, 1/h');
    refs.fps = h('span', { class: 'fps', hidden: true });
    refs.medium = h('div', { class: 'medium-list' });
    refs.mediumPanel = h('div', { class: 'floating medium' }, h('h3', {}, 'Uptake limits'), refs.medium);
    const zoom = h('div', { class: 'zoom' },
      h('button', { class: 'icon-button', type: 'button', title: 'Zoom in', 'aria-label': 'Zoom in', onclick: () => map.zoomBy(1.25) }, icon('M12 5v14M5 12h14')),
      h('button', { class: 'icon-button', type: 'button', title: 'Zoom out', 'aria-label': 'Zoom out', onclick: () => map.zoomBy(0.8) }, icon('M5 12h14')),
      h('button', { class: 'icon-button', type: 'button', title: 'Fit to window', 'aria-label': 'Fit to window', onclick: () => map.fit() }, icon('M4 9V4h5M20 9V4h-5M4 15v5h5M20 15v5h-5')));
    const exporter = h('details', { class: 'menu' },
      h('summary', { class: 'button button-quiet' }, 'Export'),
      h('div', { class: 'menu-body' },
        h('button', { type: 'button', onclick: () => exportFluxes('csv') }, 'Fluxes as CSV'),
        h('button', { type: 'button', onclick: () => exportFluxes('json') }, 'Fluxes for Escher'),
        h('button', { type: 'button', onclick: exportEscher }, 'Escher map with data')));
    el.append(
      canvas,
      h('header', { class: 'topbar' },
        h('div', { class: 'topbar-left' }, refs.title, refs.viewSelect),
        methods,
        h('div', { class: 'find-box' }, refs.find, refs.results),
        exporter,
        h('div', { class: 'growth', title: 'Value of the objective' },
          refs.growthLabel,
          refs.growth.value, refs.growth.delta, refs.growth.time)),
      refs.mediumPanel,
      zoom,
      refs.tray,
      refs.inspector,
      tooltip,
      refs.fps,
      legend(),
    );
    map = new FluxMap(canvas, {
      indexOfId: (id) => state.byId.get(id),
      insets: () => ({ left: refs.mediumPanel.hidden ? 24 : 240, top: 58, right: 62, bottom: 22 }),
      referenceGrowth: () => state.reference,
      onSelect: (hit) => select(hit),
      onCut: (id) => toggleReaction(id),
      onFps: (fps) => {
        if (state.startup.qa) invoke('qa_log', { line: `fps ${fps.toFixed(1)} particles ${map.particles.length}` });
      },
      onHover: (hit, x, y) => {
        if (!hit) {
          tooltip.hidden = true;
          return;
        }
        tooltip.hidden = false;
        const box = el.getBoundingClientRect();
        tooltip.style.left = `${x - box.left + 14}px`;
        tooltip.style.top = `${y - box.top + 14}px`;
        if (hit.kind === 'reaction') {
          const j = state.byId.get(hit.id);
          tooltip.replaceChildren(h('strong', {}, nameOf(hit.id)), h('span', {}, `${hit.id}, flux ${fmt.flux(state.solution?.fluxes[j] ?? 0)}`));
        } else {
          const met = state.model.metabolites.find((m) => m.id === hit.id);
          tooltip.replaceChildren(h('strong', {}, met?.name || hit.id), h('span', {}, hit.id));
        }
      },
    });
    on('model', async (model) => {
      refs.title.textContent = (model.name && model.name.length <= 30) ? model.name : model.id;
      refs.title.title = model.name || model.id;
      const objective = reactionOf(model.objective[0]);
      refs.growthLabel.textContent = /biomass|growth/i.test(`${objective?.id} ${objective?.name}`) ? 'Growth, 1/h' : 'Objective';
      renderMedium();
      renderTray();
      renderInspector();
      if (state.curated) await loadLayout();
    });
    on('solution', async (reply) => {
      if (!state.curated && !map.layout?.forModel) {
        await loadLayout();
        if (map.layout) map.layout.forModel = true;
      }
      map.setFluxes(reply.fluxes, reply.disabled, reply.growth);
      growthReadout(reply);
      tooltip.hidden = true;
      renderTray();
      if (state.selected) renderInspector();
    });
    on('fonts', () => map.redraw());
    on('apply-design', ({ genes, reactions }) => {
      state.genes = [...genes];
      state.reactions = [...reactions];
      emit('navigate', 'map');
      changed();
    });
    window.addEventListener('keydown', (e) => {
      if (!visible) return;
      if (e.key === 'Escape') {
        if (refs.find.value) {
          refs.find.value = '';
          refs.results.hidden = true;
        } else select(null);
      }
      if (e.target.matches('input, select, textarea')) return;
      if (e.key === '/') {
        e.preventDefault();
        refs.find.focus();
      }
      if (e.key.toLowerCase() === 'k' && state.selected?.kind === 'reaction') toggleReaction(state.selected.id);
      if (e.key === '0') map.fit();
    });
  },
  show() {
    visible = true;
    map.start();
    requestAnimationFrame(() => map.resize());
  },
  hide() {
    visible = false;
    map.stop();
  },
};

function legend() {
  const swatch = h('canvas', { class: 'legend-ramp', width: 120, height: 8 });
  const g = swatch.getContext('2d');
  const gradient = g.createLinearGradient(0, 0, 120, 0);
  gradient.addColorStop(0, 'rgb(92,132,255)');
  gradient.addColorStop(0.33, 'rgb(122,190,214)');
  gradient.addColorStop(0.66, 'rgb(196,214,128)');
  gradient.addColorStop(1, 'rgb(244,214,78)');
  g.fillStyle = gradient;
  g.fillRect(0, 0, 120, 8);
  return h('div', { class: 'legend' }, h('span', {}, 'Less flux'), swatch, h('span', {}, 'More'));
}

async function exportFluxes(kind) {
  const path = await pickSave('Export fluxes', `${state.model.id}-fluxes.${kind}`, [kind]);
  if (!path) return;
  try {
    await invoke('export_fluxes', { path, fluxes: state.solution.fluxes });
    toast(`Saved ${path}`);
  } catch (error) {
    toast(String(error));
  }
}

async function exportEscher() {
  const path = await pickSave('Export Escher map', `${state.model.id}-map.json`, ['json']);
  if (!path) return;
  try {
    await invoke('export_escher', { path, fluxes: state.solution.fluxes, subsystem: refs.viewSelect.value || null });
    toast(`Saved ${path}`);
  } catch (error) {
    toast(String(error));
  }
}
