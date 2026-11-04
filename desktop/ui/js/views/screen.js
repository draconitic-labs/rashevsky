// Knockout screens: every single or double deletion, ranked for a production target.

import { conditionsKey, deletionNote, emit, fmt, follow, h, invoke, on, perturbation, state, toast } from '../core.js';
import { scatter } from '../chart.js';

let refs = {};
let ranUnder = null;
let rows = [];
let sort = { key: 'guaranteed', down: true };
let job = null;
let pick = null;
let selectedRow = null;
const PAGE = 400;

function targetOptions() {
  const model = state.model;
  if (!model) return [];
  const exchanges = model.reactions.filter((r) => r.boundary);
  return exchanges;
}

function fillTargets() {
  const options = targetOptions();
  const chosen = options.find((r) => r.id === 'EX_succ_e') ?? options[0];
  refs.target.replaceChildren(
    h('option', { value: '' }, 'None, growth only'),
    ...options.map((r) => h('option', { value: r.id, selected: r === chosen }, r.name || r.id)),
  );
}

function summary(reference) {
  const lethal = rows.filter((r) => r.growth < 0.05 * (reference || 1)).length;
  const reduced = rows.filter((r) => r.growth >= 0.05 * (reference || 1) && r.growth < 0.999 * (reference || 1)).length;
  const coupled = rows.filter((r) => (r.guaranteed ?? 0) > 1e-6).length;
  refs.summary.replaceChildren(
    tile(rows.length, 'Designs'),
    tile(lethal, 'Lethal'),
    tile(reduced, 'Slower growth'),
    tile(coupled, 'Coupled to the target'),
  );
}

function tile(value, label) {
  return h('div', { class: 'tile' }, h('strong', {}, Number(value).toLocaleString('en')), h('span', {}, label));
}

function sorted() {
  const { key, down } = sort;
  const value = (r) => (key === 'name' ? r.ids.join(' ') : r[key] ?? -Infinity);
  return [...rows].sort((a, b) => {
    const [x, y] = [value(a), value(b)];
    if (x === y) return b.growth - a.growth;
    return (x < y ? -1 : 1) * (down ? -1 : 1);
  });
}

function renderTable(reference) {
  const list = sorted().slice(0, PAGE);
  refs.body.replaceChildren(...list.map((r) => {
    const share = reference > 0 ? Math.max(0, Math.min(1, r.growth / reference)) : 0;
    const bar = h('div', { class: 'growth-bar' }, h('div', { style: { width: `${share * 100}%` } }));
    return h('tr', {
      class: r === selectedRow ? 'is-selected' : '',
      tabindex: 0,
      onclick: () => choose(r, reference),
      onkeydown: (e) => { if (e.key === 'Enter') applyToMap(r); },
    },
    h('td', { class: 'cell-ids' }, h('code', { class: 'ident' }, r.ids.join(' + '))),
    h('td', { class: 'num' }, fmt.flux(r.growth), bar),
    h('td', { class: 'num' }, fmt.flux(r.guaranteed)),
    h('td', { class: 'num' }, fmt.flux(r.possible)),
    h('td', { class: 'num' }, String(r.disabled)),
    h('td', {}, h('button', { class: 'link', type: 'button', onclick: (e) => { e.stopPropagation(); applyToMap(r); } }, 'Show on map')));
  }));
  refs.note.textContent = rows.length > PAGE ? `Showing the first ${PAGE} of ${rows.length.toLocaleString('en')}` : '';
}

function applyToMap(row) {
  const kind = refs.kind.value;
  emit('apply-design', {
    genes: kind === 'gene' ? [...new Set([...state.genes, ...row.ids])] : state.genes,
    reactions: kind === 'reaction' ? [...new Set([...state.reactions, ...row.ids])] : state.reactions,
  });
}

function choose(row, reference) {
  selectedRow = row;
  renderTable(reference);
  drawScatter(reference);
}

function drawScatter(reference) {
  const hasTarget = rows.some((r) => r.guaranteed !== null);
  // Values below the solver's tolerance are zero; left in, they would set the scale of the axis.
  const clean = (v) => (Math.abs(v ?? 0) < 1e-9 ? 0 : v);
  const points = rows.map((r) => ({
    x: r.growth,
    y: hasTarget ? clean(r.guaranteed) : r.disabled,
    r: 2.6,
    colour: (r.guaranteed ?? 0) > 1e-6 ? '#f4d64e' : r.growth < 0.05 * (reference || 1) ? '#ef5a6f' : '#6d97ff',
    alpha: 0.7,
  }));
  pick = scatter(refs.plot, {
    points,
    xlabel: 'Growth, 1/h',
    ylabel: hasTarget ? 'Guaranteed target flux at that growth' : 'Reactions disabled',
    mark: selectedRow ? { x: selectedRow.growth, y: hasTarget ? clean(selectedRow.guaranteed) : selectedRow.disabled } : null,
  });
}

async function run() {
  if (job !== null) {
    await invoke('cancel_job', { id: job });
    return;
  }
  const request = {
    kind: refs.kind.value,
    order: Number(refs.order.value),
    target: refs.target.value || null,
    min_growth: 0.05,
    ...perturbation(),
  };
  ranUnder = conditionsKey();
  refs.stale.hidden = true;
  refs.run.textContent = 'Stop';
  refs.run.classList.add('button-cut');
  refs.progress.hidden = false;
  rows = [];
  selectedRow = null;
  try {
    job = await invoke('start_screen', { request });
    const result = await follow(job, (status) => {
      const share = status.total ? status.done / status.total : 0;
      refs.bar.style.width = `${share * 100}%`;
      refs.progressText.textContent = status.total
        ? `${status.done.toLocaleString('en')} of ${status.total.toLocaleString('en')} designs, ${status.seconds.toFixed(1)} s`
        : 'Preparing designs';
    });
    rows = result.rows;
    refs.reference = result.reference;
    summary(result.reference);
    renderTable(result.reference);
    drawScatter(result.reference);
    refs.results.hidden = false;
    refs.empty.hidden = true;
  } catch (error) {
    toast(String(error));
  } finally {
    job = null;
    refs.run.textContent = 'Run screen';
    refs.run.classList.remove('button-cut');
    refs.progress.hidden = true;
  }
}

function header(label, key, numeric = true) {
  const th = h('th', { class: numeric ? 'num' : '', scope: 'col' }, h('button', {
    class: 'th-button',
    type: 'button',
    onclick: () => {
      sort = { key, down: sort.key === key ? !sort.down : numeric };
      renderTable(refs.reference);
    },
  }, label));
  return th;
}

export default {
  id: 'screen',
  label: 'Screen',
  icon: ['M4 5h16', 'M4 12h10', 'M4 19h6', 'M17 15l4 4', 'M18 12a3 3 0 1 1-6 0a3 3 0 1 1 6 0'],
  mount(el) {
    refs.kind = h('select', { class: 'select', 'aria-label': 'Delete' }, h('option', { value: 'gene' }, 'Genes'), h('option', { value: 'reaction' }, 'Reactions'));
    refs.order = h('select', { class: 'select', 'aria-label': 'How many at once' }, h('option', { value: 1 }, 'One at a time'), h('option', { value: 2 }, 'Every pair'));
    refs.target = h('select', { class: 'select wide', 'aria-label': 'Target product' });
    refs.run = h('button', { class: 'button button-primary', type: 'button', onclick: run }, 'Run screen');
    refs.bar = h('div', { class: 'progress-fill' });
    refs.progressText = h('span', {});
    refs.progress = h('div', { class: 'progress', hidden: true }, h('div', { class: 'progress-track' }, refs.bar), refs.progressText);
    refs.progress.hidden = true;
    refs.summary = h('div', { class: 'tiles' });
    refs.plot = h('canvas', { class: 'plot', 'aria-label': 'Growth against guaranteed target flux for each design' });
    refs.body = h('tbody');
    refs.note = h('p', { class: 'hint' });
    refs.results = h('div', { class: 'screen-results', hidden: true },
      refs.summary,
      h('div', { class: 'screen-split' },
        h('div', { class: 'panel plot-panel' }, refs.plot),
        h('div', { class: 'panel table-panel' },
          h('table', { class: 'table' },
            h('thead', {}, h('tr', {},
              header('Deleted', 'name', false), header('Growth', 'growth'), header('Guaranteed', 'guaranteed'),
              header('Possible', 'possible'), header('Reactions off', 'disabled'), h('th', {}))),
            refs.body),
          refs.note)));
    refs.results.hidden = true;
    refs.empty = h('div', { class: 'empty' }, h('h2', {}, 'Rank every deletion for a product'), h('p', {}, 'Choose what to delete and a target, then run the screen.'));
    refs.stale = h('span', { class: 'stale-note', hidden: true }, 'Deletions or uptake limits changed since this run');
    refs.stale.hidden = true;
    on('solution', () => { refs.stale.hidden = ranUnder === null || ranUnder === conditionsKey(); });
    el.append(
      h('header', { class: 'view-head' },
        h('div', { class: 'view-title' }, h('h1', {}, 'Knockout screen'), deletionNote(), refs.stale),
        h('div', { class: 'controls' },
          h('label', { class: 'field' }, 'Delete', refs.kind),
          h('label', { class: 'field' }, 'At once', refs.order),
          h('label', { class: 'field' }, 'Target', refs.target),
          refs.run)),
      refs.progress,
      refs.empty,
      refs.results,
    );
    refs.plot.addEventListener('click', (e) => {
      if (!pick) return;
      const box = refs.plot.getBoundingClientRect();
      const i = pick(e.clientX - box.left, e.clientY - box.top);
      if (i !== null) {
        selectedRow = rows[i];
        sort = { key: sort.key, down: sort.down };
        renderTable(refs.reference);
        drawScatter(refs.reference);
      }
    });
    new ResizeObserver(() => { if (rows.length) drawScatter(refs.reference); }).observe(refs.plot);
    on('model', () => {
      fillTargets();
      rows = [];
      refs.results.hidden = true;
      refs.empty.hidden = false;
    });
  },
};
