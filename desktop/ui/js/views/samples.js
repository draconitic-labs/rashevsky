// Uniform samples of the feasible flux space, and whether the chains have mixed.

import { conditionsKey, deletionNote, fmt, follow, h, invoke, on, perturbation, ramp, state, toast } from '../core.js';
import { bars } from '../chart.js';

let refs = {};
let ranUnder = null;
let result = null;
let busy = false;

async function run() {
  if (busy) return;
  busy = true;
  refs.run.disabled = true;
  ranUnder = conditionsKey();
  refs.stale.hidden = true;
  refs.status.textContent = 'Sampling';
  try {
    const keep = refs.keep.checked ? Number(refs.fraction.value) : null;
    const id = await invoke('start_sampling', {
      request: {
        method: refs.method.value,
        chains: 4,
        samples: Number(refs.samples.value),
        thinning: Number(refs.thin.value),
        burn_in: 20000,
        seed: 1,
        fraction: keep,
        ...perturbation(),
      },
    });
    result = await follow(id, (s) => { refs.status.textContent = `Sampling, ${s.seconds.toFixed(1)} s`; });
    refs.status.textContent = '';
    refs.empty.hidden = true;
    refs.gridHead.hidden = false;
    renderTiles();
    renderGrid();
  } catch (error) {
    refs.status.textContent = '';
    toast(String(error));
  } finally {
    busy = false;
    refs.run.disabled = false;
  }
}

function renderTiles() {
  const d = result.diagnostics;
  const tile = (value, label, good) => h('div', { class: `tile ${good === false ? 'tile-warn' : ''}` }, h('strong', {}, value), h('span', {}, label));
  refs.tiles.replaceChildren(
    tile(String(result.dimension), 'Free dimensions'),
    tile(result.total_samples.toLocaleString('en'), 'Samples'),
    tile(d.max_split_rhat.toFixed(3), 'Largest split R-hat', d.max_split_rhat < 1.05),
    tile(Math.round(d.min_ess).toLocaleString('en'), 'Smallest effective size', d.min_ess > 100),
    tile(d.max_geweke_z.toFixed(1), 'Largest Geweke z', d.max_geweke_z < 5),
  );
}

function renderGrid() {
  if (!result) return;
  const query = refs.filter.value.trim().toLowerCase();
  const order = result.reactions
    .map((id, i) => ({ id, i, sd: result.diagnostics.sd[i] }))
    .filter((x) => x.sd > 1e-9 && (!query || x.id.toLowerCase().includes(query)))
    .sort((a, b) => b.sd - a.sd)
    .slice(0, 60);
  const top = Math.max(1e-9, ...order.map((x) => x.sd));
  refs.grid.replaceChildren(...order.map((x) => {
    const canvas = h('canvas', { class: 'hist', 'aria-hidden': 'true' });
    const hist = result.histograms[x.i];
    const card = h('div', { class: 'hist-card' },
      h('div', { class: 'hist-top' }, h('code', { class: 'ident' }, x.id), h('span', {}, `${fmt.flux(result.diagnostics.mean[x.i])} ± ${fmt.flux(x.sd)}`)),
      canvas,
      h('div', { class: 'hist-axis' }, h('span', {}, fmt.flux(hist.low)), h('span', {}, fmt.flux(hist.high))));
    requestAnimationFrame(() => bars(canvas, hist.counts, ramp(Math.log1p(x.sd) / Math.log1p(top))));
    return card;
  }));
}

export default {
  id: 'samples',
  label: 'Samples',
  icon: ['M4 20h16', 'M6 20v-4', 'M10 20V9', 'M14 20V5', 'M18 20v-8'],
  mount(el) {
    refs.method = h('select', { class: 'select', 'aria-label': 'Sampler' }, h('option', { value: 'chrr' }, 'CHRR'), h('option', { value: 'optgp' }, 'OptGP'));
    refs.samples = h('select', { class: 'select', 'aria-label': 'Samples per chain' }, ...[500, 1000, 2500, 5000].map((n) => h('option', { value: n, selected: n === 1000 }, n.toLocaleString('en'))));
    refs.thin = h('select', { class: 'select', 'aria-label': 'Moves between samples' }, ...[50, 100, 500, 1000].map((n) => h('option', { value: n, selected: n === 500 }, n.toLocaleString('en'))));
    refs.keep = h('input', { type: 'checkbox', checked: true, 'aria-label': 'Keep part of the optimal objective' });
    refs.fraction = h('input', { class: 'number', type: 'number', min: 0, max: 1, step: 0.05, value: 0.9, 'aria-label': 'Share of the optimal objective' });
    refs.run = h('button', { class: 'button button-primary', type: 'button', onclick: run }, 'Sample');
    refs.status = h('span', { class: 'hint' });
    refs.tiles = h('div', { class: 'tiles' });
    refs.filter = h('input', { class: 'find', type: 'search', placeholder: 'Filter reactions', 'aria-label': 'Filter reactions', oninput: renderGrid });
    refs.grid = h('div', { class: 'hist-grid' });
    refs.gridHead = h('div', { class: 'grid-head', hidden: true }, refs.filter);
    refs.gridHead.hidden = true;
    refs.empty = h('div', { class: 'empty' }, h('h2', {}, 'Sample what the network can do'), h('p', {}, 'Each reaction gets a distribution of fluxes over the feasible space, with convergence checks.'));
    refs.stale = h('span', { class: 'stale-note', hidden: true }, 'Deletions or uptake limits changed since this run');
    refs.stale.hidden = true;
    on('solution', () => { refs.stale.hidden = ranUnder === null || ranUnder === conditionsKey(); });
    el.append(
      h('header', { class: 'view-head' },
        h('div', { class: 'view-title' }, h('h1', {}, 'Flux samples'), deletionNote(), refs.stale),
        h('div', { class: 'controls' },
          h('label', { class: 'field' }, 'Sampler', refs.method),
          h('label', { class: 'field' }, 'Per chain', refs.samples),
          h('label', { class: 'field' }, 'Moves between', refs.thin),
          h('label', { class: 'field field-check' }, h('span', {}, 'Keep objective'), h('span', { class: 'inline' }, refs.keep, refs.fraction)),
          refs.run, refs.status)),
      refs.tiles,
      refs.empty,
      refs.gridHead,
      h('div', { class: 'scroll' }, refs.grid),
    );
    on('model', () => { result = null; refs.tiles.replaceChildren(); refs.grid.replaceChildren(); refs.empty.hidden = false; refs.gridHead.hidden = true; });
  },
};
