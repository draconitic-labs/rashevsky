// What a network can do: production envelopes and flux variability.

import { fmt, h, invoke, on, perturbation, state, toast } from '../core.js';
import { lineChart } from '../chart.js';

let refs = {};
let last = null;
let lastRanges = [];

function fill() {
  const options = state.model?.reactions ?? [];
  const pick = (node, preferred) => {
    const chosen = options.find((r) => r.id === preferred) ?? options[0];
    node.replaceChildren(...options.map((r) => h('option', { value: r.id, selected: r === chosen }, r.name ? `${r.name} (${r.id})` : r.id)));
  };
  pick(refs.control, state.model?.objective[0]);
  pick(refs.target, 'EX_ac_e');
}

async function runEnvelope() {
  try {
    const points = await invoke('envelope', { control: refs.control.value, target: refs.target.value, points: 40, perturbation: perturbation() });
    last = points;
    draw();
  } catch (error) {
    toast(String(error));
  }
}

function draw() {
  if (!last) return;
  const x = last.map((p) => p.x);
  const low = last.map((p) => p.minimum);
  const high = last.map((p) => p.maximum);
  lineChart(refs.plot, {
    bands: [{ x, low, high, colour: '#f4d64e' }],
    series: [{ x, y: low, colour: '#6d97ff', width: 2 }, { x, y: high, colour: '#f4d64e', width: 2 }],
    xlabel: optionText(refs.control),
    ylabel: optionText(refs.target),
  });
}

function optionText(select) {
  return select.options[select.selectedIndex]?.textContent ?? '';
}

async function runRanges() {
  try {
    const fraction = Number(refs.fraction.value);
    refs.fractionOut.textContent = `${Math.round(fraction * 100)} %`;
    lastRanges = await invoke('variability', { fraction, perturbation: perturbation() });
    renderRanges();
  } catch (error) {
    toast(String(error));
  }
}

function renderRanges() {
  const query = refs.filter.value.trim().toLowerCase();
  const scale = Math.max(1e-9, ...lastRanges.map((r) => Math.max(Math.abs(r.minimum), Math.abs(r.maximum))).filter(Number.isFinite));
  const rows = lastRanges
    .filter((r) => !query || r.id.toLowerCase().includes(query))
    .sort((a, b) => (b.maximum - b.minimum) - (a.maximum - a.minimum))
    .slice(0, 300);
  refs.rows.replaceChildren(...rows.map((r) => {
    const left = 50 + (50 * r.minimum) / scale;
    const width = Math.max(0.4, (50 * (r.maximum - r.minimum)) / scale);
    return h('tr', {},
      h('td', {}, h('code', { class: 'ident' }, r.id)),
      h('td', { class: 'num' }, fmt.flux(r.minimum)),
      h('td', { class: 'num' }, fmt.flux(r.maximum)),
      h('td', { class: 'range-cell' }, h('div', { class: 'range-track' }, h('div', { class: 'range-zero' }), h('div', { class: 'range-span', style: { left: `${left}%`, width: `${width}%` } }))));
  }));
}

export default {
  id: 'ranges',
  label: 'Ranges',
  icon: ['M4 20V4', 'M4 20h16', 'M8 16c3-1 5-4 6-9', 'M8 20c5-1 9-5 11-13'],
  mount(el) {
    refs.control = h('select', { class: 'select wide', 'aria-label': 'Flux along the horizontal axis', onchange: runEnvelope });
    refs.target = h('select', { class: 'select wide', 'aria-label': 'Flux along the vertical axis', onchange: runEnvelope });
    refs.plot = h('canvas', { class: 'plot', 'aria-label': 'Production envelope' });
    refs.fraction = h('input', { type: 'range', min: 0, max: 1, step: 0.05, value: 1, 'aria-label': 'Share of the optimal objective kept', onchange: runRanges, oninput: (e) => { refs.fractionOut.textContent = `${Math.round(e.target.value * 100)} %`; } });
    refs.fractionOut = h('output', {}, '100 %');
    refs.filter = h('input', { class: 'find', type: 'search', placeholder: 'Filter reactions', 'aria-label': 'Filter reactions', oninput: renderRanges });
    refs.rows = h('tbody');
    el.append(
      h('header', { class: 'view-head' },
        h('h1', {}, 'Ranges'),
        h('div', { class: 'controls' },
          h('label', { class: 'field' }, 'Horizontal', refs.control),
          h('label', { class: 'field' }, 'Vertical', refs.target))),
      h('div', { class: 'ranges-split' },
        h('div', { class: 'panel plot-panel' }, refs.plot),
        h('div', { class: 'panel table-panel' },
          h('div', { class: 'panel-head' },
            h('h2', {}, 'Flux variability'),
            h('label', { class: 'field field-inline' }, 'Objective kept', refs.fraction, refs.fractionOut),
            refs.filter),
          h('div', { class: 'scroll' }, h('table', { class: 'table' },
            h('thead', {}, h('tr', {}, h('th', { scope: 'col' }, 'Reaction'), h('th', { class: 'num', scope: 'col' }, 'Minimum'), h('th', { class: 'num', scope: 'col' }, 'Maximum'), h('th', { scope: 'col' }, 'Range'))),
            refs.rows)))),
    );
    new ResizeObserver(draw).observe(refs.plot);
    on('model', () => { fill(); last = null; lastRanges = []; refs.rows.replaceChildren(); });
    refs.active = false;
    on('solution', () => { if (refs.active) { runEnvelope(); runRanges(); } });
  },
  show() {
    refs.active = true;
    runEnvelope();
    runRanges();
  },
  hide() {
    refs.active = false;
  },
};
